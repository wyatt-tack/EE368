function y_samples = digital_filter(input_wavfile,output_wavfile)
%function yout = digital_filter(input_wavfile,output_wavfile)
%   Inplements a digital filter and uses it to process the contents of a
%   .wav file

% Global Variables
global x_samples;
global y_samples;
global sample_index; sample_index=1;

% Read in the input waveform and get needed info
[x_samples,Fs]=audioread(input_wavfile);
% If input waveform is stereo, make it mono
if (size(x_samples,2)>1), x_samples = x_samples(:,1);end;
number_of_samples = length(x_samples)
y_samples = zeros(size(x_samples));  % Initialize Output array



%***********************************************************************
% STUDENTS DEFINE FILTER PARAMETERS HERE:
% Create and initialize variables or arrays for your filter coefficients
%  and previous Inputs or Outputs needed for you difference equation.

B0 = 0.0354;% * 0.0245;
B1 = -0.0000;% * 0.955;
B2 = -0.0455;% * 0.2061;
B3 = 0.0000;% * 0.3455;
B4 = 0.0637;% * 0.5000;
B5 = -0.0000;% * 0.6545;
B6 = -0.1061;% * 0.7939;
B7 = 0.0000;% * 0.9045;
B8 = 0.3183;% *  0.9755;
B9 = 0.5000;% * 1.000;
B10 = 0.3183;% * 0.9755;
B11 = 0.0000;% * 0.9045;
B12 = -0.1061;% * 0.7939;
B13 = -0.0000;% * 0.6545;
B14 = 0.0637;% * 0.5000;
B15 = 0.0000;% * 0.3455;
B16 = -0.0455;% * 0.2061;
B17 = -0.0000;% * 0.0955;
B18 = 0.0354;% *0.0245;

xnm1=0;
xnm2=0;
xnm3=0;
xnm4=0;
xnm5=0;
xnm6=0;
xnm7=0;
xnm8=0;
xnm9=0;
xnm10=0;
xnm11=0;
xnm12=0;
xnm13=0;
xnm14=0;
xnm15=0;
xnm16=0;
xnm17=0;
xnm18=0;
%*********************************************************************
%   (End of student variable definitions)



% Loop to run through all the samples in the input file, 
%  one sample at a time; pretending to be an A/D sample timer interrupt:
for sample_num = 1:number_of_samples
    
    % New sample arrived for processing...
    %  (this would usually be a sample timer interrupt service routine):
    
    %--------------------------------------------------------------
    % STUDENTS ENTER YOUR DIFFERENCE EQUATION COMPUTATION AND ANY OTHER
    % OPERATIONS NEEDED FOR COMPUTING EACH OUTPUT SAMPLE VALUE
    % AND PREPARING FOR THE NEXT SAMPLE TO ARRIVE
    
    xn = Read_ADC(); % Get a new input sample
    yn = (B0.*xn) + (B1.* xnm1) + (B2.*xnm2) + (B3.*xnm3) + (B4.* xnm4) + (B5.*xnm5) + (B6.* xnm6);
    xnm6 = xnm5;
    xnm5 = xnm4;
    xnm4 = xnm3;
    xnm3 = xnm2;
    xnm2 = xnm1;
    xnm1=xn;
    Write_DAC(yn); % Output your result
    
    %----END OF STUDENT PROGRAM AREA--------------------------------
    %  (Do not change anything below this line)
    
    
end  % End of loop for getting each sample

% After all samples are processed, write the results to the output file
audiowrite(output_wavfile, y_samples, Fs);

% Look at the spectrum of the output samples
Analyze_Data(y_samples, Fs,100);
end   % End of the Main Program


% SUPPORTING FUNCTIONS

function xn = Read_ADC()
global sample_index;
global x_samples;
% Get next input sample to process
xn=x_samples(sample_index);
sample_index = sample_index + 1;
end

function Write_DAC(yn)
global sample_index;
global y_samples;
y_samples(sample_index-1) = yn;
end

function Analyze_Data(samples,fs, fig_num)
npts = length(samples);
figure(fig_num)
subplot(2,2,1)
%spectrogram(samples,hamming(256),0,256,fs);
spectrogram(samples,4096,4096-256,4096,fs);
colormap gray
[cmin, cmax]=caxis;
caxis([cmax- 60, cmax]);
title('Power Spectral Density Spectrogram')

subplot(2,2,2)
plot([0:npts-1]/npts*fs,abs(fft(samples))/npts)
xlim([0,fs/2])
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of All Data')
grid on;

subplot(2,2,3)
% Plot an averaged magnitude spectrum using sqrt of power spectral density
% by Welch's method (for noisy data spectra)
%[Pxx,F] = pwelch(samples,hamming(256),128,256,fs);
[Pxx,F] = pwelch(samples,4096,4096-256,4096,fs);
plot(F,sqrt(Pxx));
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Estimated Average Spectrum of All Data')
grid on;

subplot(2,2,4)
% Show a spectrum for every 0.1 sec of data
spectrum_time_interval = 0.1; % seconds for each spectrum shown
samples_per_spectrum_time = round(spectrum_time_interval * fs);
block_size=2.^nextpow2(samples_per_spectrum_time); % round up block size to 2^N for faster FFT
sp=spectrogram(samples,hamming(block_size),0,block_size,fs);
[spec_length,num_spectra]=size(sp);
peak_mag = max(max(abs(sp)))/spec_length;
for ii=1:num_spectra
    plot([0:spec_length-1]/spec_length*fs/2,abs(sp(:,ii))/spec_length)
    ylim([0,peak_mag])
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title(['Block Spectrum ', int2str(ii), ' of ',int2str(num_spectra)])
    grid on;
    
    pause(0.9*spectrum_time_interval);
end

end



