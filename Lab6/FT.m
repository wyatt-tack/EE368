function [DFTx,fcDFT,FdDFT]=FT(x_samples,Fs,noPlot)
%function [DFTx,fcDFT,FdDFT]=FT(x_samples,Fs,noPlot)
% Function computes the discrete Fourier Transform and an approximation of the
% continuous-time Fourier Transform and DTFT of a sampled signal (by zero-padding)
% using an FFT.  Returns equivalent continuous time and digital frequency
% axes for plotting, with the arrays limited to the Nyquist frequency and
% axes covering the principal range -Fs/2 -> +Fs/2.
% Inputs:
%   x_samples = samples of the time-domain signal
%   Fs = the sampling rate (in Hz)
%   noPlot = any value will stop the function from plotting results
% Outputs:
%   DFTx = Sampled DFT spectrum values (complex) 
%           [same # samples as x_samples]
%   fcDFT = Analog frequency values for each DFT sample
%   FdDFT = Digital frequency values for each DFT sample
% Plots:
%   Figure 100: DFT Stem Plot with equivalent analog frequency axis
% THESE SUPPRESSED TO REDUCE STUDENT CONFUSION:
%   Figure 200: DFT Stem Plot with digital frequency axis
%   Figure 300: Approximated C-T FT Plot with equivalent analog frequencies
%   Figure 400: Approx. DTFT Plot with digital frequency axis

% Handle missing parameters
if nargin<3,  % If showing plots or not is not specified
    noPlot= false; % Do plots (Not told NOT to plot results)
else
     noPlot= true; % Any input will turn OFF plotting
end;
if nargin<2,  % If not Fs specified
    disp('No sampling rate provided.  Analog frequencies will be wrong.');
    Fs=1; % Make the analog freq axis match the digital freq 
end;


% Set sizes for the number of samples used in each spectrum.
% Discrete DFT spectrum uses same number of samples as input signal
num_discrete_spectrum_samples=length(x_samples);
% Continuous spectra (DTFT, C-TFT) uses more samples for spectral interp.
% Use power of 2 for fast FFT, but at least 256 samples for continuous FT's
num_cont_spectrum_samples=max([2.^(floor(log2(length(x_samples)))+1),256]);

%Compute the discrete spectrum
% fftshifted to show -1/2 -> + 1/2 principal range of freq.
% with normalized values
DFTx=fftshift(fft(x_samples,num_discrete_spectrum_samples))./length(x_samples);
% Set frequency axis to handle odd or even number of samples
FdDFT=(([-floor(num_discrete_spectrum_samples/2):...
    -floor(num_discrete_spectrum_samples/2)+num_discrete_spectrum_samples-1]...
    ./num_discrete_spectrum_samples));
fcDFT=FdDFT.*Fs;

% Compute an approximation of the C-TFT and DTFT
CFTx=fftshift(fft(x_samples,num_cont_spectrum_samples))./length(x_samples);
% Make analog frequency axis for this number of samples
fcCFT=(([-floor(num_cont_spectrum_samples/2):...
    -floor(num_cont_spectrum_samples/2)+num_cont_spectrum_samples-1]...
    ./num_cont_spectrum_samples)).*Fs;
% COnvert analog frequencies to digital frequencies for DTFT
FdCFT=fcCFT./Fs;

if noPlot==false,
figure(100)
hs=stem(fcDFT,abs(DFTx),'.');
%set(hs,'Marker','.')
grid on
xlabel('Equivalent Analog Frequency (Hz)')
xlim([-Fs/2 Fs/2])
ylabel('Magnitude Response (Normalized)')
title('Discrete Fourier Transform','FontWeight','bold')
end;
return;

if noPlot==false,
figure(200)
hs=stem(FdDFT,abs(DFTx),'.');
%set(hs,'Marker','.')
grid on
xlim([-.5 .5])
xlabel('Digital Frequency (cycles/sample)')
ylabel('Magnitude Response (Normalized)')
title('Discrete Fourier Transform','FontWeight','bold')

figure(300)
plot(fcCFT,abs(CFTx))
grid on
xlabel('Analog Frequency (Hz)')
xlim([-Fs/2 Fs/2])
ylabel('Magnitude Response')
title('Continuous-Time Fourier Transform (Approx)','FontWeight','bold')

figure(400)
plot(FdCFT,abs(CFTx));
grid on
xlabel('Digital Frequency (cycles/sample)')
xlim([-.5 .5]);
ylabel('Magnitude Response (Normalized)')
title('Discrete-Time Fourier Transform','FontWeight','bold')

figure(100)
end;
return;