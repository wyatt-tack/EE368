function [sampleTime,sinusoid,ftSin, dftSin]=sinAnalyze(N_samples,Ts,f0, time_domain_plot, freq_domain_plot, soundflag)
% function takes in number of samples, sampling period, and sinusoidal
% frequency to sample a sine at that frequency for the given samples. 
% Inputs:
%   N_samples = samples of the time-domain signal
%   Ts = the sampling period (in S)
%   f0 = frequency of 2Vpp sine wave to be analyzed  
%   time_domain_plot = if any argument sin will be discretely plotted
%   freq_domain_plot = if any argument, dft will be plotted
%   sound = if any argument, sine will be played out of speakers
% Outputs:
%   sampleTime = returned array of timestamps for samples at period Ts
%   sampleSin = returned sineusoid array of sampled points
%   ftSin = returns estimated continuous time fourier series of sinusoid
%   dftSin = returns dft of sampled sinusoid
% Plots:
%   Figure 200: Continous time stem plot of sampled sinusoid
%   Figure 100: DFT Stem Plot with equivalent analog frequency axis
noTDP = true; %set automatic reset for plots to be off
noFDP = true;
noSound = true;
switch nargin  %argument case by case, if too few arguments display, parse which plots wanted
    case 1
        disp('No sampling rate provided');
    case 2
        disp('No sinusoidal frequency');
    case 4
        noTDP = false;
    case 5
        noTDP = false;
        noFDP = false;
    case 6
        noTDP = false;
        noFDP = false;
        noSound = false;
    otherwise 
        disp('Too many arguments provided');
end
% use sample number and sampling period to get time sample array

for sampleIndex = 0:N_samples
    sampleTime(sampleIndex+1) = sampleIndex * Ts;
end
sinusoid = sin(2*pi*f0.*sampleTime); %use time sample array to sample sinusoid at f0
%FT function courtosy of Prof. Pilkington
[dftSin, ftSin] = FT(sinusoid, 1/Ts, 'no plot'); %use FT to recieve DFt, Cft
% plot time domain sinusoid if requested
if noTDP == false,
    figure(300)
    newplot;
    stem(sampleTime,sinusoid);
    xlabel('Time (s)');
    ylabel('x(t) = sin(2pi*f0*n*Ts)');
    title('Discrete f0 Sinusoid Sampled at Ts','FontWeight','bold');
end;    
% plot frequency domain sinusoid if requested, taken from FT function
% courtosy  of Prof. Pilkington
if noFDP == false,
    FT(sinusoid, 1/Ts);
end;    
% play sound if sound requested
if noSound == false,
    sound(sinusoid, 1/Ts);
end;