function waveform =myCos(frequency, phase, duration, Fs)
% function takes in frequency, phase shift, duration, and 
% discrete sampling window to construct a cosine. 
% Inputs:
%   frequency = cosine frequency (Hz)
%   phase = cosine phase shift (radians)
%   duration = total signal time (s)
%   Fs = sampling rate (Hz)
% Outputs:
%   waveform = cosine at frequency and phase sampled at Fs for duration
samplingInterval = 1/Fs;
numSamples = mod(duration/samplingInterval, 0);
sampleTime = linspace(0, numSamples*samplingInterval, numSamples+1);
%return sampled cosine waveform
waveform = cos(frequency*2*pi*sampleTime + phase);
