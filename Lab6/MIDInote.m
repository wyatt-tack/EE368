function myWav = MIDInote(MIDI, T,Fs)
% Function accepts MIDI musical note number, Time for note to be played,
% and frequency for note to be sampled/played at. piano coefficients 
% for magnitudes and phase compiled from DFS function and piano.wav file
estMag = [1.0000, 0.6697, 0.2179, 0.1188, 0.1590, 0.1589, 0.0112, 0.0054, 0.0030, 0.0005];
estPhase = [-0.3902, -0.1212, 2.7223, 1.0472, 0.0450, 1.4863, 0.8726, 1.5655, -2.7004, 2.9044];
Fnot = (2^((MIDI-69)/12))*(440);
% initialize myWav
myWav = linspace(0,0,T*Fs);
for n = 1:10
    myWav = myWav + estMag(n)*myCos(Fnot*n, estPhase(n), T, (Fs) - 1);
end    
sound(myWav, Fs)
end