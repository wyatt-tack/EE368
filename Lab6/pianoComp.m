function [myWave] = pianoComp()
% Function runs and compiles estimated values into piano key based on
% piano.wav file. 
[piano, Fpiano] = audioread("piano.wav");
% Following coefficients estimated from the function below
[estMag, estPhase] = DFS(piano, Fpiano);
Fnot = 260.94
estFreq = linspace(Fnot, 10*Fnot, 10);
% initialize myWav
myWav = linspace(0,0,44.1e3);
for n = 1:length(estFreq)
    myWav = myWav + myCos(estFreq(n), estPhase(n), 1, (44.1e3) - 1);
end    
%sound(piano, 44.1e3)
sound(myWav, 44.1e3)
end