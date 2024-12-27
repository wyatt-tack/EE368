function [S,DATA,PLAYER] = livespectrogram(path,varargin)
% LIVESPECTROGRAM
% K. Clay McKell
% v 2.0 2020-08-07
% SYNTAX
% livespectrogram(path)     'path' is a string containing the path to the
% audio file.
% livespectrogram(...,'fig',NFIG)   Paints to figure NFIG.  If figure NFIG
% does not exist, it is created.  If NFIG is set to 0, then plotting is
% suppressed.
% S = livespectrogram(...)  Returns a cell array of length T containing FFT 
% values for nonnegative frequency where T indexes time.
% [S,DATA] = livespectrogram(...)   Returns the time-domain audio data in
% DATA.
% [S,DATA,PLAYER] = livespectrogram(...)    Returns the audio player object
% in PLAYER.

%% Input checking
yesplot = true; nfig = -1;
if (nargin == 2*ceil(nargin/2)) || (nargin > 3)
    error('livespectrogram:nargin','Improper number of input arguments');
elseif nargin > 1
    if strcmp(varargin{1},'fig')
        nfig = varargin{2};
        yesplot = nfig>0;
    else
        error('livespectrogram:name','Invalid input name.');
    end
end    
%% Actual good stuff
imported_data = importdata(path);
DATA = imported_data.data;
fs = imported_data.fs;
Ts = 1/fs;
frame_min = 1/24;
c = 1;
S = cell(1);
if yesplot && nfig>0
    fig = figure(nfig);
elseif yesplot
    fig = figure;
end
if yesplot
    plot3(0,0,0); hold on; view(98,29);
    xlabel('Time [s]'); ylabel('Frequency [Hz]'); zlabel('Power');
end
PLAYER = audioplayer(DATA,fs);
PLAYER.play;
pause(frame_min);
while PLAYER.isplaying
    tic;
    C = PLAYER.CurrentSample;
    L = C-c+1;
    s = fft(DATA(c:C)/L);
    s = s(1:ceil(L/2));
    if yesplot
        figure(fig);
        freqs = fs*(0:ceil(L/2)-1)/L;
        plot3(C*Ts*ones(size(freqs)),freqs,abs(s));
    end
    if isempty(S{1})
        S{1} = s;
    else
        S{end+1} = s;
    end
    c = C+1;
    t = toc;
    % This microsecond pause is necessary for the loop to terminate:
    % https://www.mathworks.com/matlabcentral/answers/46603-audioplayer-isplaying-won-t-exit-tight-loop
    pause(max(1e-6,frame_min-t));
end
if yesplot
    hold off;
end