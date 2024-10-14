function plspow = pulsepower(x,tr)

% Function to take a DIMAC timeseries and calculate the ratio of average
% power over cardiac pulse frequencies (40-120 bpm) to total power averaged
% over all frequencies
%
% IDD 20/07/2023
%
% USAGE:
%       plspow  = pulsepower(x,tr);
%
%       x       = DIMAC timeseries
%       tr      = temporal resolution (in seconds)

fs = 1/tr; % sampling frequency = 1/repetition time
NyqFreq = fs/2; % Nyquist frequency

% Low-pass filter required to remove DC and low-order offsets that could
% exaggerate total power and underestimate the pulse power
X=fourier_design_matrix(numel(x),7,0);
B=regress(double(x),X);
x_lp=X*B;

Fx=fftshift(fft(double(x)-x_lp)); % FFT of the detrended timeseries
kx = linspace(-NyqFreq,+NyqFreq,numel(x)); % Frequency vector (-/+ Nyquist)

% Pulse frequency band defined in the range 40-120 beats per minute
% (0.66-2 Hz)
pulsePower = mean(abs(Fx(((abs(kx)>40/60).*(abs(kx)<2))==1)).^2);
totPower = mean(abs(Fx).^2); % Total power averaged over all frequencies

plspow = pulsePower./totPower;