clear all
N = 9;                       % Order
USAMPR=8;                       %Intended Upsampling rate
F = [[0 0.2 0.4 0.6 0.8]/USAMPR 1];    % Frequency Vector
A = [1 1 1 0.4 0.6 0.8];  % Amplitude Vector
W = [1 1 1];                  % Weight Vector

% Calculate the coefficients using the FIRLS function.
b  = firls(N, F, A, W);
Hd = dfilt.dffir(b);
phase_to_mag=20;
channel=Hd.Numerator*((phase_to_mag-1)/phase_to_mag)-j*Hd.Numerator./phase_to_mag
fvtool(channel)