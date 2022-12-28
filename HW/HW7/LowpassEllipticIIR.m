function Hd = LowpassEllipticIIR
%LOWPASSELLIPTICIIR Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.9 and Signal Processing Toolbox 8.5.
% Generated on: 15-Apr-2022 08:26:19

% Elliptic Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

Fpass = 250;     % Passband Frequency
Fstop = 300;     % Stopband Frequency
Apass = 1;       % Passband Ripple (dB)
Astop = 50;      % Stopband Attenuation (dB)
match = 'both';  % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'ellip', 'MatchExactly', match);

% [EOF]