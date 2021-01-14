% MEANFREQUENCY This computes the mean frequency of a power spectrum.
%
% mf = meanfrequency(f,p)
%
% Author Adrian Chan
%
% This computes the mean fequency of a power spectrum.
%
% Reference: Öberg T, Sandsjö L, Kadefors R, "EMG mean power frequency: 
% Obtaining a reference value", Clinical Biomechanics, vol. 9, pp. 253-257,
% 1994.
%
% Inputs
%    f: frequencies (Hz)
%    p: power spectral density values
%
% Outputs
%    mf: mean frequency
%
% Modifications
% 09/09/21 AC First created.
function mf = meanfrequency(f,p)

mf = sum(f.*p)/sum(p);
