%
% prd Computes percentage root mean square (rms) / residual difference (PRD).
%
% error = prd(ideal, approx, submean)
%
% Author: Vesal Badee
%
% Computes percentage root mean square (rms) / residual difference (PRD).
% Independent of the amplitude scale of the inputs.
%
% Inputs
%   ideal: column vector containing the ideal signal. Represents the
%           original signal.
%   approx: column vector containing the approximate signal. Represents
%           the observed outcome.
%   submean: 1 to subtract mean during normalization (default); otherwise
%           normalize without subtracting the mean.
%   Note: Rows represent data. Columns represent channels.
%
% Outputs
%   error: Scalar value representing the PRD between ideal and approx.
%
% Reference
%   Y. Zigel, A. Cohen & A. Katz, 'The weighted diagnostic distortion (WDD)
%   measure for ECG signal compression', IEEE Trans. on Biomedical
%   Engineering, vol. 47, no. 11, Nov. 2000, Equations (1) & (2)
%
% Modifications
% August 7, 2004 VB First created.
% 05/09/27 AC Changed to enable multiple columns for the approx input
%
% Version 0.1

function error = prd(ideal, approx, submean)

if (nargin < 3)
    submean = 1;
end

Nsignal = size(approx,2);
Ndata = size(ideal,1);
ideal = repmat(ideal,1,Nsignal);

% Compute numerator (residuals).
num = sum((ideal - approx).^2,1);

% Compute denominator (used for normalization).
if (submean == 1)
    den = sum((ideal - repmat(mean(ideal,1),Ndata,1)).^2,1);    % Equation (2)
else
    den = sum(ideal.^2,1);    % Equation (1)
end

% Compute error.
error = sqrt(num./den) * 100;