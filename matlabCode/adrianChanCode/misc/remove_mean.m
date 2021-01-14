%
% REMOVE_MEAN Removes the mean.
%
% y = remove_mean(x)
%
% Author Adrian Chan
%
% This function removes the mean of the columns.
%
% Inputs
%    x: columns of signals
%
% Outputs
%    y: zero mean signals
%
% Modifications
% 04/08/09 AC First created.
function y = remove_mean(x)

N = size(x,1);

y = x - repmat(mean(x),N,1);
