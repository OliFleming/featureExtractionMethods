%
% FIND_DELAY Estimates delay between two signals
%
% d = find_delay(x,y,maxdelay)
%
% Author Adrian Chan
%
% Find the shift between x and y by performing a cross-correlation
% to find T where y = x(t - T) and t and T are units are the
% array indices.
%
% Inputs
%    x: signal 1
%    y: signal 2
%    maxdelay: maximum abs delay between the two signals
%                default: T/2 (may not work for larger maxdelay)
%
% Outputs
%    d: delay between x and y
%
% Modifications
% 04/08/09 AC Added maxdelay input
% 03/02/21 AC First created.
function d = find_delay(x,y,maxdelay)

if length(x) > length(y)
    x = x(1:length(y));
elseif length(y) > length(x)
    y = y(1:length(x));
end

if nargin < 3
    maxdelay = length(x)/2;
end

z = xcorr(x,y,'unbiased');
t = -(length(x)-1):(length(x)-1);
mask = (abs(t) < maxdelay);
indx = find(z == max(z(mask)));
d = min(t(indx));
