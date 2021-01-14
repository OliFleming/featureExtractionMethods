%
% APPROXEQUAL	Logic function that checks if two numbers are equal
%    within a specified error.
%
% eq = approxequal(x,y,err)
%
% Author Adrian Chan
%
% Checks if y - err <= x <= y + err
%
% Inputs
%    x		first number to compare
%    y		second number to compare
%    err	accuracy of comparison (default 1e-10)
%
% Outputs
%    eq	returns 1 if the two numbers are approximately equal
%			returns 0 if the two numbers are not approximately equal
%
% Modifications
% 01/05/20 AC First created.

function eq = approxequal(x,y,err)

if nargin < 3
   err = 1e-10;
end

eq = ( (x >= (y-err)) & (x <= (y+err)) );
