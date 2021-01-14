%
% LOGGAUSSPDF Computes the probability density function for a Gaussian distribution.
%
% p = loggausspdf(x,mu,C)
%
% Author Adrian Chan
%
% Computes the PDF using the formula
% f(x) = 1/sqrt((2*pi)^n*det(C)) * exp(-1/2 * X * inv(C) * X')
% see page 197 equation 8-58
% Papooulis A, Probability, random variables, and stochastic processes, 3rd ed, McGraw-Hill New York NY, 1991.
% but returns the natural log of this answer
%
% Inputs
%    x: pdf to be evaluated here (column vector, multiple column can be used to indicate multiple inputs)
%    mu: mean vector of the Gaussian distribution (default = 0)
%    C: covariance matrix of the Gaussian distribution (default = identity matrix)
%
% Outputs
%    p: pdf evaluated at x (columns indicate multiple outputs)
%
% Modifications
% 03/04/08 AC First created.

function p=gausspdf(x, mu, C)

n = size(x,1); % vector size
M = size(x,2); % number of inputs

if nargin < 3
    C = eye(n);
    if nargin < 2;
        mu = ones(n,1);
    end
end

x = x - mu*ones(1,M);

logdenom = -log((2*pi)^(n/2)*sqrt(abs(det(C))));
invC = inv(C);
mahal=sum((x'*invC).*x',2);
p = logdenom-0.5*mahal;
p = p';