%
% anc_rls Perform adaptive noise cancellation using the RLS algorithm.
%
% [output, error, mse, lambda] = anc_rls(u, d, u_trigger, d_trigger, lambda, order, tv_lambda, reload_state)
%
% Author: Vesal Badee
%
% Perform adaptive noise cancellation using the RLS algorithm. System state
% (i.e. tap weights, gain factor and inverse correlation matrix) are saved 
% to rls_state.mat.
%
% Inputs
%   u: Reference vector (eg. thoracic maternal ECG).
%   d: Primary/Desired vector (eg. abdominal composite ECG).
%   u_trigger: Vector of ones and zeros. Ones indicate event location for 
%               adapt enable on u input (eg. during maternal ECG).
%               Default = [] (no triggering).
%   d_trigger: Vector of ones and zeros. Ones indicate event location for 
%               adapt disable on d input (eg. during fetal ECG). 
%               Default = [] (no triggering).
%   lambda: Forgetting factor (0 < lambda < 1). Default = 0.999.
%   order: Order of transversal filter. Default = 28.
%   tv_lambda: Time-varying lambda. Three element vector with minimum desired lambda,
%               a small, positive learning rate parameter, and maximum desired lambda
%               for variable forgetting factor option (eg. [0.999 0.001 0.9999]).
%               Default = [] (constant lambda).
%   reload_state: Possible values: 0 (do NOT reload state) and 1 (reload state).
%               The state is automatically saved by the function and can be
%               reloaded to be used as an input to run the system from its 
%               previous state. It contains w_hat (weights), k (gain factor),
%               and P (inverse correlation matrix). Default = 0.
%
% Outputs
%   output: RLS filter output.
%   error: Estimation error vector.
%   mse: Mean Square Error vector.
%   lambda: Forgetting factor vector.
%
% Reference
% Haykin, S., "Adaptive filter theory", Prentice-Hall, 1996
%
% EXAMPLE: Adaptive filtering of a square wave to a sine wave
% fs = 100; % sampling rate
% t = (0:1/fs:10)'; % time scale
% N = length(t); % number of samples
% fc = 2; % frequency of input square wave and sine wave desired output
% input = square(2*pi*fc*t) + 0.1*rand(size(t)); % square wave
% desired = sin(2*pi*fc*t) + 0.1*rand(size(t)); % sine wave (can be fc, 3*fc, 5*fc, etc.)
% Nfilter = 32; % filter length
% 
% lambda = 0.99; % forgetting factor
% 
% [output, error, mse, lambda] = anc_rls(input, desired, [], [], lambda, Nfilter)
% 
% plot(t,input,t,desired,t,output,':');
% legend('input','desired','output')
%
% Modifications
% June 13, 2005: VB, Created
% June 24, 2005: Added PCG trigger
% August 9, 2005: Added mse
% December 7, 2005: Added variable lambda and reload state
%
% Version 0.6

function [output, error, mse, lambda] = anc_rls(u, d, u_trigger, d_trigger, lambda, order, tv_lambda, reload_state)

if (nargin<3) u_trigger = []; end
if (nargin<4) d_trigger = []; end
if (nargin<5) lambda = 0.999; end
if (nargin<6) order = 28; end
if (nargin<7) tv_lambda = []; end
if (nargin<8) reload_state = 0; end

if (isempty(lambda)) lambda = 0.999; end
if (isempty(order)) order = 28; end

% **************  A few error checks  **************

if (size(u,1)==1 && size(u,2)>1) % row vector
    u = u'; % convert to column vector
end

if (size(u,1)>1 && size(u,2)>1) % matrix
    output = NaN; error = NaN; mse = NaN; w_final = NaN; lambda = NaN;
    errordlg('Error: Input u must be a vector');
    return;
end

if (length(u) ~= length(d))
    output = NaN; error = NaN; mse = NaN; w_final = NaN; lambda = NaN;
    errordlg('Error: Inputs u and d must have equal lengths');
    return;
end

if (~isempty(u_trigger))
    if (length(u) ~= length(u_trigger))
        output = NaN; error = NaN; mse = NaN; w_final = NaN; lambda = NaN;
        errordlg('Error: u_trigger and u inputs must have equal lengths');
        return;
    end
end

if (~isempty(d_trigger))
    if (length(u) ~= length(d_trigger))
        output = NaN; error = NaN; mse = NaN; w_final = NaN; lambda = NaN;
        errordlg('Error: d_trigger and d inputs must have equal lengths');
        return;
    end
end

% **************  End of error checks  **************

% Initialization of RLS algorithm

% k(n): M-by-1 gain vector
% P(n): inverse of M-by-M correlation matrix
% w_hat: tap-weight vector
% S: derivative of P wrt lambda
% psi_hat: derivative of w wrt lambda

if (reload_state)
    load 'rls_state.mat';
else
    w_hat = zeros(order, 1);
    k = zeros(order, 1);
    P = eye(order, order)/(0.01*0.1e-1); % Using var(data) = 0.1, delta = 10^-4
end

% Initialize S and psi_hat only if lambda is time-varying.
if (~isempty(tv_lambda))
    S = zeros(order,order);
    psi_hat = zeros(order, 1);
end

lambda_all = [];
output = [];
error = [];
mse = [];

% The next few lines are for the string to be displayed in the waitbar.
cr = sprintf('\n');
if (~isempty(d_trigger) & ~isempty(u_trigger)) triggerString = 'u and d';
elseif (~isempty(d_trigger)) triggerString = 'd';
elseif (~isempty(u_trigger)) triggerString = 'u';
else triggerString = 'None'; end
if (isempty(tv_lambda)) lambdaString = num2str(lambda);
else lambdaString = [num2str(tv_lambda(1)) ' - ' num2str(tv_lambda(3)) ...
        ', \alpha = ' num2str(tv_lambda(2))]; end
waitbarMsg = ['Performing RLS ANC:' cr cr ...
        '\lambda = ' lambdaString cr 'Order = ' num2str(order) cr ...
        'Triggering: ' triggerString cr];

h = waitbar(0,waitbarMsg); 

for i=1:length(u)
    waitbar(i/length(u));
    
    lambda_all = [lambda_all; lambda];
    
    if (i >= order)
        kNum = (1/lambda)*P*u(i-order+1:i);
        kDen = 1 + ((1/lambda)*u(i-order+1:i)'*P*u(i-order+1:i));
        k = kNum/kDen;
        
        output = [output; w_hat'*u(i-order+1:i)];
        xi = d(i) - w_hat'*u(i-order+1:i);
        
    else
        tempU = [u(1:i); zeros(order-i,1)];
        kNum = (1/lambda)*P*tempU;
        kDen = 1 + ((1/lambda)*tempU'*P*tempU);
        k = kNum/kDen;
        
        output = [output; w_hat'*tempU];
        xi = d(i) - w_hat'*tempU;        
    end    
    
    error = [error; xi];
    mse = [mse; mean(error.^2)];
    
    % Adapt enable on u_trigger except during d_trigger
    if (isempty(u_trigger) | u_trigger(i))
        if (isempty(d_trigger) | ~d_trigger(i))
            w_hat = w_hat + k*conj(xi);
        end
    end
    
    if (i >= order)
        P = (1/lambda)*P - (1/lambda)*k*u(i-order+1:i)'*P;
%        figure(1)
%        mesh(P)
%        pause(0.1)
        
        % Update lambda if it is time-varying.
        if (~isempty(tv_lambda))
            lambda = lambda + tv_lambda(2)*real(psi_hat'*u(i-order+1:i)*conj(xi));
            if (lambda>tv_lambda(3)) lambda=tv_lambda(3); end
            if (lambda<tv_lambda(1)) lambda=tv_lambda(1); end
            S = ((1/lambda)*(eye(order,order) - (k*u(i-order+1:i)'))*S*(eye(order,order) - ...
                (u(i-order+1:i)*k'))) + ((1/lambda)*k*k') - ((1/lambda)*P);
            psi_hat = ((eye(order,order) - (k*u(i-order+1:i)'))*psi_hat) + ...
                (S*u(i-order+1:i)*conj(xi));
        end
    else
        tempU = [u(1:i); zeros(order-i,1)];
        P = (1/lambda)*P - (1/lambda)*k*tempU'*P;
        
        % Update lambda if it is time-varying.
        if (~isempty(tv_lambda))
            lambda = lambda + tv_lambda(2)*real(psi_hat'*tempU*conj(xi));
            if (lambda>tv_lambda(3)) lambda=tv_lambda(3); end
            if (lambda<tv_lambda(1)) lambda=tv_lambda(1); end
            S = ((1/lambda)*(eye(order,order) - (k*tempU'))*S*(eye(order,order) ...
                - (tempU*k'))) + ((1/lambda)*k*k') - ((1/lambda)*P);
            psi_hat = ((eye(order,order) - (k*tempU'))*psi_hat) + ...
                (S*tempU*conj(xi));
        end
    end
end
close(h);

lambda = lambda_all;

save rls_state.mat w_hat k P