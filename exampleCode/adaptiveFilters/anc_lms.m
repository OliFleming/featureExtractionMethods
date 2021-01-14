%
% anc_lms Perform adaptive noise cancellation using the LMS algorithm.
%
% [output, error, mse, mu] = anc_lms(u, d, u_trigger, d_trigger, mu, order, tv_mu, reload_state)
%
% Author: Vesal Badee
%
% Perform adaptive noise cancellation using the LMS algorithm. System state
% (i.e. tap weights) are saved to lms_state.mat.
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
%   mu: step-size (0 < mu < 1). Default = 0.01.
%   order: Order of transversal filter. Default = 28.
%   tv_mu: 3 element vector with minimum desired mu, a small,
%               positive learning rate parameter, and maximum desired mu
%               for variable step-size option (eg. [0.0001 0.001 0.01]).
%               Default = Constant mu.
%   reload_state: Possible values: 0 (do NOT reload state) and 1 (reload state).
%               The state is automatically saved by the function and can be
%               reloaded to be used as an input to run the system from its
%               previous state. It contains w_hat (weights). Default = 0.
%
% Outputs
%   output: LMS filter output.
%   error: Estimation error vector.
%   mse: Mean Square Error vector.
%   mu: Step-size vector.
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
% mu = 0.01; % learning rate
% 
% [output, error, mse, lambda] = anc_lms(input, desired, [], [], mu, Nfilter)
% 
% plot(t,input,t,desired,t,output,':');
% legend('input','desired','output')
%
% Modifications
% June 14, 2005: VB, Created
% June 24, 2005: Added PCG trigger
% August 9, 2005: Added mse
% December 7, 2005: Added variable mu and reload state
%
% Version 0.6

function [output, error, mse, mu] = anc_lms(u, d, u_trigger, d_trigger, mu, order, tv_mu, reload_state)

if (nargin<3) u_trigger = []; end
if (nargin<4) d_trigger = []; end
if (nargin<5) mu = 0.01; end
if (nargin<6) order = 28; end
if (nargin<7) tv_mu = []; end
if (nargin<8) reload_state = 0; end

if (isempty(mu)) mu = 0.01; end
if (isempty(order)) order = 28; end

% **************  A few error checks  **************

% if (mu>maxmu(u, order))
%     error = NaN;
%     errordlg('Error: Mu is larger than maximum mu allowed.');
%     return;
% end

if (size(u,1)==1 && size(u,2)>1) % row vector
    u = u'; % convert to column vector
end

if (size(u,1)>1 && size(u,2)>1) % matrix
    output = NaN; error = NaN; mse = NaN; w_final = NaN; mu = NaN;
    errordlg('Error: Input u must be a vector');
    return;
end

if (length(u) ~= length(d))
    output = NaN; error = NaN; mse = NaN; w_final = NaN; mu = NaN;
    errordlg('Error: Inputs u and d must have equal lengths');
    return;
end

if (~isempty(u_trigger))
    if (length(u) ~= length(u_trigger))
        output = NaN; error = NaN; mse = NaN; w_final = NaN; mu = NaN;
        errordlg('Error: u_trigger and u inputs must have equal lengths');
        return;
    end
end

if (~isempty(d_trigger))
    if (length(u) ~= length(d_trigger))
        output = NaN; error = NaN; mse = NaN; w_final = NaN; mu = NaN;
        errordlg('Error: d_trigger and d inputs must have equal lengths');
        return;
    end
end

% **************  End of error checks  **************

% Initialization of LMS algorithm

% w_hat: tap-weight vector
% psi_hat: derivative of w wrt mu

if (reload_state)
    load 'lms_state.mat';
else
    w_hat = zeros(order, 1);
end

% Initialize psi_hat only if mu is time-varying.
if (~isempty(tv_mu))
    psi_hat = zeros(order, 1);
end

mu_all = [];
output = [];
error = [];
mse = [];

% The next few lines are for the string to be displayed in the waitbar.
cr = sprintf('\n');
if (~isempty(d_trigger) & ~isempty(u_trigger)) triggerString = 'u and d';
elseif (~isempty(d_trigger)) triggerString = 'd';
elseif (~isempty(u_trigger)) triggerString = 'u';
else triggerString = 'None'; end
if (isempty(tv_mu)) muString = num2str(mu);
else muString = [num2str(tv_mu(1)) ' - ' num2str(tv_mu(3)) ...
        ', \alpha = ' num2str(tv_mu(2))]; end
waitbarMsg = ['Performing LMS ANC:' cr cr ...
        '\mu = ' muString cr 'Order = ' num2str(order) cr ...
        'Triggering: ' triggerString cr];

h = waitbar(0,waitbarMsg); 

for i=1:length(u)
    waitbar(i/length(u));   
    
    mu_all = [mu_all; mu];
    
    if (i >= order)
        e = d(i) - w_hat'*u(i-order+1:i);
        output = [output; w_hat'*u(i-order+1:i)];
    else
        tempU = [u(1:i); zeros(order-i,1)];
        e = d(i) - w_hat'*tempU;
        output = [output; w_hat'*tempU];
    end
    
    error = [error; e];
    mse = [mse; mean(error.^2)];
    
    % Adapt enable on u_trigger except during d_trigger
    if (isempty(u_trigger) | u_trigger(i))
        if (isempty(d_trigger) | ~d_trigger(i))
            if (i >= order)
                w_hat = w_hat + mu*u(i-order+1:i)*conj(e);

                % Update mu if it is time-varying.
                if (~isempty(tv_mu))
                    mu = mu + (tv_mu(2)*real(psi_hat'*u(i-order+1:i)*conj(e)));
                    if (mu>tv_mu(3)) mu=tv_mu(3); end
                    if (mu<tv_mu(1)) mu=tv_mu(1); end
                    psi_hat = ((eye(order,order) - (mu*u(i-order+1:i)*u(i-order+1:i)'))*psi_hat) + ...
                        (u(i-order+1:i)*conj(e));
                end
            else
                tempU = [u(1:i); zeros(order-i,1)];
                w_hat = w_hat + mu*tempU*conj(e);
                
                % Update mu if it is time-varying.
                if (~isempty(tv_mu))
                    mu = mu + (tv_mu(2)*real(psi_hat'*tempU*conj(e)));
                    if (mu>tv_mu(3)) mu=tv_mu(3); end
                    if (mu<tv_mu(1)) mu=tv_mu(1); end
                    psi_hat = ((eye(order,order) - (mu*tempU*tempU'))*psi_hat) + ...
                        (tempU*conj(e));
                end
            end            
        end
    end
end
close(h);

mu = mu_all;

save lms_state.mat w_hat