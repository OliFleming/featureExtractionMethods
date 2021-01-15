
% testBench.m
%
% Signal Feature extraction test bench
%
%
%
% Author: Lachlan Smith
% Work address: 8 Little Queen Street, Chippendale NSW 2008.
% email: lsmi5655@uni.sydney.edu.au
% Website: https://www.sydney.edu.au/engineering
% Janurary 2021; Last revision: 14-1-2021

%------------- BEGIN CODE --------------

%initial Matlab setup
clc;
clear;
clf;
close;

%Waveform Generation variables
time = 10;
interval = 0.01;
amp = 5;
frequency = 1;
phase = 0;
dc = 0;
noise = 0.1;
type = 'sin';

%generate waveform.
%insert waveform generation here
A = readmatrix('smallnoise_sin.csv');

t = A(:,1);
s = A(:,2);
plot(t, s)
grid

%Load test data;

%Initial feature Variables
winsize = 2000;
wininc = 2000;
datawin = 1;
dispstatus = 1;
deadzone = 0.01;

%---- RMS ----
rms_feat = getrmsfeat(s,winsize,wininc,datawin,dispstatus);


%---- SSC ----
ssc_feat = getsscfeat(s,deadzone,winsize,wininc,datawin,dispstatus);

zc_feat = getzcfeat(s,deadzone,winsize,wininc,1);

disp('rms_feat')
disp(rms_feat)

disp('ssc_feat')
disp(ssc_feat)

disp('zcc_feat')
disp(zc_feat)



%------------- END OF CODE --------------