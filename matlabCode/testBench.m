
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
phase = 0
dc = 0
noise = 0.1;
type = 'sin';

%generate waveform.
[t,s] = waveformGen(time,interval,amp,frequency,phase,dc,noise,type)


%Load test data;

%Initial feature Variables
winsize = 100;
wininc = 50;
datawin = 1;
dispstatus = 1;

%---- RMS ----
rms_feat = getrmsfeat(s,winsize,wininc,datawin,dispstatus);


%---- SSC ----
deadzone = 0.01;
ssc_feat = getsscfeat(s,deadzone,winsize,wininc,datawin,dispstatus);

disp('rms_feat\n')
disp(rms_feat)

disp('ssc_feat\n')
disp(ssc_feat)




%------------- END OF CODE --------------