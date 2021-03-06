function [t,s] = waveformGen(time,interval,amp,frequency,phase,dc,noise,type)

%FUNCTION_NAME - Generates a Specified waveform
%
% Syntax:  [t,s] = WaveformGen(time,interval,amp,frequency,phase,dc,noise,type)
%
% Inputs:
%    time - number of data points to producs
%    amp - signal amplitude
%    frequency - frequency of wave if applicable
%    phase - Phase offset of wave in radians
%    dc - dc offset aplied to wave
%    noise - noise amplitude to super impose over signal
%    type - Type pf wave to be produced
%           <rand> - 
%           <sin> - 
%           <cos> - 
%           <square> - 
%           <triangle> - 
%    
%
% Outputs:
%    t - time stamps
%    s - signal value
%
% Example: 
%    WaveformGen(10,0.01,1,1,0,0,0.1,'triangle');
%    WaveformGen(10,0.01,1,1,0,0,0,'square');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: N/A

% Author: Lachlan Smith
% Work address: 8 Little Queen Street, Chippendale NSW 2008.
% email: lsmi5655@uni.sydney.edu.au
% Website: https://www.sydney.edu.au/engineering
% Janurary 2021; Last revision: 14-1-2021

%------------- BEGIN CODE --------------
t = 0:interval:time;

switch type

case 'rand'
  s = -amp + (2*amp)*rand(1,length(t));
  
case 'sin'
  s = dc + amp * (sin( (2*pi*frequency*t + phase) ));

case 'cos'
  s = dc + amp * (sin( (2*pi*frequency*t + phase) ));
  
case 'square' 
  s = dc + amp * square( (2*pi*frequency*t + phase), 50 );
  
case 'triangle'
  s = dc + amp * sawtooth( (2*pi*frequency*t + phase), 0.5);
  
 
otherwise
  disp('Invalid Type');
  
end

%noise addition
s = s + (-noise + (2*noise)*rand(1,length(t)));
s = s';
t = t';

%Plot and save waveform
inpt = input ('View Waveform? (Y/N): ', 's');
if strcmp(inpt, 'y') || strcmp(inpt, 'Y')
  plot(t, s)
  grid
end
 

%Save saveform
inpt = input ('Save Waveform? (Y/N): ', 's');
if strcmp(inpt, 'y') || strcmp(inpt, 'Y')
  inpt = input ('File Name: ', 's');
  
  A = [t,s];  
  writematrix(A,inpt);
 
end 
  
  
%------------- END OF CODE --------------