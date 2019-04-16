function [ out ] = lowPassFilter( data,  lowPassCutOff , sampleRate )
%LOWPASSFILTER Filters the data using a lowpass butter filter of specified cut out
%   
%   INPUTS
%   data - trace to be filtered
%   lowPassCutOff - value (Hz) that will be the top limit of the filter
%   sampleRate - rate data is sampled at to allow correct conversion in to
%   Hz.
% 
%   OUTPUT
%   out - filtered version of data
%   
%   Yvette Fisher 1/2018

% low pass filter the data trace :
% fprintf('\nLow-pass filtering at %d Hz\n',lowPassCutOff);

% build butter function
[b,a] = butter(1 , lowPassCutOff / ( sampleRate / 2),'low');

% filter data using butter function
out = filtfilt( b ,a , data);

end