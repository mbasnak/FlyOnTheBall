function decodedYPos = processPanelDataY( data , panelParams )
% PROCESSPANELDATAY
%
%   "The DAC0, and DAC1 voltages will be values between 0V – 10V. 
%   The size of the voltage steps is set by the number of frames in X, Y
%   for the current pattern (e.g. 96 frames in one channel would lead to 
%   voltage steps of 10/96 V; frame index 48 would be roughly 5 V). These 
%   outputs are consistent and can be used to recover the exact frame position
%   for patterns of at least 500 frames (can only be done approximately for larger patterns). 
%   These instantaneous pattern positions are essential for off-line analysis of closed-loop behaviors, 
%   and are useful for validating open-loop experimental protocols.  
% store aquired Panel x and y information into data" - Panel creators
%
% Y Fisher 3/2017
VOLTAGE_RANGE = 10; % Volts

% Default y dim number 
maxValY = 3;% This is true as of 3/17 for all patterns: 1 = OFF, 2 = pattern or 3 = ON for y

% EXCEPTION FOR PATTERN 10: CONTRAST DATA with 9 channels!!
 if( panelParams.patternNum == 10) % diff contrast channel
     maxValY = 9;% This is true as of 3/17 for all patterns: 1 = OFF, 2 = pattern con1  or 3 = ON, 4 - 9 (other contasts)
 end
 
 % Exception for pattern 11:  bar diff width, 7 channels
 if( panelParams.patternNum == 11 || panelParams.patternNum == 12) % diff stripe width channel
     maxValY = 7;% stripe width 0, 1, 3, 7, 13, 25, 120
 end
 
decodedYPos = round (( data  * maxValY) /VOLTAGE_RANGE);


% 
% %% REMEMBER
% maxVal = 3;
% minVal = 1;
% frames = 2;
% voltsPerStep = (maxVal - minVal)/(frames);
% fr = round((rawData - minVal)./voltsPerStep);
end