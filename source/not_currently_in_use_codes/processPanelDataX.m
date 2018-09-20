function [ decodedXPos ] = processPanelDataX( data , panelParams )
% PROCESSPANELDATAX
%
%  WARNING currently slight off by ~ 1 pos value for this X decode
%   "The DAC0, and DAC1 voltages will be values between 0V – 10V. 
%   The size of the voltage steps is set by the number of frames in X, Y
%   for the current pattern (e.g. 96 frames in one channel would lead to 
%   voltage steps of 10/96 V; frame index 48 would be roughly 5 V). These 
%   outputs are consistent and can be used to recover the exact frame position
%   for patterns of at least 500 frames (can only be done approximately for larger patterns). 
%   These instantaneous pattern positions are essential for off-line analysis of closed-loop behaviors, 
%   and are useful for validating open-loop experimental protocols.  
% store aquired Panel x and y information into data" - Panel creators
% Y Fisher 3/2017, updated 5/2017
VOLTAGE_RANGE = 10; % Volts

%maxValX = 56;% This should be true as of 3/17 for all patterns
%maxValX = 72;%Updated when switched to 230 deg arena 10/2017
maxValX = 96; % updated when switch to having patterns for whole 360 degree world 12/2017


% edited to handled 270 wrapped world again....all new patterns 14+ (and 2) are for
% 270 world:
if ( (panelParams.patternNum >= 14) || (panelParams.patternNum == 2)  )
    maxValX = 72;
end

 % Exception for pattern 24:  dark bar for open loop
 if( panelParams.patternNum == 24) % diff stripe width channel
     maxValX =  96 ;% pattern.x_num
 end

%  % Exception for pattern 13:  dot at diff width
%  if( panelParams.patternNum == 13) && (panelParams.positionFuncNumX < 28) % diff stripe width channel
%      maxValX =  840 ;% pattern.x_num 
%  end
%  
%   % Exception for pattern 14:  dot at diff width
%  if( panelParams.patternNum == 14 ) % diff stripe width channel
%      maxValX = 728 ;% pattern.x_num 
%  end

decodedXPos = round (( data  * maxValX ) /VOLTAGE_RANGE);

end

