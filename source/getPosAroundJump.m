function [angPos] = getPosAroundJump(data,jumpFrames,aroundJumpSec)

%This function outputs a data structure for position and velocity around
%the bar jumps

%data is the original data you get in voltage
%jumpFrames will contain the frames when the jumps occurred
%aroundJumpSec is how many seconds before and after the jump you're
%interested in getting.

for i = 1:length(jumpFrames)
      
    %angular position
    aroundJumpa(:,i) = data(jumpFrames(i)-aroundJumpSec*1000:jumpFrames(i)+aroundJumpSec*1000,1);
    downsampled.aJa(:,i) = downsample(aroundJumpa(:,i),1000/25);
    downsRad.aJa(:,i) = downsampled.aJa(:,i) .* 2 .* pi ./ 10;
    unwrapped.aJa(:,i) = unwrap(downsRad.aJa(:,i));
    smoothed.aJa(:,i) = smoothdata(unwrapped.aJa(:,i),'rlowess',10); 
    angPosition(:,i) = (smoothed.aJa(:,i) / (2*pi)) * 360;
    
     %remap the angular position acording to our set-up
 flyPosPixels = angPosition/3.75;
 angPos = zeros(size(flyPosPixels,1),size(flyPosPixels,2));

% Convert from xpos to degrees, knowing that xpos 70 = 0 deg
for j = 1:numel(angPos)
        if flyPosPixels(j) == 70
        angPos(j) = 0;
        elseif flyPosPixels(i) >70 
        angPos(j) = (flyPosPixels(j)-70)*3.75; % Correct the offset and multiply by factor to get deg
        else
        angPos(j) = (flyPosPixels(j)+27)*3.75; % Correct the offset and multiply by factor to get deg
        end
end
           
end