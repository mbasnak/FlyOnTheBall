function [DATA] = getDataAroundJump(data,jumpFrames,aroundJumpSec,sizeBall,samplingRate)

%This function outputs a data structure for position and velocity around
%the bar jumps

%data is the original data you get in voltage
%jumpFrames will contain the frames when the jumps occurred
%aroundJumpSec is how many seconds before and after the jump you're
%interested in getting.
%sizeBall is the size of the ball I used, in mm, to calculate the movement
%in metrix units aproprietaly

for i = 1:length(jumpFrames)
    
    %x position and forward walking
    aroundJumpIntx(:,i) = data(jumpFrames(i)-aroundJumpSec*samplingRate:jumpFrames(i)+aroundJumpSec*samplingRate,3);
    %The - sign is necessary to invert the x axis.
    downsampled.aJ(:,i) = downsample(aroundJumpIntx(:,i),samplingRate/25);
    downsRad.aJ(:,i) = downsampled.aJ(:,i) .* 2 .* pi ./ 10;
    unwrapped.aJ(:,i) = unwrap(downsRad.aJ(:,i));
    smoothed.aJ(:,i) = smoothdata(unwrapped.aJ(:,i),'rlowess',25); 
    deg.aJ(:,i) = smoothed.aJ(:,i) * (sizeBall/2);
    diff.aJ(:,i) = gradient(deg.aJ(:,i)).* 25;
    DATA.forwardVel(:,i) = smoothdata(diff.aJ(:,i),'rlowess',15);
    
    %angular velocity
    aroundJumpa(:,i) = data(jumpFrames(i)-aroundJumpSec*samplingRate:jumpFrames(i)+aroundJumpSec*samplingRate,1);
    downsampled.aJa(:,i) = downsample(aroundJumpa(:,i),samplingRate/25);
    downsRad.aJa(:,i) = downsampled.aJa(:,i) .* 2 .* pi ./ 10;
    unwrapped.aJa(:,i) = unwrap(downsRad.aJa(:,i));
    smoothed.aJa(:,i) = smoothdata(unwrapped.aJa(:,i),'rlowess',25); 
    smoothed.angPos(:,i) = (smoothed.aJa(:,i) / (2*pi)) * 360;
    diff.aJa(:,i) = gradient(smoothed.angPos(:,i)).* 25;
    DATA.angVel(:,i) = smoothdata(diff.aJa(:,i),'rlowess',15);
    
end