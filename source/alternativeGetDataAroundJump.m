function [DATA] = alternativeGetDataAroundJump(data,jumpFrames,aroundJumpSec,sizeBall)

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
    aroundJumpIntx(:,i) = -data(jumpFrames(i)-aroundJumpSec*1000:jumpFrames(i)+aroundJumpSec*1000,3);
    %The - sign is necessary to invert the x axis.
    downsampled.aJ(:,i) = downsample(aroundJumpIntx(:,i),1000/25);
    downsRad.aJ(:,i) = downsampled.aJ(:,i) .* 2 .* pi ./ 10;
    unwrapped.aJ(:,i) = unwrap(downsRad.aJ(:,i));
    deg.aJ(:,i) = unwrapped.aJ(:,i) * (sizeBall/2);
    DATA.forwardVel(:,i) = filter(-smooth_diff(10),1,deg.aJ(:,i))/0.04;
    
    %angular velocity
    aroundJumpa(:,i) = data(jumpFrames(i)-aroundJumpSec*1000:jumpFrames(i)+aroundJumpSec*1000,1);
    downsampled.aJa(:,i) = downsample(aroundJumpa(:,i),1000/25);
    downsRad.aJa(:,i) = downsampled.aJa(:,i) .* 2 .* pi ./ 10;
    unwrapped.aJa(:,i) = unwrap(downsRad.aJa(:,i));
    smoothed.angPos(:,i) = (unwrapped.aJa(:,i) / (2*pi)) * 360;
    DATA.angVel(:,i) = filter(-smooth_diff(10),1,smoothed.angPos(:,i))/0.04;

end