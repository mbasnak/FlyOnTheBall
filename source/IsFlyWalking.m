
function [percentMoving,moving] = IsFlyWalking(rawData,voltThresh)
% This function determines the percentage of the trial that the fly spends
% walking. For this, it takes as inputs the raw voltage data (rawData), as
% well as a voltThresh above which mvt will be considered (this threshold
% will if possible be determined by a previous assessment of the noise of
% the voltage reading, using the function assesNoise).
% It outputs both the percentage of the trial that the fly spends moving,
% as well as the frame by frame assessment of whether the fly is moving or
% not.



%we need to compute a vector that is the difference of each pos with the
%previous one
for j = 1:size(rawData,2)
    for i = 2:length(rawData)
        mvtData(i-1,j) = rawData(i,j)-rawData(i-1,j);
    end
end


%we need to see whether in any of the 3 FicTrac outputs the voltage is
%above that threshold, and in that case set the "moving" variable to 1
%(otherwise set it to 0)

for i = 1:length(mvtData)
    if mvtData(i,1) > voltThresh | mvtData(i,2) > voltThresh | mvtData(i,3) > voltThresh
        moving(i) = 1;
    else
        moving(i) = 0;
    end
end

% Determine what is the percentage of the trial during which the fly is
% moving and save it

percentMoving = 100*(sum(moving)/length(moving));


end