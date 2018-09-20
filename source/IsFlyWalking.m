
function [percentMoving,moving] = IsFlyWalking(rawData)

%we need to compute a vector that is the difference of each pos with the
%previous one
for j = 1:size(rawData,2)
    for i = 2:length(rawData)
        mvtData(i-1,j) = rawData(i,j)-rawData(i-1,j);
    end
end

%we need to set a minimum voltage change to imply the animal is walking
voltThresh = 0.002;

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