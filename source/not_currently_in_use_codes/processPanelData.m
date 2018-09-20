function [fr, angle] = processPanelData(rawData,stim)
% PROCESSPANELDATA
%
%Obtained by Jenny Lu 3/2017

%% REMEMBER
% Pattern x=1 starts eastward most.
% Pattern goes counterclockwise thereafter.
%5 panels = 75;
%7 panels = 105;
if ~isfield(stim, 'initialAngle')
    initialAngle = 105; %CHECK!!!!
else
    initialAngle = stim.initialAngle;
end
maxVal = 5;
minVal = 0;
frames = stim.frames;
voltsPerStep = (maxVal - minVal)/(frames-1);
fr = round((rawData - minVal)./voltsPerStep);
pixelAngle = 360./96;
arenaAngle = frames*pixelAngle;
angle = (initialAngle-((fr-1)+stim.barWidth/2).*pixelAngle); % need to subtract for the bar width
%% CHECK
angle = wrapTo180(angle);
if arenaAngle < 360
    halfArena = arenaAngle./2;
    indexOver = angle < -halfArena;
    angle = angle + indexOver.*arenaAngle;
end