%% Analysis of the panel/FicTrac data
% To run it, you should be standing inside the folder of the fly whose data
% you're trying to analyze.

%%% This code analyses the output of the panels and the FicTrac data
close all; clear all;

% prompt the user to select the file to open and load it.
[file,path] = uigetfile();
load([path,file]);

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;

%% Subset acquisition of x and y pos, as well as FicTrac data

data.xPanelVolts =  rawData (:,xPanels); 
VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
data.xPanelPos = round ((data.xPanelVolts  * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels

data.yPanelVolts =  rawData (:, yPanels);
VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
maxValY = 1;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2
data.yPanelPos = round ((data.yPanelVolts  * maxValY) /VOLTAGE_RANGE);

%FicTrac data
data.ficTracAngularPosition = rawData ( : , headingFly); 
data.ficTracIntx = rawData ( : , xFly); 
data.ficTracInty = rawData ( : , yFly); 


%% Downsample, unwrap and smooth position data, then get velocity and smooth

[smoothed] = singleTrialVelocityAnalysis(data,1000);
% for most experiments I have used 1000 Hz as a sample rate, and it is what
% I will use from now on, so that's how I'll leave it, but this could be
% changed in case of need

%% Output in degrees of the Panels position

% Pos x=92 is 0 deg (ie facing the fly), I measured this empirically

pxToDeg = 360/97; % There are 97 possible positions (the last one = first one) and this represents 360 deg
posToDeg = zeros(1,length(data.xPanelPos));


% Convert from xpos to degrees, knowing that xpos 5 = 0 deg
for i=1:length(data.xPanelPos)
    if data.xPanelPos(i) ==93 | data.xPanelPos(i) ==94 | data.xPanelPos(i) ==95 | data.xPanelPos(i) ==96 | data.xPanelPos(i) ==97
        posToDeg(i) = (data.xPanelPos(i)-92)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        posToDeg(i) = (data.xPanelPos(i)+4)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end


% Create a time vector in sec
time = linspace(0,(size(rawData,1)/1000),size(rawData,1)); %1000 is our sampling rate

% Plot the position of the stimulus in degrees for the trial
figure,
plot(time,posToDeg)
ylim([0 365]);
ylabel('Position of the stimulus (deg)');
xlabel('Time (s)');
title({'Position of the stimulus in deg as a function of time', typeOfStim}, 'Interpreter', 'none');

saveas(gcf,strcat(path,'PanelOutput_ExpNum', file(end-4), '.png'))
%% Probability density function of the stimulus position

% Remapping the positions to span -180 to 180 deg
remapPosToDeg = posToDeg;
for i = 1:length(remapPosToDeg)   
    if remapPosToDeg(i) > 180
        remapPosToDeg(i) = remapPosToDeg(i) - 360;  
    end   
end

% Plot the histogram and probability density
[counts] = histcounts(remapPosToDeg,20);
probabilities = counts./sum(counts);
degs = linspace(-180,180,length(counts));

figure,
subplot(1,2,1)
histogram(remapPosToDeg,20,'Normalization','probability')
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Histogram of the stimulus position');
ylabel('Probability'); xlabel('Stimulus position (deg)');

subplot(1,2,2),
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')
suptitle(typeOfStim);
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Probability density of the stimulus position');
ylabel('Probability density'); xlabel('Stimulus position (deg)');
if contains(typeOfStim,'closed_loop')
%if isequal(typeOfStim, 'closed_loop_bar') | isequal(typeOfStim, 'dark_closed_loop_bar') 
   hold on 
   %add line showing the start pos. of the stim. to see whether the fly
   %just kept it there
    remapStartPos = 0;
     if startPos(1) ==93 | startPos(1) ==94 | startPos(1) ==95 | startPos(1) ==96 | startPos(1) ==97
        startingPos = (startPos(1)-92)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        startingPos = (startPos(1)+4)*pxToDeg; % Correct the offset and multiply by factor to get deg
     end
    
    if startingPos > 180
       remapStartPos = startingPos - 360;
    else
       remapStartPos = startingPos;
    end 
       line([remapStartPos remapStartPos],[0 max(probabilities)+0.05],'Color',[1 0 0])

%add rectangle showing where the panel is off (and the bar can't be seen)
noPanel = [17:23]; %xpos where a 2 px bar can't be seen.
noPanelDeg = (noPanel+4)*pxToDeg;
l1 = line([noPanelDeg(1) noPanelDeg(1)],[0 max(probabilities)+0.05]);
l2 = line([noPanelDeg(7) noPanelDeg(7)],[0 max(probabilities)+0.05]);
set([l1 l2],'Color',[.5 .5 .5]);
patch([noPanelDeg(1) noPanelDeg(7) noPanelDeg(7) noPanelDeg(1)], [0 0 max(probabilities)+0.05 max(probabilities)+0.05],[.5 .5 .5],'FaceAlpha',0.3)

end

saveas(gcf,strcat(path,'ProbabilityDensityStimPosition_ExpNum', file(end-4), '.png'))

%%  How much is the fly moving?

[percentMoving, moving] = IsFlyWalking(rawData,0.002);

%add a zero before moving start, for the frame 1 to have a "is not moving"
%assigned
moving = [0,moving];

%change to boolean
Moving = true(size(moving,1),size(moving,2));

for i = 1:length(moving)
    if moving(1,i) == 0
        Moving(i) = false;
    else
        Moving(i) = true;
    end
end

%% Probability density of the fly's position (for open loop trials we need it because the stimulus one is not informative of the fly's behavior)
%There is something I have to check because the fly's position seems like a
%mirror of that of the stimulus.


% flyPosToDeg = data.ficTracAngularPosition.*36; %if 10 V = 360 deg, then xV = 36x deg
% 
% remapFlyPosToDeg = flyPosToDeg;
% for i = 1:length(remapFlyPosToDeg)   
%     if remapFlyPosToDeg(i) > 180
%         remapFlyPosToDeg(i) = remapFlyPosToDeg(i) - 360;  
%     end   
% end
% 
% % Plot the histogram and probability density
% [countsFly] = histcounts(remapFlyPosToDeg,20);
% probabilitiesFly = countsFly./sum(countsFly);
% degsFly = linspace(-180,180,length(countsFly));
% 
% figure,
% subplot(1,2,1)
% histogram(remapFlyPosToDeg,20,'Normalization','probability')
% xlim([-180 180]); ylim([0 max(probabilitiesFly)+0.05]);
% title('Histogram of the fly heading');
% ylabel('Probability'); xlabel('Fly heading (deg)');
% 
% subplot(1,2,2),
% plot(degsFly,probabilitiesFly,'k')
% xlim([-180 180]); ylim([0 max(probabilitiesFly)+0.05]);
% title('Probability density of the fly heading');
% ylabel('Probability density'); xlabel('Fly heading (deg)');
% %Add the starting pos of the bar if the stimulus was a closed-loop bar
% if isequal(typeOfStim, 'closed_loop_bar') | isequal(typeOfStim, 'dark_closed_loop_bar') 
%    hold on
%    line([remapStartPos remapStartPos],[0 max(probabilitiesFly)+0.05],'Color',[1 0 0])
% end



%% Probability density keeping only those frames when the fly was moving

flyPosToDegMoving = data.ficTracAngularPosition(Moving).*36; %if 10 V = 360 deg, then xV = 36x deg

remapFlyPosToDegMoving = flyPosToDegMoving;
for i = 1:length(remapFlyPosToDegMoving)   
    if remapFlyPosToDegMoving(i) > 180
        remapFlyPosToDegMoving(i) = remapFlyPosToDegMoving(i) - 360;  
    end   
end

% Plot the histogram and probability density
[countsFlyMoving] = histcounts(remapFlyPosToDegMoving,20);
probabilitiesFlyMoving = countsFlyMoving./sum(countsFlyMoving);
degsFlyMoving = linspace(-180,180,length(countsFlyMoving));

figure,
subplot(1,2,1)
histogram(remapFlyPosToDegMoving,20,'Normalization','probability')
xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
title('Histogram of the fly heading');
ylabel('Probability'); xlabel('Fly heading (deg)');

subplot(1,2,2),
plot(degsFlyMoving,probabilitiesFlyMoving,'k')
suptitle(typeOfStim)
xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
title('Probability density of the fly heading');
ylabel('Probability density'); xlabel('Fly heading (deg)');
%Add the starting pos of the bar if the stimulus was a closed-loop bar
if contains(typeOfStim,'closed_loop')
%if isequal(typeOfStim, 'closed_loop_bar') | isequal(typeOfStim, 'dark_closed_loop_bar') 
   hold on
   line([remapStartPos remapStartPos],[0 max(probabilitiesFlyMoving)+0.05],'Color',[1 0 0])
   patch([noPanelDeg(1) noPanelDeg(7) noPanelDeg(7) noPanelDeg(1)], [0 0 max(probabilitiesFlyMoving)+0.05 max(probabilitiesFlyMoving)+0.05],[.5 .5 .5],'FaceAlpha',0.3)
end

saveas(gcf,strcat(path,'ProbabilityDensityFlyHeading_ExpNum', file(end-4), '.png'))
%% Forward velocity

forwardVelocity = smoothed.xVel;
meanVelocity = mean(forwardVelocity);
time = linspace(0,(length(rawData)/1000),length(forwardVelocity));

figure,
plot(time,forwardVelocity,'k')
hold on
hline = refline([0 meanVelocity]);
hline.Color = 'r'; hline.LineStyle = '--';
title({'Forward velocity of the fly', typeOfStim}, 'Interpreter', 'none');
xlabel('Time (s)')
ylabel('Velocity (mm/s)')

saveas(gcf,strcat(path,'ForwardVelocity_ExpNum', file(end-4), '.png'))