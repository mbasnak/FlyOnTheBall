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
data.ficTracIntx = -rawData ( : , xFly); 
%data.ficTracIntx = -(rawData ( : , xFly)); 
data.ficTracInty = rawData ( : , yFly); 


%% Downsample, unwrap and smooth position data, then get velocity and smooth

[smoothed] = singleTrialVelocityAnalysis(data,1000);
% for most experiments I have used 1000 Hz as a sample rate, and it is what
% I will use from now on, so that's how I'll leave it, but this could be
% changed in case of need


%% Forward velocity analysis

% The forward velocity is a good indicative of whether the fly is walking
% well in that trial or not

forwardVelocity = smoothed.xVel;
meanVelocity = mean(forwardVelocity);
time = linspace(0,(length(rawData)/1000),length(forwardVelocity));

figure,
subplot(2,1,1)
plot(time,forwardVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
hold on
hline = refline([0 meanVelocity]);
hline.Color = 'r'; hline.LineStyle = '--';
rline = refline([0 0]);
rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
title('Forward velocity of the fly', 'Interpreter', 'none');
xlabel('Time (s)')
ylabel('Velocity (mm/s)')
legend('Mean forward velocity');

subplot(2,1,2)
histogram(forwardVelocity,'FaceColor',[.4 .2 .6])
title('Distribution of forward velocities');
xlabel('Forward velocity (mm/s)');
ylabel('Frequency');

%saveas(gcf,strcat(path,'ForwardVelocity_ExpNum', file(11:end-4), '.png'))
%saveas(gcf,strcat(path,'ForwardVelocity_ExpNum', file(11:end-4), '.svg'))


%% Angular velocity analysis

angularVelocity = smoothed.angularVel;
meanAngVelocity = mean(angularVelocity);
time = linspace(0,(length(rawData)/1000),length(angularVelocity));

figure,
subplot(2,1,1)
plot(time,angularVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
hold on
hline = refline([0 meanAngVelocity]);
hline.Color = 'r'; hline.LineStyle = '--';
rline = refline([0 0]);
rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
title('Angular velocity of the fly', 'Interpreter', 'none');
xlabel('Time (s)')
ylabel('Velocity (deg/s)')
legend('Mean angular velocity');

subplot(2,1,2)
histogram(angularVelocity,'FaceColor',[.2 .8 .6])
title('Distribution of angular velocities');
xlabel('Angular velocity (deg/s)');
ylabel('Frequency');

%saveas(gcf,strcat(path,'AngularVelocity_ExpNum', file(11:end-4), '.png'))
%saveas(gcf,strcat(path,'AngularVelocity_ExpNum', file(11:end-4), '.svg'))

%%  Keep the frames during which the fly is moving

% We are going to decide whether a fly is moving or not based on the
% forward velocity. If it's above 1 mm/s we will consider it is moving
% We will work with the downsampled data

downsampled.xPanelPos = downsample(data.xPanelPos,1000/25); %downsample the panels position
dataMoving.xPanelPos = downsampled.xPanelPos(forwardVelocity>1); %keep the position frames during which the fly moved
moving = smoothed.angularPosition(forwardVelocity>1); %keep the angular position frames during which the fly moved

percentageActivity = 100*size(moving)/size(smoothed.angularPosition);
activity = zeros(length(forwardVelocity),1);

for i = 1:length(forwardVelocity)
    if forwardVelocity(i,1) > 1
        activity(i,1) = 1;
    else
        activity(i,1) = 0;
    end
end

figure,
set(gcf, 'Position', [500, 500, 1000, 100])
plot(time,activity,'k');
title('Activity raster plot');
ylabel('Activity');
xlabel('Time (s)');

%saveas(gcf,strcat(path,'ActivityRP_ExpNum', file(11:end-4), '.png'))

%% Output in degrees of the Panels position

% Pos x=92 is 0 deg (ie facing the fly), I measured this empirically

pxToDeg = 360/97; % There are 97 possible positions (the last one = first one) and this represents 360 deg
posToDeg = dataMoving.xPanelPos*pxToDeg;

%% Probability density function of the stimulus position

% Remapping the positions to span -180 to 180 deg
remapPosToDeg = wrapTo180(posToDeg);

% Plot the histogram and probability density
edges = [-180:20:180];
[counts] = histcounts(remapPosToDeg,edges);
probabilities = counts./sum(counts);
degs = linspace(-180,180,length(counts));

figure,
subplot(1,2,1)
h = histogram(remapPosToDeg,edges,'Normalization','probability');
h.FaceColor = [0.2,0.5,1];
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Histogram of the stimulus position');
ylabel('Probability'); xlabel('Stimulus position (deg)');

subplot(1,2,2),
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')

xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Probability density of the stimulus position');
ylabel('Probability density'); xlabel('Stimulus position (deg)');
if contains(typeOfStim,'closed_loop')
   hold on 
   %add line showing the start pos. of the stim.
     if startPos(1)~=1 
        startingPos = (98-startPos(1))*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        startingPos = startPos(1)*pxToDeg; % Correct the offset and multiply by factor to get deg
     end  
    remapStartPos = wrapTo180(startingPos);
    line([remapStartPos remapStartPos],[0 max(probabilities)+0.05],'Color',[1 0 0])

%add rectangle showing where the panel is off (and the bar can't be seen)
noPanel = [71:77]; %xpos where a 2 px bar is partially seen
noPanelDeg = -wrapTo180((noPanel)*pxToDeg);
l1 = line([noPanelDeg(1) noPanelDeg(1)],[0 max(probabilities)+0.05]);
l2 = line([noPanelDeg(7) noPanelDeg(7)],[0 max(probabilities)+0.05]);
set([l1 l2],'Color',[.5 .5 .5]);
patch([noPanelDeg(1) noPanelDeg(7) noPanelDeg(7) noPanelDeg(1)], [0 0 max(probabilities)+0.05 max(probabilities)+0.05],[.5 .5 .5],'FaceAlpha',0.3)
legend('Missing panel','Starting position');

end

%saveas(gcf,strcat(path,'ProbabilityDensityStimPosition_ExpNum', file(11:end-4), '.png'))

%% Polar coordinates analysis of the stimulus position
posToRad = deg2rad(posToDeg);
% some statistics...
CircularStats = circ_stats(posToRad);
[pval,z] = circ_rtest(posToRad);
circLength = circ_r(posToRad,[],2);

%Plot the histogram in polar coordinates

circedges = [0:20:360];
circedges = deg2rad(circedges);

figure,
polarhistogram(posToRad,circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
title('Probability density of the stimulus position');
ax = gca;
d = ax.ThetaDir;
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
hold on
% Add the starting position
points = linspace(0,max(probabilities),1000);
if remapStartPos<0
    starts = repelem(deg2rad(startingPos),1,1000);
else
    starts = repelem(((2*pi)-deg2rad(startingPos)),1,1000);
end
polarplot(starts,points,'r','LineWidth',1.5) %add a line that shows the startPos
% Add the circular mean
circMean = repelem(CircularStats.mean,1,1000);
points2 = linspace(0,max(probabilities)*circLength,1000);
polarplot(circMean,points2,'k','LineWidth',2)
legend('Start position','Circular mean');

%saveas(gcf,strcat(path,'PolarHistogramStimPosition_ExpNum', file(11:end-4), '.png'))


%% Angular position of the stimulus in time

time = linspace(0,(length(rawData)/1000),length(posToDeg));

figure,
plot(time,wrapTo360(posToDeg),'k','HandleVisibility','off')
ylabel('Angular position of the stimulus (deg)'); xlabel('Time (s)');
ylim([0 360]);
title('Angular position of the stimulus as a function of time');;
hline = refline([0 wrapTo360(rad2deg(CircularStats.median))]);
set(hline,'Color',[1,0,0])
legend('Median position');

%saveas(gcf,strcat(path,'AngulaPosStimInTime_ExpNum', file(11:end-4), '.png'))
%% Probability density of the fly heading

flyPosToDegMoving = rad2deg(moving); 
remapFlyPosToDegMoving = wrapTo180(flyPosToDegMoving);

% Plot the histogram and probability density
[countsFlyMoving] = histcounts(remapFlyPosToDegMoving,edges);
probabilitiesFlyMoving = countsFlyMoving./sum(countsFlyMoving);
degsFlyMoving = linspace(-180,180,length(countsFlyMoving));

figure,
subplot(1,2,1)
h = histogram(remapFlyPosToDegMoving,edges,'Normalization','probability')
h.FaceColor = [1,0.2,0.7];
xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
title('Histogram of the fly heading');
ylabel('Probability'); xlabel('Fly heading (deg)');

subplot(1,2,2),
plot(degsFlyMoving,probabilitiesFlyMoving,'k')
%suptitle(typeOfStim)
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
%saveas(gcf,strcat(path,'ProbabilityDensityFlyHeading_ExpNum', file(11:end-4), '.png'))

%% In polar coordinates

FlyPosToRad = deg2rad(flyPosToDegMoving);

CircularStatsFly = circ_stats(FlyPosToRad);

figure,
polarhistogram(FlyPosToRad,circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
title({'Probability density of the fly heading';typeOfStim});
ax = gca;
d = ax.ThetaDir;
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
hold on
points = linspace(0,max(probabilitiesFlyMoving),1000);
starts = repelem(deg2rad(startingPos),1,1000);
polarplot(starts,points,'b','LineWidth',1.5) %add a line that shows the startPos
circMean = repelem(CircularStatsFly.mean,1,1000);
circMeanMinErr = repelem(CircularStatsFly.mean-CircularStatsFly.std,1,1000);
circMeanMaxErr = repelem(CircularStatsFly.mean+CircularStatsFly.std,1,1000);
meanLength = circ_r(FlyPosToRad);
scaledMeanLength = meanLength*max(probabilitiesFlyMoving);
points = linspace(0,scaledMeanLength,1000);
polarplot(circMean,points,'k','LineWidth',1.5)
% polarplot(circMeanMinErr,points,'k','LineWidth',1)
% polarplot(circMeanMaxErr,points,'k','LineWidth',1)
legend('Start position','Circular mean');

%saveas(gcf,strcat(path,'PolarHistFlyHeading_ExpNum', file(11:end-4), '.png'))
%saveas(gcf,strcat(path,'PolarHistFlyHeading_ExpNum', file(11:end-4), '.svg'))


%% Polar analysis taking 60s windows to compute a vector of heading

% Based on Jonathan's 2018 paper, I am taking 60s windows to compute the
% mean of unit vectors, sliding it over the heading time series by 1 s
% increments, and visualizing the distribution in a polar plot.

% I'M NOT SURE IF WHAT I'M DOING IS CORRECT

% circMean60s = movmean(FlyPosToRad,1000);
% figure, subplot(1,2,1)
% polarhistogram(circMean60s,20,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
% title({'Probability density of the fly heading';typeOfStim});
% ax = gca;
% d = ax.ThetaDir;
% ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
% 
% 
% subplot(1,2,2), plot(360-rad2deg(circMean60s),'LineWidth',1.5)
% title('Mean heading over time');
% xlabel('Time');ylabel('Mean heading (deg)');
% ylim([0 360]);

%% Angular position of the fly in time

figure,
plot(time,wrapTo360(flyPosToDegMoving),'k','HandleVisibility','off')
ylabel('Heading angle (deg)'); xlabel('Time (s)');
title('Angular position of the Fly as a function of time');
ylim([0 360]);
hline = refline([0 wrapTo360(rad2deg(CircularStatsFly.median))]);
set(hline,'Color',[1,0,0])
legend('Median heading');

%saveas(gcf,strcat(path,'AngulaPosFlyInTime_ExpNum', file(11:end-4), '.png'))
%saveas(gcf,strcat(path,'AngulaPosFlyInTime_ExpNum', file(11:end-4), '.svg'))

%% Plot 2D virtual trajectory of the fly

dataMoving.Intx = smoothed.Intx(forwardVelocity>1);
dataMoving.Inty = smoothed.Inty(forwardVelocity>1);

figure,
c = linspace(1,(length(data.xPanelPos)/1000),length(smoothed.Intx)); %create a color vector with the time
scatter(smoothed.Intx,smoothed.Inty,[],c)
hold on
plot(smoothed.Intx,smoothed.Inty,'k')
c = colorbar; c.Label.String = 'Time (s)'; %add the colorbar
title('2D trajectory of the fly');
xlabel('x pos (mm)'); ylabel('y pos (mm)');
axis tight equal; %scale the axes with respect to one another

%saveas(gcf,strcat(path,'2DTrajectory_ExpNum', file(11:end-4), '.png'));