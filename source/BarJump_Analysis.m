% Bar jump experiment analysis

%%% This code analyses the output of the panels and the FicTrac data
close all; clear all;

% prompt the user to select the file to open and load it.
[file,path] = uigetfile();
load([path,file]);

rawData = daq_data';

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;
PanelStatus = 6; %this signal tells whether the panels are on or off.

JumpTime = 200; %how long was the time between bar jumps (in sec)
%% Subset acquisition of x and y pos, as well as FicTrac data

data.xPanelVolts =  rawData (:,xPanels); 
VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
data.xPanelPos = round ((data.xPanelVolts  * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels

data.yPanelVolts =  rawData (:, yPanels);
VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
maxValY = 96;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2
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

%% Determine bar jumps

% Plot Panel acquisition time
figure,
subplot(2,1,1)
plot(rawData(:,6))
ylabel('Voltage signal (V)');
title('Panels ON and OFF');
xlim([0 size(daq_data,2)]);
ylim([-0.1 10.1]);

% Define when the panels turn on and off
panelON = diff(rawData(:,6));
subplot(2,1,2)
plot(panelON)
xlim([0 size(daq_data,2)]);
xlabel('Time (frames)');
ylabel('Diff of voltage signal (V)');
%Find the frame when the panels turn on
[M, I] = max(panelON);
%Find the frame when they turn off.
[M2, I2] = min(panelON);
startFrame = I(1);
endFrame = I2(1);
%add text indicating start and end frames
text(startFrame+35000,5,strcat('startFrame_',num2str(startFrame)))
text(endFrame-600000,5,strcat('endFrame_',num2str(endFrame)))

% Define the frames with bar jumps
%1) Using the start signal and the trial duration
barjumpFrames = [startFrame:JumpTime*1000:endFrame];
barjumpSeconds = barjumpFrames/1000;

%check that the number of trials is the same as requested
if size(barjumpFrames,2)~=TrialNum
    warning('Incorrect trial number');
else
    display('Correct trial number');
end

%plot the data from the yPanels and add lines of the previously determined
%bar jumps
figure, plot(data.yPanelVolts)
title('Bar jumps');
xlabel('Time (frames)'); ylabel('Voltage (V)');
hold on
for i = 1:length(barjumpFrames)
    plot([barjumpFrames(i) barjumpFrames(i)],[0 10],'r');
end

%2) Using the signal from the yPanels channel
% jumps = diff(data.yPanelVolts);
% jumps(abs(jumps)>0.08 & abs(jumps)<1)=1;
% jumps = round(jumps);
% 
% figure,
% suptitle('Bar jumps');
% subplot(1,2,1)
% plot(data.yPanelVolts)
% ylabel('Voltage (V)');xlabel('Time');
% subplot(1,2,2)
% plot(jumps);
% ylabel('Voltage difference (V)');xlabel('Time');
% 

% j = find(jumps); %indices of the actual bar jumps, taken from the y signal
% j = j(2:end);
% jsec = j/1000;

%% Forward velocity analysis

% The forward velocity is a good indicative of whether the fly is walking
% well in that trial or not

forwardVelocity = smoothed.xVel;
meanVelocity = mean(forwardVelocity);
time = linspace(0,(length(rawData)/1000),length(forwardVelocity));

%Global forward velocity analysis
figure,
subplot(2,1,1)
plot(time,forwardVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
ylim([min(forwardVelocity)-10 max(forwardVelocity)+10]);
hold on
for i = 1:length(barjumpSeconds)
     plot([barjumpSeconds(i) barjumpSeconds(i)],[min(forwardVelocity)-10 max(forwardVelocity)+10],'g');
end
hline = refline([0 meanVelocity]);
hline.Color = 'r'; hline.LineStyle = '--';
rline = refline([0 0]);
rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
title('Forward velocity of the fly', 'Interpreter', 'none');
xlabel('Time (s)')
ylabel('Velocity (mm/s)')

subplot(2,1,2)
histogram(forwardVelocity,'FaceColor',[.4 .2 .6])
title('Distribution of forward velocities');
xlabel('Forward velocity (mm/s)');
ylabel('Frequency');

saveas(gcf,strcat(path,'ForwardVelocity_ExpNum', file(11:end-4), '.png'))

%% Velocity around the jumps

for i = 1:length(barjumpSeconds)-1
    jumpVel(:,i) = forwardVelocity(barjumpSeconds(i+1)-50:barjumpSeconds(i+1)+50);
end

time = linspace(-2,2,length(jumpVel));

figure, plot(time,jumpVel,'.')
hold on
plot(time,jumpVel)
title('Forward velocity around the bar jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');
ylim([min(min(jumpVel))-10 max(max(jumpVel))+10]);
line([0 0], [min(min(jumpVel))-10 max(max(jumpVel))+10],'Color','black');

%It seems like there is a change in velocity about 0.25 sec before the
%jump,but I think that time difference might be accounted for by the smoothing,
%and gradient

saveas(gcf,strcat(path,'AroundJumpVelocity_ExpNum', file(11:end-4), '.png'))


% Determine the extent of the bar jumps
%I'm not sure that this can be done with the data I have now.


% Looking at them individually
for i = 1:length(barjumpSeconds)-1
   figure,
   plot(time,jumpVel(:,i),'.')
   hold on
   plot(time,jumpVel(:,i))
   title('Forward velocity around the bar jumps');
   xlabel('Time(s)');
   ylabel('Velocity (mm/s)');
   ylim([min(min(jumpVel))-10 max(max(jumpVel))+10]);
   line([0 0], [min(min(jumpVel))-10 max(max(jumpVel))+10],'Color','black');     
end



%% Angular velocity analysis

angularVelocity = smoothed.angularVel;
meanAngVelocity = mean(angularVelocity);
time = linspace(0,(length(rawData)/1000),length(angularVelocity));

figure,
subplot(2,1,1)
plot(time,angularVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
ylim([min(angularVelocity)-50 max(angularVelocity)+50]);
hold on
for i = 1:length(barjumpSeconds)
     plot([barjumpSeconds(i) barjumpSeconds(i)],[min(angularVelocity)-50 max(angularVelocity)+50],'g');
end
hline = refline([0 meanAngVelocity]);
hline.Color = 'r'; hline.LineStyle = '--';
rline = refline([0 0]);
rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
title('Angular velocity of the fly', 'Interpreter', 'none');
xlabel('Time (s)')
ylabel('Velocity (deg/s)')

subplot(2,1,2)
histogram(angularVelocity,'FaceColor',[.2 .8 .6])
title('Distribution of angular velocities');
xlabel('Angular velocity (deg/s)');
ylabel('Frequency');

saveas(gcf,strcat(path,'AngularVelocity_ExpNum', file(11:end-4), '.png'))


%%  Keep the frames during which the fly is moving

% We are going to decide whether a fly is moving or not based on the
% forward velocity. If it's above 0.7 mm/s we will consider it is moving
% We will work with the downsampled data

downsampled.xPanelPos = downsample(data.xPanelPos,1000/25); %downsample the panels position
dataMoving.xPanelPos = downsampled.xPanelPos(forwardVelocity>1.5); %keep the position frames during which the fly moved
moving = smoothed.angularPosition(forwardVelocity>1.5); %keep the angular position frames during which the fly moved

percentageActivity = 100*size(moving)/size(smoothed.angularPosition);
activity = zeros(length(forwardVelocity),1);

for i = 1:length(forwardVelocity)
    if forwardVelocity(i,1) > 1.5
        activity(i,1) = 1;
    else
        activity(i,1) = 0;
    end
end

figure,
set(gcf, 'Position', [500, 500, 1000, 100])
plot(time,activity,'k');
xlim([0 time(end)]);
title('Activity raster plot');
ylabel('Activity');
xlabel('Time (s)');

saveas(gcf,strcat(path,'ActivityRP_ExpNum', file(11:end-4), '.png'))

%% Output in degrees of the Panels position

% Pos x=92 is 0 deg (ie facing the fly), I measured this empirically

pxToDeg = 360/97; % There are 97 possible positions (the last one = first one) and this represents 360 deg
posToDeg = zeros(1,length(dataMoving.xPanelPos));

% Convert from xpos to degrees, knowing that xpos 92 = 0 deg
for i=1:length(dataMoving.xPanelPos)
    if dataMoving.xPanelPos(i) ==93 | dataMoving.xPanelPos(i) ==94 | dataMoving.xPanelPos(i) ==95 | dataMoving.xPanelPos(i) ==96 | dataMoving.xPanelPos(i) ==97
        posToDeg(i) = (dataMoving.xPanelPos(i)-92)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        posToDeg(i) = (dataMoving.xPanelPos(i)+4)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end


%% Probability density function of the stimulus position throughout the experiment

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
hold on 
%add rectangle showing where the panel is off (and the bar can't be seen)
noPanel = [17:23]; %xpos where a 2 px bar can't be seen.
noPanelDeg = (noPanel+4)*pxToDeg;
l1 = line([noPanelDeg(1) noPanelDeg(1)],[0 max(probabilities)+0.05]);
l2 = line([noPanelDeg(7) noPanelDeg(7)],[0 max(probabilities)+0.05]);
set([l1 l2],'Color',[.5 .5 .5]);
patch([noPanelDeg(1) noPanelDeg(7) noPanelDeg(7) noPanelDeg(1)], [0 0 max(probabilities)+0.05 max(probabilities)+0.05],[.5 .5 .5],'FaceAlpha',0.3)

%saveas(gcf,strcat(path,'ProbabilityDensityStimPosition_ExpNum', file(11:end-4), '.png'))

% I don't think this is useful, given that I think the x position is being
% taken without taking into account the jumps.


%% Global Polar coordinates analysis of the stimulus position
posToRad = deg2rad(posToDeg);
% some statistics...
CircularStats = circ_stats(posToRad);
[pval,z] = circ_rtest(posToRad);

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
% Add the circular mean
circMean = repelem(CircularStats.mean,1,1000);
points = linspace(0,max(probabilities),1000);
polarplot(circMean,points,'k','LineWidth',1.5)
legend('Circular mean');

%saveas(gcf,strcat(path,'PolarHistogramStimPosition_ExpNum', file(11:end-4), '.png'))


%% Angular position of the stimulus in time

time = linspace(0,(length(rawData)/1000),length(posToDeg));

figure,
plot(time,wrapTo360(posToDeg),'k','HandleVisibility','off')
ylabel('Angular position of the stimulus (deg)'); xlabel('Time (s)');
ylim([0 360]);xlim([0 max(time)]);
title('Angular position of the stimulus as a function of time');;
hline = refline([0 wrapTo360(rad2deg(CircularStats.median))]);
set(hline,'Color',[1,0,0])
hold on 
for i = 1:length(barjumpSeconds)
     plot([barjumpSeconds(i) barjumpSeconds(i)],[min(angularVelocity) max(angularVelocity)],'g');
end
%legend('Median position');

%this is weird. The jumps don't seem to show in the data.
%It might be because the data is smoothed at this point.

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
xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
title('Probability density of the fly heading');
ylabel('Probability density'); xlabel('Fly heading (deg)');
%Add the starting pos of the bar if the stimulus was a closed-loop bar
hold on
patch([noPanelDeg(1) noPanelDeg(7) noPanelDeg(7) noPanelDeg(1)], [0 0 max(probabilitiesFlyMoving)+0.05 max(probabilitiesFlyMoving)+0.05],[.5 .5 .5],'FaceAlpha',0.3)

%saveas(gcf,strcat(path,'ProbabilityDensityFlyHeading_ExpNum', file(11:end-4), '.png'))

FlyPosToRad = deg2rad(flyPosToDegMoving);

CircularStatsFly = circ_stats(FlyPosToRad);

figure,
polarhistogram(FlyPosToRad,circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
title('Probability density of the fly heading');
ax = gca;
d = ax.ThetaDir;
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
hold on
circMean = repelem(CircularStatsFly.mean,1,1000);
polarplot(circMean,points,'k','LineWidth',1.5)
legend('Circular mean');

saveas(gcf,strcat(path,'PolarHistFlyHeading_ExpNum', file(11:end-4), '.png'))

meanLength = circ_r(FlyPosToRad);


%% Angular position of the fly in time

time = linspace(0,(length(rawData)/1000),length(posToDeg));
angMoving = wrapTo360(flyPosToDegMoving);

figure,
plot(time,angMoving,'k','HandleVisibility','off')
ylabel('Heading angle (deg)'); xlabel('Time (s)');
title('Angular position of the Fly as a function of time');
ylim([0 360]); xlim([0 max(time)]);
hline = refline([0 wrapTo360(rad2deg(CircularStatsFly.median))]);
set(hline,'Color',[1,0,0])
hold on
for i = 1:length(barjumpSeconds)
     plot([barjumpSeconds(i) barjumpSeconds(i)],[min(angularVelocity) max(angularVelocity)],'g');
end
%legend('Median heading');


% Angular position around the bar jumps

for i = 1:length(barjumpSeconds)-1
    jumpAngMoving(:,i) = angMoving(barjumpFrames(i+1)/70-5000/70:barjumpFrames(i+1)/70+5000/70);
end

% This is a bit weird...

time = linspace(-200,200,length(jumpAngMoving));

figure, plot(time,jumpAngMoving,'.')
hold on
plot(time,jumpAngMoving)
title('Angular position of the fly around the bar jumps');
xlabel('Time(s)');
ylabel('Heading angle (deg)');
ylim([0 360]);
line([0 0], [0 360],'Color','black');

%saveas(gcf,strcat(path,'AngulaPosFlyInTime_ExpNum', file(11:end-4), '.png'))


%% Plot 2D virtual trajectory of the fly during the full experiment

dataMoving.Intx = smoothed.Intx(forwardVelocity>1.5);
dataMoving.Inty = smoothed.Inty(forwardVelocity>1.5);

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

%% Per-"trial" analysis

% Separte the data by "trial"

for i = 1:TrialNum-1  
    trials(i).xPanelVolts = data.xPanelVolts(barjumpFrames(i):barjumpFrames(i+1)-1);
    trials(i).yPanelVolts = data.yPanelVolts(barjumpFrames(i):barjumpFrames(i+1)-1);
    trials(i).xPanelPos = data.xPanelPos(barjumpFrames(i):barjumpFrames(i+1)-1);
    trials(i).yPanelPos = data.yPanelPos(barjumpFrames(i):barjumpFrames(i+1)-1);
    trials(i).ficTracAngularPosition = data.ficTracAngularPosition(barjumpFrames(i):barjumpFrames(i+1)-1);
    trials(i).ficTracIntx = data.ficTracIntx(barjumpFrames(i):barjumpFrames(i+1)-1);
    trials(i).ficTracInty = data.ficTracInty(barjumpFrames(i):barjumpFrames(i+1)-1);
end

%add last trial

trials(TrialNum).xPanelVolts = data.xPanelVolts(barjumpFrames(end):I2);
trials(TrialNum).yPanelVolts = data.yPanelVolts(barjumpFrames(end):I2);
trials(TrialNum).xPanelPos = data.xPanelPos(barjumpFrames(end):I2);
trials(TrialNum).yPanelPos = data.yPanelPos(barjumpFrames(end):I2);
trials(TrialNum).ficTracAngularPosition = data.ficTracAngularPosition(barjumpFrames(end):I2);
trials(TrialNum).ficTracIntx = data.ficTracIntx(barjumpFrames(end):I2);
trials(TrialNum).ficTracInty = data.ficTracInty(barjumpFrames(end):I2);

% downsample, smooth as subset the moving frames
for i = 1:size(trials,2)
    [smoothedTrials(i)] = singleTrialVelocityAnalysis(trials(i),1000);
    forwardVelocityTrials{i} = smoothedTrials(i).xVel;
    downsampledTrials(i).xPanelPos = downsample(trials(i).xPanelPos,1000/25); %downsample the panels position
    dataMovingTrials(i).xPanelPos = downsampledTrials(i).xPanelPos(forwardVelocityTrials{i}>1.5); %keep the position frames during which the fly moved
    movingTrials{i} = smoothedTrials(i).angularPosition(forwardVelocityTrials{i}>1.5); %keep the angular position frames during which the fly moved
end

%% Determine activity percentage and plots

figure,
suptitle('Activity raster plots');
xlabel('Time');
ylabel('Activity');
for i = 1:size(trials,2)
percentageActivity(i) = 100*size(movingTrials{i})/size(smoothedTrials(i).angularPosition);
activityTrials{i} = zeros(length(forwardVelocityTrials{i}),1);

    for p = 1:length(forwardVelocityTrials{i})
        if forwardVelocityTrials{i}(p) > 1.5
        activityTrials{i}(p) = 1;
        else
        activityTrials{i}(p) = 0;
        end
    end

subplot(5,TrialNum/5,i)
plot(activityTrials{i},'k');

end

% Weibull fit

% %1) Determine the first frame with activity per trial:
% 
% Act = cellfun(@find, activityTrials, 'UniformOutput', 0);
% 
% for i = 1:length(Act)
%     firstAct(i) = Act{i}(1);
% end
% 
% %2) Calculate the probabilities of not being active for a certain time
% 
% frames = [0:300];
% 
% for i = 1:length(frames)
%     probInactive(i) = sum(firstAct>=frames(i))/size(firstAct,2);
% end
% 
% figure, plot(probInactive)
% xlabel('Time (frames)'); ylabel('Probability of IBI >= time');
% title('Survival analysis of IBI');

%It looks bad with soo few trials! This is probably better suited for a
%global analysis...

%3) Fit a Weibull distribution
% 
% %keep the probabilities above zero to be able to perform the fit
% probInactivePositive = probInactive(probInactive>0);
% parmhat = wblfit(probInactivePositive);
% 
% y = exp(-(probInactivePositive/parmhat(1)).^parmhat(2));
%% Determine fly heading per trial

figure,
%suptitle('Fly heading angle');

for i = 1:size(trials,2)
    flyPosToDegMovingTrials{i} = rad2deg(movingTrials{i}); 
    remapFlyPosToDegMovingTrials{i} = wrapTo180(flyPosToDegMovingTrials{i});

% Plot the histogram and probability density
    [countsFlyMovingTrials{i}] = histcounts(remapFlyPosToDegMovingTrials{i},edges);
    probabilitiesFlyMovingTrials{i} = countsFlyMovingTrials{i}./sum(countsFlyMovingTrials{i});
    degsFlyMovingTrials{i} = linspace(-180,180,length(countsFlyMovingTrials{i}));

    subplot(5,TrialNum/5,i)

    plot(degsFlyMovingTrials{i},probabilitiesFlyMovingTrials{i},'k')
    xlim([-180 180]);
    ylim([0 max(probabilitiesFlyMovingTrials{i})+0.05]);

%Add the starting pos of the bar 
    hold on
    startPos(i) =  data.yPanelPos(barjumpFrames(i)+1);
     if startPos(i) ==93 | startPos(i) ==94 | startPos(i) ==95 | startPos(i) ==96 | startPos(i) ==97
        startingPos(i) = (startPos(i)-92)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        startingPos(i) = (startPos(i)+4)*pxToDeg; % Correct the offset and multiply by factor to get deg
     end  
    remapStartPos(i) = wrapTo180(startingPos(i));
    line([remapStartPos(i) remapStartPos(i)],[0 max(probabilitiesFlyMovingTrials{i})+0.05],'Color',[1 0 0])

end

suptitle('Heading of the fly');
%xlabel('Heading angle (deg)');
%ylabel('Density');

%saveas(gcf,strcat(path,'ProbabilityDensityFlyHeading_ExpNum', file(11:end-4), '.png'))

%% Plot the heading in circular coodinates

figure,
suptitle('Probability density of the fly heading');

for i = 1:size(trials,2)
    
    FlyPosToRadTrials{i} = deg2rad(flyPosToDegMovingTrials{i});
    CircularStatsFlyTrials(i) = circ_stats(FlyPosToRadTrials{i});

    subplot(TrialNum/5,5,i)
    polarhistogram(FlyPosToRadTrials{i},circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
    ax = gca;
    d = ax.ThetaDir;
    ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

    hold on
    
    pointsTrials{i} = linspace(0,max(probabilitiesFlyMovingTrials{i}),1000);
    starts{i} = repelem(deg2rad(startingPos(i)),1,1000);
    polarplot(starts{i},pointsTrials{i},'r','LineWidth',1.5) %add a line that shows the startPos
    %circMeanTrials{i} = repelem(CircularStatsFlyTrials(i).mean,1,1000);
    %polarplot(circMeanTrials{i},pointsTrials{i},'k','LineWidth',1.5)
    %legend('Start position','Circular mean');
    
    
    Ax = gca; % current axes
    Ax.ThetaGrid = 'on';
    Ax.RGrid = 'on';
    Ax.RTickLabel = []; 
    Ax.ThetaTickLabel = [];
end

%saveas(gcf,strcat(path,'PolarHistFlyHeading_ExpNum', file(11:end-4), '.png'))
%% Plot the heading angle as a function of the trial number


for i = 1:size(trials,2)
    circMeanTrials(i) = CircularStatsFlyTrials(i).mean;
    circStdTrials(i) = CircularStatsFlyTrials(i).std;
end

figure,
bar(rad2deg(circMeanTrials))
hold on
errorbar(rad2deg(circMeanTrials),rad2deg(circStdTrials),'.')
ylim([-250 250]);
ylabel('Angle (deg)'); xlabel('Trial #');
title('Mean heading across trials');

saveas(gcf,strcat(path,'MeanHeadingAcrossTrials_ExpNum', file(11:end-4), '.png'))
