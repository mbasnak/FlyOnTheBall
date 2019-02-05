% Experiment 3 analysis

%%% This code analyses the output of the panels and the FicTrac data
%for experiment 3, in which a bar jumps every 200 sec by a predetermined
%magnitude
close all; clear all;

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3'
[file,path] = uigetfile();
load([path,file]);

rawData = daq_data'; %transpose matrix for future use

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;
PanelStatus = 6; %this signal tells whether the panels are on or off.

JumpTime = 200; %how long was the time between bar jumps (in sec)
%% Subset acquisition of x and y pos, as well as FicTrac data

data.xpanelVolts =  rawData (:,xPanels); 
VOLTAGE_RANGE_x = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)

data.yPanelVolts =  rawData (:, yPanels);
VOLTAGE_RANGE_y = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
maxValY = 96;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2

%FicTrac data
data.fictracAngularPosition = rawData ( : , headingFly); 
data.ficTracIntx = rawData ( : , xFly); 
data.ficTracInty = rawData ( : , yFly); 

%% Determine bar jumps
%This is going to be done in two ways:
%1) using the signal from channel 6 to determine when the panels were on
%and off, and adding jumps every 200 sec in between
%2) taking the derivative of the voltage signal from the ypanels channel to identify jumps

% Plot Panel acquisition time
figure,
subplot(3,1,1)
plot(rawData(:,6))
ylabel('Voltage signal (V)');
title('Panels ON and OFF');
xlim([0 size(daq_data,2)]);
ylim([-0.1 10.1]);

% Define when the panels turn on and off
%take the derivative of the panels' signal
panelON = diff(rawData(:,6));
subplot(3,1,2)
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
text(startFrame+35000,5,strcat('startFrame',num2str(startFrame)))
text(endFrame-100000,5,strcat('endFrame',num2str(endFrame)))

% Define the frames with bar jumps
%1) Using the start signal and the trial duration
barjumpFrames = [startFrame:JumpTime*1000:endFrame]; %1000 is the sample rate from the nidaq
barjumpFrames = barjumpFrames(2:end);
barjumpSeconds = barjumpFrames/1000;

%check that the number of trials is the same as requested
if size(barjumpFrames,2)~=TrialNum
    warning('Incorrect trial number');
else
    display('Correct trial number');
end

%plot the data from the yPanels and add lines of the previously determined
%bar jumps
%plot the panels y dimension signal
subplot(3,1,3), plot(data.yPanelVolts)
title('Bar jumps');
xlabel('Time (frames)'); ylabel('Voltage (V)');
hold on
%add the bar jumps
for i = 1:length(barjumpFrames)
    plot([barjumpFrames(i) barjumpFrames(i)],[0 10],'r');
end


%2) Using the signal from the yPanels channel
Jumps = diff(data.yPanelVolts);
Jumps(abs(Jumps)>0.4 & abs(Jumps)<1)=1;
Jumps = round(Jumps);

figure,
suptitle('Bar jumps');
subplot(3,1,1)
plot(data.yPanelVolts)
ylabel('Voltage (V)');xlabel('Time');
subplot(3,1,2)
plot(Jumps);
ylabel('Voltage difference (V)');xlabel('Time');

j = find(Jumps); %indices of the actual bar jumps, taken from the y signal
j = j(2:end); %I leave out the 1st cause it's taking the on signal
jsec = j/1000;

%plot the data from the yPanels and add lines of the previously determined
%bar jumps
%plot the panels y dimension signal
subplot(3,1,3), plot(data.yPanelVolts)
title('Bar jumps');
xlabel('Time (frames)'); ylabel('Voltage (V)');
hold on
%add the bar jumps
for i = 1:length(j)
    plot([j(i) j(i)],[0 10],'r');
end
%% Fixing the data relative to the bar jumps

%The x position of the panels and the FicTrac output are "ignorant" of the
%bar jumps. The x position of the bar will move with the angular position
%of the fly, but the coordinate system changes every time the bar jumps.
%We need to make sure the xpos and heading of the fly are corrected to take
%this coordinate change into account.

yVoltsBJ = data.yPanelVolts(j-1);
yVoltsAJ = data.yPanelVolts(j+1);
ydiff = yVoltsAJ-yVoltsBJ; %this is the offset that I need to adjust by

%Correct the data after every jump except for the last
data.xPanelVoltsUW = data.xpanelVolts;
data.ficTracAngularPositionUW = data.fictracAngularPosition;
for i = 1:size(j)-1
   data.xPanelVoltsUW(j(i)+1:j(i+1)) = (data.xpanelVolts(j(i)+1:j(i+1)))+sum(ydiff(1:i));  
   data.ficTracAngularPositionUW(j(i)+1:j(i+1)) = (data.fictracAngularPosition(j(i)+1:j(i+1)))+sum(ydiff(1:i));
end

%Correct the data after the last jump
data.xPanelVoltsUW(j(end)+1:end) = data.xpanelVolts(j(end)+1:end)+sum(ydiff(1:end));
data.ficTracAngularPositionUW(j(end)+1:end) = data.fictracAngularPosition(j(end)+1:end)+sum(ydiff(1:end));

%I now have to wrap this data to get it to be between 0 and 10 V. 
data.xPanelVolts = data.xPanelVoltsUW;
data.ficTracAngularPosition = data.ficTracAngularPositionUW;

for i = 1:size(data.xPanelVolts)
    
    if data.xPanelVoltsUW(i) > 10
    data.xPanelVolts(i) = data.xPanelVoltsUW(i)-10;
    else
    data.xPanelVolts(i) = data.xPanelVoltsUW(i);
    end
    
    if data.ficTracAngularPositionUW(i) > 10
    data.ficTracAngularPosition(i) = data.ficTracAngularPositionUW(i)-10;
    else
    data.ficTracAngularPosition(i) = data.ficTracAngularPositionUW(i);
    end

end

% Getting the data in x and y position from the voltage info.
data.xPanelPos = round ((data.xPanelVolts  * maxValX ) /VOLTAGE_RANGE_x); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels
data.yPanelPos = round ((data.yPanelVolts  * maxValY) /VOLTAGE_RANGE_y);

%check the wrapping graphically
figure,
subplot(2,1,1),
plot(data.yPanelVolts,'r')
hold on
plot(data.xPanelVoltsUW,'b')
plot(data.xPanelVolts,'g')
ylabel('Voltage (V)');
title('x panels voltage signals before and after wrapping');
ylim([0 max(data.xPanelVoltsUW)]);

subplot(2,1,2),
plot(data.yPanelVolts,'r')
hold on
plot(data.ficTracAngularPositionUW,'b')
plot(data.ficTracAngularPosition,'g')
xlabel('Time (frame)');ylabel('Voltage (V)');
title('angular position voltage signals before and after wrapping');
ylim([0 max(data.ficTracAngularPositionUW)]);

% Getting the degrees that each jump represents and checking they are
% correct

figure,

jumpPos = data.yPanelPos(j+1)-data.yPanelPos(j-1);
degJumps = jumpPos*(360/97);
subplot(1,3,1), plot(degJumps,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from yPanels');
ylabel('deg');xlabel('Trial #');
% compare that to the jump function we have stored
hold on
plot(jumps(1:TrialNum),'b')

%Check if the jump magnitude appears ok in the x panel data
jumpMag = data.xPanelPos(j+1)-data.xPanelPos(j-1); 
degMag = wrapTo180(jumpMag*(360/97));
subplot(1,3,2), plot(degMag,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from xPanels');
ylabel('deg');xlabel('Trial #');
% compare that to the jump function we have stored
hold on
plot(jumps(1:TrialNum),'b')


%check if the jump magnitude appears ok in the angular position data
jumpMag2 = data.ficTracAngularPosition(j+1)-data.ficTracAngularPosition(j-1); 
radMag2 = jumpMag2.* 2 .* pi ./ 10; %convert from voltage to radians
degMag2 = wrapTo180(rad2deg(radMag2)); %convert from radians to degrees and wrap 180
subplot(1,3,3), plot(degMag2,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from angular position');
ylabel('deg');xlabel('Trial #');
% compare that to the jump function we have stored
hold on
plot(jumps(1:TrialNum),'b')

%% Downsample, unwrap and smooth position data, then get velocity and smooth

[smoothed] = singleTrialVelocityAnalysis(data,1000);
% for most experiments I have used 1000 Hz as a sample rate, and it is what
% I will use from now on, so that's how I'll leave it, but this could be
% changed in case of need


%% global analysis
%% Forward velocity analysis

% The forward velocity is a good indicative of whether the fly is walking
% well in that trial or not

%To look at the velocity, I should do so without adding the offsets of the
%jumps, because the offsets will give me weird jumps in velocity that the
%fly might not actually have made. For the forward velocity, I can use the
%one calculates by the singleTrialVelocityAnalysis function, because it
%uses the variable Intx, that I haven't corrected for the offset

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
for i = 1:length(jsec)
     plot([jsec(i) jsec(i)],[min(forwardVelocity)-10 max(forwardVelocity)+10],'g');
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

%saveas(gcf,strcat(path,'ForwardVelocity_ExpNum', file(11:end-4), '.png'))

%% Angular velocity

%To look at the velocity, I should do so without adding the offsets of the
%jumps, because the offsets will give me weird jumps in velocity that the
%fly might not actually have made. For the angular velocity, I'm using a new 
%function that I made to downsample, unwrap, smooth and get the velocity
%with the angular position uncorrected for the offset

AngularPosition = rawData ( : , headingFly);
angularVelocity = getAngVel(AngularPosition);

meanAngVelocity = mean(angularVelocity);
time = linspace(0,(length(rawData)/1000),length(angularVelocity));

figure,
subplot(2,1,1)
plot(time,angularVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
ylim([min(angularVelocity)-50 max(angularVelocity)+50]);
hold on
for i = 1:length(j) %add the bar jumps
     plot([jsec(i) jsec(i)],[min(angularVelocity)-50 max(angularVelocity)+50],'g');
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

%saveas(gcf,strcat(path,'AngularVelocity_ExpNum', file(11:end-4), '.png'))


%%  Activity levels

% We are going to decide whether a fly is moving or not based on the
% forward velocity. If it's above 1 mm/s we will consider it is moving
% We will work with the downsampled data

downsampled.xPanelPos = downsample(data.xPanelPos,1000/25); %downsample the panels position
dataMoving.xPanelPos = downsampled.xPanelPos(forwardVelocity>1); %keep the position frames during which the fly moved
moving = smoothed.angularPosition(forwardVelocity>1); %keep the angular position frames during which the fly moved
% I need to think more carefully about whether there is a problem in using
% the forward velocity obtained from smoothed data to choose the frames
% during which the fly is or not moving in unsmoothed data...

percentageActivity = 100*size(moving)/size(smoothed.angularPosition);
activity = zeros(length(forwardVelocity),1);

for i = 1:length(forwardVelocity)
    if forwardVelocity(i,1) > 1
        activity(i,1) = 1;
    else
        activity(i,1) = 0;
    end
end

time = linspace(0,(length(rawData)/1000),length(activity));

figure,
set(gcf, 'Position', [500, 500, 1000, 100])
plot(time,activity,'k');
title('Activity raster plot');
ylabel('Activity');
xlabel('Time (s)');
xlim([0 time(end)]);

%saveas(gcf,strcat(path,'ActivityRP_ExpNum', file(11:end-4), '.png'))


%% Fly's heading thoughout the experiment

%I think for the fly's heading I don't need to remap anything, cause it
%should be based on fictrac's calibration and have the front be 0?

figure('Position', [100 100 1200 900]),

% Heading in time, with every frame
flyPos180 = wrapTo180(rad2deg(smoothed.angularPosition));
time = linspace(0,(length(rawData)/1000),length(flyPos180));
subplot(3,2,1)
scatter(flyPos180, time, [], forwardVelocity); %add the forward velocity as a color
hold on
plot(flyPos180, time,'k','HandleVisibility','off')
xlabel('Heading angle (deg)'); ylabel('Time (s)');
title('Angular position of the Fly as a function of time');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [jsec(i) jsec(i)],'r--');
end
ax = gca;
ax.YDir = 'reverse';

 % Heading in time, with only moving frames.
%time = linspace(0,(length(rawData)/1000),length(flyPos180(forwardVelocity>1)));
subplot(3,2,2)
scatter(flyPos180(forwardVelocity>1), time(forwardVelocity>1));
xlabel('Heading angle (deg)'); ylabel('Time (s)');
title('Angular position of the Fly as a function of time');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [jsec(i) jsec(i)],'r--');
end
ax = gca;
ax.YDir = 'reverse'; 
 

% Plot the histogram and probability density
edges = [-180:20:180];
[countsFly] = histcounts(flyPos180,edges);
[countsFlyMoving] = histcounts(flyPos180(forwardVelocity>1),edges);
probabilitiesFly = countsFly./sum(countsFly);
probabilitiesFlyMoving = countsFlyMoving./sum(countsFlyMoving);
degsFly = linspace(-180,180,length(countsFly));
degsFlyMoving = linspace(-180,180,length(countsFlyMoving));

subplot(3,2,3) %with every frame
plot(degsFly,probabilitiesFly,'k')
xlim([-180 180]); ylim([0 max(probabilitiesFly)+0.05]);
title('Probability density of the fly heading with every frame');
ylabel('Probability density'); xlabel('Fly heading (deg)');
%hold on
%patch([noPanelDeg(1) noPanelDeg(7) noPanelDeg(7) noPanelDeg(1)], [0 0 max(probabilitiesFlyMoving)+0.05 max(probabilitiesFlyMoving)+0.05],[.5 .5 .5],'FaceAlpha',0.3)
subplot(3,2,4) %with moving frames only
plot(degsFlyMoving,probabilitiesFlyMoving,'k')
xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
title('Probability density of the fly heading with moving frames');
ylabel('Probability density'); xlabel('Fly heading (deg)');

% In polar coordinates...
%Taking every frame into account
CircularStatsFly = circ_stats(smoothed.angularPosition);
circedges = [0:20:360];
circedges = deg2rad(circedges);
subplot(3,2,5)
polarhistogram(smoothed.angularPosition,circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
ax = gca;
d = ax.ThetaDir;
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%with only moving frames
CircularStatsFlyMoving = circ_stats(smoothed.angularPosition(forwardVelocity>1));
subplot(3,2,6)
polarhistogram(smoothed.angularPosition(forwardVelocity>1),circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
ax = gca;
d = ax.ThetaDir;
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top


%% Look at the goal and calculate the distance to it...


%Taking all the experiment
figure,
%with every frame
goal = mean(wrapTo180(rad2deg(smoothed.angularPosition)));
dist2goal = flyPos180-goal;
[countsDist] = histcounts(dist2goal,edges);
probabilitiesDist = countsDist./sum(countsDist);
degsFlyDist = linspace(-180,180,length(countsDist));
subplot(1,2,1), plot(degsFlyDist,probabilitiesDist,'r')
title('Distance to the goal with every frame');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

%taking only 'moving frames'
goalMoving = mean(wrapTo180(rad2deg(smoothed.angularPosition(forwardVelocity>1))));
dist2goalMoving = flyPos180(forwardVelocity>1)-goalMoving;
[countsDistMoving] = histcounts(dist2goalMoving,edges);
probabilitiesDistMoving = countsDistMoving./sum(countsDistMoving);
degsFlyDistMoving = linspace(-180,180,length(countsDistMoving));
subplot(1,2,2), plot(degsFlyDistMoving,probabilitiesDistMoving,'r')
title('Distance to the goal with moving frames');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');


% Moving window for the goal


%% 2D trajectories

figure,
c = linspace(1,(length(data.xPanelPos)/1000),length(smoothed.Intx)); %create a color vector with the time
scatter(smoothed.Intx,smoothed.Inty,[],c)
hold on
plot(smoothed.Intx,smoothed.Inty,'k')
c = colorbar; c.Label.String = 'Time (s)'; %add the colorbar
title('2D trajectory of the fly');
xlabel('x pos (mm)'); ylabel('y pos (mm)');
axis tight equal; %scale the axes with respect to one another



%% Per 'trial'

%I'm using a function to get and smooth data around the jumps
sec = 10; %how many sec before and after the jump I want to look at
perTrialData = getDataAroundJump(rawData,j,sec);

%apend jump data and save for using in the pooled flies analysis
perTrialData.jumpMag = jumps;

save(strcat(path,'perTrialData.mat'),'perTrialData');

%% Velocity around the jumps

time = linspace(-sec,sec,length(perTrialData.forwardVel));

%Individually

for i = 1:length(j)
    figure,
    subplot(1,2,1)
    plot(time,perTrialData.forwardVel(:,i),'.')
    hold on
    plot(time,perTrialData.forwardVel(:,i))
    title('Forward velocity around the bar jumps');
    xlabel('Time(s)');
    ylabel('Velocity (mm/s)');
    
    subplot(1,2,2)
    plot(time,perTrialData.angVel(:,i),'.')
    hold on
    plot(time,perTrialData.angVel(:,i))
    title('Angular velocity around the bar jumps');
    xlabel('Time(s)');
    ylabel('Velocity (mm/s)');
end


%Pooling results from similar magnitude jumps
%1) make a jump vector using the preloaded jump vector from the experiment
%and taking only as many elements as trials there were
trials = jumps(1:TrialNum);
%2) identify elements in that vector belonging to the 4 different groups,
%and put the perTrialData into those groups
Data45.forwardVel = perTrialData.forwardVel(:,trials == 45);
Data45.angVel = perTrialData.angVel(:,trials == 45);
DataNeg45.forwardVel = perTrialData.forwardVel(:,trials == -45);
DataNeg45.angVel = perTrialData.angVel(:,trials == -45);
Data90.forwardVel = perTrialData.forwardVel(:,trials == 90);
Data90.angVel = perTrialData.angVel(:,trials == 90);
DataNeg90.forwardVel = perTrialData.forwardVel(:,trials == -90);
DataNeg90.angVel = perTrialData.angVel(:,trials == -90);

%plot forward and angular velocity for every group
figure,
subplot(1,2,1)
plot(time,Data45.forwardVel,'.')
hold on
plot(time,Data45.forwardVel)
title('Forward velocity for 45 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');    
subplot(1,2,2)
plot(time,Data45.angVel,'.')
hold on
plot(time,Data45.angVel)
title('Angular velocity for 45 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');

figure,
subplot(1,2,1)
plot(time,DataNeg45.forwardVel,'.')
hold on
plot(time,DataNeg45.forwardVel)
title('Forward velocity for -45 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');    
subplot(1,2,2)
plot(time,DataNeg45.angVel,'.')
hold on
plot(time,DataNeg45.angVel)
title('Angular velocity for -45 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');


figure,
subplot(1,2,1)
plot(time,Data90.forwardVel,'.')
hold on
plot(time,Data90.forwardVel)
title('Forward velocity for 90 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');    
subplot(1,2,2)
plot(time,Data90.angVel,'.')
hold on
plot(time,Data90.angVel)
title('Angular velocity for 90 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');


figure,
subplot(1,2,1)
plot(time,DataNeg90.forwardVel,'.')
hold on
plot(time,DataNeg90.forwardVel)
title('Forward velocity for -90 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');    
subplot(1,2,2)
plot(time,DataNeg90.angVel,'.')
hold on
plot(time,DataNeg90.angVel)
title('Angular velocity for -90 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');


%plot mean forward and angular velocity per group
meanForwardVel45 = mean(Data45.forwardVel,2);
meanForwardVelNeg45 = mean(DataNeg45.forwardVel,2);
meanForwardVel90 = mean(Data90.forwardVel,2);
meanForwardVelNeg90 = mean(DataNeg90.forwardVel,2);

meanAngVel45 = mean(Data45.angVel,2);
meanAngVelNeg45 = mean(DataNeg45.angVel,2);
meanAngVel90 = mean(Data90.angVel,2);
meanAngVelNeg90 = mean(DataNeg90.angVel,2);

figure,
subplot(1,2,1)
plot(time,meanForwardVel45)
hold on
plot(time,meanForwardVelNeg45)
plot(time,meanForwardVel90)
plot(time,meanForwardVelNeg90)
title('Mean forward velocity');
legend({'45','-45','90','-90'});
ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
subplot(1,2,2)
plot(time,meanAngVel45)
hold on
plot(time,meanAngVelNeg45)
plot(time,meanAngVel90)
plot(time,meanAngVelNeg90)
title('Mean angular velocity');
legend({'45','-45','90','-90'});
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');



%% Output in degrees of the Panels position

% Pos x=73 is 0 deg (ie facing the fly), I measured this empirically

pxToDeg = 360/96; % There are 97 possible positions (the last one = first one) and this represents 360 deg
posToDeg = zeros(1,length(data.xPanelPos));

% Convert from xpos to degrees, knowing that xpos 70 = 0 deg
for i=1:length(data.xPanelPos)
    if data.xPanelPos(i) == 70
        posToDeg(i) = 0;
    elseif data.xPanelPos(i) >70 
        posToDeg(i) = (data.xPanelPos(i)-70)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        posToDeg(i) = (data.xPanelPos(i)+27)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end

% Probability density function of the stimulus position throughout the experiment

% Remapping the positions to span -180 to 180 deg
remapPosToDeg = wrapTo180(posToDeg);

% Plot the histogram and probability density
edges = [-180:20:180];
[counts] = histcounts(remapPosToDeg,edges);
probabilities = counts./sum(counts);
degs = linspace(-180,180,length(counts));

figure ('Position', [100 100 1200 900]),
subplot(1,3,1),
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Probability density of the stimulus position');
ylabel('Probability density'); xlabel('Stimulus position (deg)');

% Global Polar coordinates analysis of the stimulus position
posToRad = deg2rad(posToDeg);
% some statistics...
CircularStats = circ_stats(posToRad);
[pval,z] = circ_rtest(posToRad);
%Plot the histogram in polar coordinates
circedges = [0:20:360];
circedges = deg2rad(circedges);
subplot(1,3,2),
polarhistogram(posToRad,circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
title('Probability density of the stimulus position');
ax = gca;
d = ax.ThetaDir;
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

% Angular position of the stimulus in time

time = linspace(0,(length(rawData)/1000),length(posToDeg));
subplot(1,3,3)
plot(fliplr(wrapTo180(posToDeg)),fliplr(time),'k','HandleVisibility','off')
xlabel('Angular position of the stimulus (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
title('Angular position of the stimulus as a function of time');
hold on 
for i = 1:length(j)
     plot([jsec(i) jsec(i)],[0 360],'g');
end
ax = gca;
ax.YDir = 'reverse';

