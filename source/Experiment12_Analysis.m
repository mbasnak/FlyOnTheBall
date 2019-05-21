% Experiment 12 analysis

%%% This code analyses the output of the panels and the FicTrac data
%for experiment 12, 
close all; clear all;

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment12'
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


%% Determine trial change

% Plot Panel acquisition time
figure,
subplot(2,1,1)
plot(rawData(:,6))
ylabel('Voltage signal (V)');
title('Panels signal');
xlim([0 size(daq_data,2)]);
ylim([-0.1 10.1]);

% Define when the panels turn on and off
%take the derivative of the panels' signal
panelON = diff(rawData(:,6));
subplot(2,1,2)
plot(panelON)
xlim([0 size(daq_data,2)]);
xlabel('Time (frames)');
ylabel('Diff of voltage signal (V)');


% get the trial changes using the signal from the 6th channel
panelsON = find(panelON>1.5); %these are the frames when the closed-loop resumes
panelsOFF = find(panelON<-1.5); %these are the frames when the 3 sec open-loop bouts start
jsec = panelsOFF/1000; %when the open-loop bouts start in sec
j2sec = panelsON/1000; %when the closed-loop bouts resume in sec


%Plot to check for issues
figure,
plot(panelsON,'ko')
hold on
plot(panelsOFF,'ro')
xlabel('Trial number');ylabel('Frame number');
title('Panel signal pattern');
xlim([0 length(panelsOFF)+1]);


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
data.ficTracIntx = -rawData ( : , xFly); %the negative sign is necessary under my current conditions (z axis facing up)
data.ficTracInty = rawData ( : , yFly); %I think if I wanted to look into this one, I probably need to invert it too



%% See if the panel signal correctly points to when the loop changed
figure('Position', [200 200 1400 800]), 
subplot(1,2,1)
plot(data.xPanelPos,'k.')
hold on
for i = 1:length(panelsOFF)
plot([panelsOFF(i) panelsOFF(i)],[0 96],'r')
end
title('Change to open-loop');

subplot(1,2,2)
plot(data.xPanelPos,'k.')
hold on
for i = 1:length(panelsOFF)-1
plot([panelsON(i+1) panelsON(i+1)],[0 96],'g')
end
title('Back to closed-loop');

%It seems like the jumps are occuring about 10 frames after the A0 signal
%both for the open-loop bouts than for the closed-loop bouts.
%I'm correcting them:

panelsON = panelsON + 10;
panelsOFF = panelsOFF +10;
%% Downsample, unwrap and smooth position data, then get velocity and smooth

%[smoothed] = singleTrialVelocityAnalysis9mm(data,1000);
[smoothed] = posDataDecoding(data,1000);

% for most experiments I have used 1000 Hz as a sample rate, and it is what
% I will use from now on, so that's how I'll leave it, but this could be
% changed in case of need


%% Forward velocity analysis

% The forward velocity is a good indicative of whether the fly is walking
% well in that trial or not

forwardVelocity = smoothed.xVel;
meanVelocity = mean(forwardVelocity);
time = linspace(0,(length(rawData)/1000),length(forwardVelocity));

figure('Position', [300, 200, 1000, 800]),
subplot(2,1,1)
plot(time,forwardVelocity,'k','HandleVisibility','off')%plot(time,forwardVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
hold on 
for i = 1:TrialNum
   plot([jsec(i) jsec(i)],[min(forwardVelocity)-1 max(forwardVelocity)+1],'g','HandleVisibility','off')    
end
ylim([min(forwardVelocity)-1 max(forwardVelocity)+1]);
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



%% Angular velocity analysis

% The forward velocity is a good indicative of whether the fly is walking
% well in that trial or not

angularVelocity = smoothed.angularVel;
meanVelocity = mean(angularVelocity);
time = linspace(0,(length(rawData)/1000),length(forwardVelocity));

figure('Position', [300, 200, 1000, 800]),
subplot(2,1,1)
plot(time,angularVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]); ylim([min(angularVelocity)-5 max(angularVelocity)+5]);
hold on 
for i = 1:TrialNum
   plot([jsec(i) jsec(i)],[min(angularVelocity)-5 max(angularVelocity)+5],'g','HandleVisibility','off')    
end
hline = refline([0 meanVelocity]);
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


figure
set(gcf, 'Position', [500, 500, 1000, 100])
newMap = flipud(gray);
xaxis = time;
trials = [0:1];
imagesc(xaxis,trials,activity')
colormap(newMap)
title(strcat('Activity raster plot, percentage activity:', num2str(percentageActivity), '%'));
ylabel('Activity');
xlabel('Time (s)');


%% Output in degrees of the Panels position

pxToDeg = 360/96; % There are 96 possible positions and this represents 360 deg
posToDeg = zeros(1,length(downsampled.xPanelPos));

% Convert from xpos to degrees, knowing that xpos 70 = 0 deg
for i=1:length(downsampled.xPanelPos)
    if downsampled.xPanelPos(i) == 70
        posToDeg(i) = 0;
    elseif downsampled.xPanelPos(i) >70 
        posToDeg(i) = (downsampled.xPanelPos(i)-70)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        posToDeg(i) = (downsampled.xPanelPos(i)+27)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end

%% Probability density function of the stimulus position

% Remapping the positions to span -180 to 180 deg
remapPosToDeg = wrapTo180(posToDeg);

figure('Position', [100 100 1200 900]),

% Stimulus position in time, with every frame
time = linspace(0,(length(rawData)/1000),length(remapPosToDeg));
subplot(2,4,[1,5])
scatter(remapPosToDeg, time, [], forwardVelocity); %add the forward velocity as a color
hold on
plot(remapPosToDeg, time,'k','HandleVisibility','off')
for i = 1:length(panelsOFF)
    plot([-180 180], [jsec(i) jsec(i)],'r')
end
xlabel('Heading angle (deg)'); ylabel('Time (s)');
title('Angular position of the stimulus');
xlim([-180 180]); ylim([0 max(time)]);
ax = gca;
ax.YDir = 'reverse';

 % Heading in time, with only moving frames.
%time = linspace(0,(length(rawData)/1000),length(flyPos180(forwardVelocity>1)));
subplot(2,4,[4,8])
scatter(remapPosToDeg(forwardVelocity>1), time(forwardVelocity>1));
hold on
for i = 1:length(panelsOFF)
    plot([-180 180], [jsec(i) jsec(i)],'r')
end
xlabel('Heading angle (deg)'); ylabel('Time (s)');
title('Angular position of the stimulus');
xlim([-180 180]); ylim([0 max(time)]);
ax = gca;
ax.YDir = 'reverse'; 


% Plot the histogram and probability density
%With every frame
edges = [-180:20:180];
[counts] = histcounts(remapPosToDeg,edges);
probabilities = counts./sum(counts);
degs = linspace(-180,180,length(counts));

subplot(2,4,2)
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title({'Pb density of the stimulus position','every frame'});
ylabel('Probability density'); xlabel('Stimulus position (deg)');

%With only frames with velocity>1
[countsMoving] = histcounts(remapPosToDeg(forwardVelocity>1),edges);
probabilitiesMoving = countsMoving./sum(countsMoving);
degsMoving = linspace(-180,180,length(countsMoving));

subplot(2,4,3)
plot(degsMoving,probabilitiesMoving,'k')
set(0, 'DefaulttextInterpreter', 'none')
xlim([-180 180]); ylim([0 max(probabilitiesMoving)+0.05]);
title({'Pb density of the stimulus position',' with only moving frames'});
ylabel('Probability density'); xlabel('Stimulus position (deg)');


% Polar coordinates analysis of the stimulus position
posToRad = deg2rad(posToDeg);

%Plot the histogram in polar coordinates
circedges = [0:20:360];
circedges = deg2rad(circedges);
subplot(2,4,6)
polarhistogram(posToRad,circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
title({'Pb density of the stimulus position' , 'every frame'});
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

subplot(2,4,7)
polarhistogram(posToRad(forwardVelocity>1),circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
title({'Pb density of the stimulus position', 'moving frames'});
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';


%% Fly's heading thoughout the experiment

%I think for the fly's heading I don't need to remap anything, cause it
%should be based on fictrac's calibration and have the front be 0?

figure('Position', [100 100 1200 900]),

% Heading in time, with every frame
flyPosPx = rad2deg(smoothed.angularPosition)/pxToDeg;
FlyPosDeg = zeros(1,length(smoothed.angularPosition));

% Convert from xpos to degrees, knowing that xpos 70 = 0 deg
for i=1:length(smoothed.angularPosition)
    if flyPosPx(i) == 70
        FlyPosDeg(i) = 0;
    elseif flyPosPx(i) >70 
        FlyPosDeg(i) = (flyPosPx(i)-70)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        FlyPosDeg(i) = (flyPosPx(i)+27)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end

FlyPos360 = wrapTo360(FlyPosDeg);
flyPos180 = wrapTo180(FlyPos360);

time = linspace(0,(length(rawData)/1000),length(flyPos180));
subplot(2,4,[1,5])
scatter(flyPos180, time, [], forwardVelocity); %add the forward velocity as a color
hold on
plot(flyPos180, time,'k','HandleVisibility','off')
for i = 1:length(panelsOFF)
    plot([-180 180], [jsec(i) jsec(i)],'r')
end
xlabel('Heading angle (deg)'); ylabel('Time (s)');
title('Angular position of the Fly');
xlim([-180 180]); ylim([0 max(time)]);
ax = gca;
ax.YDir = 'reverse';

 % Heading in time, with only moving frames.
%time = linspace(0,(length(rawData)/1000),length(flyPos180(forwardVelocity>1)));
subplot(2,4,[4,8])
scatter(flyPos180(forwardVelocity>1), time(forwardVelocity>1));
hold on
for i = 1:length(panelsOFF)
    plot([-180 180], [jsec(i) jsec(i)],'r')
end
xlabel('Heading angle (deg)'); ylabel('Time (s)');
title('Angular position of the Fly');
xlim([-180 180]); ylim([0 max(time)]);
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

subplot(2,4,2) %with every frame
plot(degsFly,probabilitiesFly,'k')
xlim([-180 180]); ylim([0 max(probabilitiesFly)+0.05]);
title('Fly heading with every frame');
ylabel('Probability density'); xlabel('Fly heading (deg)');

subplot(2,4,3) %with moving frames only
plot(degsFlyMoving,probabilitiesFlyMoving,'k')
xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
title('Fly heading with moving frames');
ylabel('Probability density'); xlabel('Fly heading (deg)');

% In polar coordinates...
%Taking every frame into account
posToRadFly = deg2rad(FlyPos360);
%CircularStatsFly = circ_stats(posToRadFly);
circedges = [0:20:360];
circedges = deg2rad(circedges);
subplot(2,4,6)
polarhistogram(posToRadFly,circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%with only moving frames
%CircularStatsFlyMoving = circ_stats(posToRadFly(forwardVelocity>1));
subplot(2,4,7)
polarhistogram(posToRadFly(forwardVelocity>1),circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

saveas(gcf,strcat(path,'FlyHeading',file(1:end-4),'.png'));

%% Distance to the goal

%Taking all the experiment
figure('Position', [100 100 1200 900]),
%with every frame
goal = circ_mean(posToRadFly,[],2);
dist2goal2 = circ_dist(posToRadFly,goal);
dist2goal = wrapTo180(rad2deg(dist2goal2));
[countsDist] = histcounts(dist2goal,edges);
probabilitiesDist = countsDist./sum(countsDist);
degsFlyDist = linspace(-180,180,length(countsDist));
subplot(1,2,1), plot(degsFlyDist,probabilitiesDist,'r')
title('Distance to the goal with every frame');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

%taking only 'moving frames'
goalMoving = circ_mean(posToRadFly(forwardVelocity>1),[],2);
dist2goalMoving2 = circ_dist(posToRadFly(forwardVelocity>1),goalMoving);
dist2goalMoving = wrapTo180(rad2deg(dist2goalMoving2));
[countsDistMoving] = histcounts(dist2goalMoving,edges);
probabilitiesDistMoving = countsDistMoving./sum(countsDistMoving);
degsFlyDistMoving = linspace(-180,180,length(countsDistMoving));
subplot(1,2,2), plot(degsFlyDistMoving,probabilitiesDistMoving,'r')
title('Distance to the goal with moving frames');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

saveas(gcf,strcat(path,'Dist2Goal',file(1:end-4),'.png'));
%% Plot 2D trajectory

[posx,posy]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.angularPosition);

figure, plot(posx,posy);
axis equal
axis tight
title('2D trajectory of the fly');

%plot with a different color the parts in the trajectory that correspond to
%the 3 sec in open loop

%% Around the jumps for the jump followed by 3 sec in open-loop

%I'm getting the data for the 3 sec before and after each of the jumps that go from closed-loop to open-loop
jumpData = getDataAroundJumpFiltered(rawData,panelsOFF,3,9); 
jumpData.angPos = getPosAroundJumpFiltered(data.ficTracAngularPosition,panelsOFF,3);

%Correct the xPanelPos data for the full non-downsampled raw-data, then
%change to degrees and wrap to 180
fullPosToDeg = data.xPanelPos;

for i=1:length(fullPosToDeg)
    if data.xPanelPos(i) == 70
        fullPosToDeg(i) = 0;
    elseif data.xPanelPos(i) >70 
        fullPosToDeg(i) = (data.xPanelPos(i)-70)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        fullPosToDeg(i) = (data.xPanelPos(i)+27)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end
fullPosToDeg = wrapTo180(fullPosToDeg);
%Get the jump magnitude as the difference in degrees between the panel pos
%20 frames before the jumps and 20 sec after the jumps
jumpData.jumpMag = wrapTo180(fullPosToDeg(panelsOFF+20)-fullPosToDeg(panelsOFF-20));

%plot the aJ data.
time = linspace(-3,3,length(jumpData.angVel));

positiveJumps = find(jumpData.jumpMag>0);
negativeJumps = find(jumpData.jumpMag<0);

figure('Position', [100 100 1200 900]),
subplot(1,2,1)
plot(time,jumpData.angVel(:,positiveJumps))
hold on
plot([0 0],[-300 300],'k')
plot(time,mean(jumpData.angVel(:,positiveJumps),2),'k','lineWidth',2)
title('Angular velocity around the jumps for positive jumps');
xlabel('Time (s)'); ylabel('Angular velocity (deg/s)');
subplot(1,2,2)
plot(time,jumpData.angVel(:,negativeJumps))
hold on
plot([0 0],[-300 300],'k')
plot(time,mean(jumpData.angVel(:,negativeJumps),2),'k','lineWidth',2)
title('Angular velocity around the jumps for negative jumps');
xlabel('Time (s)'); ylabel('Angular velocity (deg/s)');


%All together
positiveMean = mean(jumpData.angVel(:,positiveJumps),2);
positiveStd = std(jumpData.angVel(:,positiveJumps),[],2);

negativeMean = mean(jumpData.angVel(:,negativeJumps),2);
negativeStd = std(jumpData.angVel(:,negativeJumps),[],2);

figure('Position', [100 100 1200 900]),
p1 = boundedline(time,positiveMean,positiveStd/sqrt(length(positiveJumps)),'k','alpha');
hold on
p2= boundedline(time,negativeMean,negativeStd/sqrt(length(negativeJumps)),'r','alpha');
plot([0 0],[-30 30],'color',[.5 .5 .5])
plot([-3 3],[0 0],'color',[.5 .5 .5])
title('Angular velocity around the jumps');
xlabel('Time around the jumps (s)');
ylabel('Angular velocity (deg/s)');
legend([p1,p2],'Positive jumps', 'Negative jumps');

saveas(gcf,strcat(path,'AJangularvelocity',file(1:end-4),'.png'));


%Taking the angular speed (absolute value)
figure('Position', [100 100 1200 900]),
subplot(1,2,1)
plot(time,abs(jumpData.angVel))
hold on
plot([0 0],[-300 300],'k')
ylim([0 300])
title('Angular speed around the jumps')
xlabel('Time (s)'); ylabel('Angular speed (deg/s)');

subplot(1,2,2)
absMean = mean(abs(jumpData.angVel),2);
stdAbsMean = std(abs(jumpData.angVel),[],2);
p3= boundedline(time,absMean,stdAbsMean/sqrt(length(jsec)),'r','alpha');
hold on
plot([0 0],[0 50],'color',[.5 .5 .5]);

saveas(gcf,strcat(path,'AbsAJangularvelocity',file(1:end-4),'.png'));



%If I substract the values before the jump from the values after the jump,
%and then plot that...
respMag = abs(jumpData.angVel - mean(jumpData.angVel(1:75,:)));
figure,
plot(time,mean(respMag,2))
title('Mean response magnitude');

saveas(gcf,strcat(path,'AbsAJangularvelocity2',file(1:end-4),'.png'));

%% Forward velocity around the jumps

figure,
plot(time,jumpData.forwardVel)
hold on
plot([0 0],[min(min(jumpData.forwardVel))-1 max(max(jumpData.forwardVel))+1],'k')
plot(time,mean(jumpData.forwardVel,2),'k','lineWidth',2)
ylim([min(min(jumpData.forwardVel))-1 max(max(jumpData.forwardVel))+1]);
title('Forward velocity around the jumps')
xlabel('Time (s)'); ylabel('Forward velocity (mm/s)');


%% Response magnitude

%Get the magnitude of the response
%I'll take the response between frames 100 and 140.
%I'll take the baseline between frames 30 and 70.
response = mean(jumpData.angVel(100:140,:))-mean(jumpData.angVel(30:70,:));
responseMag = abs(response);

%plot a scatter of the speed vs the jump magnitude (using the difference
%between the startPosition and the pre-jump position).
absJumpMag = abs(jumpData.jumpMag);
figure,
plot(absJumpMag,responseMag,'ro');
xlabel('Jump magnitude (deg)'); ylabel('Response magnitude (deg/s)');
title('Response magnitude as a function of jump magnitude');


saveas(gcf,strcat(path,'ResponseMagvsJumpMag',file(1:end-4),'.png'));

%% Evolution of the responses in time

figure('Position', [100 100 1200 900]),
% subplot(1,2,1)
% colormap(hot)
% imagesc(time,[1:length(panelsOFF)],abs(jumpData.angVel'))
% colorbar
% hold on
% plot([0, 0], [1,length(panelsOFF)],'k','LineWidth',2);
% xlim([-5,5]);
% title('Angular speed in time');
% xlabel('Time around jump (s)'); ylabel('Trial number');
% 
% subplot(1,2,2)
plot(responseMag,'ko')
title('Response magnitude in time');
xlabel('Trial number'); ylabel('Response magnitude (deg/s)');



%% Around the jumps for the jump going back to closed loop

panelsON2 = panelsON(2:end);
jumpData2 = getDataAroundJump(rawData,panelsON2,3,9);
jumpData2.angPos = getPosAroundJump(data.ficTracAngularPosition,panelsON2,3);

%Get the jump magnitude as the difference in degrees between the panel pos
%20 frames before the jumps and 20 sec after the jumps
jumpData2.jumpMag = wrapTo180(fullPosToDeg(panelsON2+20)-fullPosToDeg(panelsON2-20));

%plot the aJ data.
positiveJumps2 = find(jumpData2.jumpMag>0);
negativeJumps2 = find(jumpData2.jumpMag<0);

figure('Position', [100 100 1200 900]),
subplot(1,2,1)
plot(time,jumpData2.angVel(:,positiveJumps2))
hold on
plot([0 0],[-300 300],'k')
plot(time,mean(jumpData2.angVel(:,positiveJumps2),2),'k','lineWidth',2)
title('Angular velocity around the jumps for positive jumps');
xlabel('Time (s)'); ylabel('Angular velocity (deg/s)');
subplot(1,2,2)
plot(time,jumpData2.angVel(:,negativeJumps2))
hold on
plot([0 0],[-300 300],'k')
plot(time,mean(jumpData2.angVel(:,negativeJumps2),2),'k','lineWidth',2)
title('Angular velocity around the jumps for negative jumps');
xlabel('Time (s)'); ylabel('Angular velocity (deg/s)');


%All together
positiveMean2 = mean(jumpData2.angVel(:,positiveJumps2),2);
positiveStd2 = std(jumpData2.angVel(:,positiveJumps2),[],2);

negativeMean2 = mean(jumpData2.angVel(:,negativeJumps2),2);
negativeStd2 = std(jumpData2.angVel(:,negativeJumps2),[],2);

figure('Position', [100 100 1200 900]),
p12 = boundedline(time,positiveMean2,positiveStd2/sqrt(length(positiveJumps2)),'k','alpha');
hold on
p22= boundedline(time,negativeMean2,negativeStd2/sqrt(length(negativeJumps2)),'r','alpha');
plot([0 0],[-30 30],'color',[.5 .5 .5])
plot([-3 3],[0 0],'color',[.5 .5 .5])
title('Angular velocity around the jumps');
xlabel('Time around the jumps (s)');
ylabel('Angular velocity (deg/s)');
legend([p12,p22],'Positive jumps', 'Negative jumps');

saveas(gcf,strcat(path,'AJBackClosedAngularvelocity',file(1:end-4),'.png'));


%Taking the angular speed (absolute value)
figure('Position', [100 100 1200 900]),
subplot(1,2,1)
plot(time,abs(jumpData2.angVel))
hold on
plot([0 0],[-300 300],'k')
ylim([0 300])
title('Angular speed around the jumps')
xlabel('Time (s)'); ylabel('Angular speed (deg/s)');

subplot(1,2,2)
absMean2 = mean(abs(jumpData2.angVel),2);
stdAbsMean2 = std(abs(jumpData2.angVel),[],2);
p32= boundedline(time,absMean2,stdAbsMean2/sqrt(length(jsec)),'r','alpha');
hold on
plot([0 0],[0 50],'color',[.5 .5 .5]);

saveas(gcf,strcat(path,'AbsAJBacktoClosedangularvelocity',file(1:end-4),'.png'));



%If I substract the values before the jump from the values after the jump,
%and then plot that...
respMag2 = abs(jumpData2.angVel - mean(jumpData2.angVel(1:75,:)));
figure,
plot(time,mean(respMag2,2))
title('Mean response magnitude');

saveas(gcf,strcat(path,'AbsAJBacktoClosedangularvelocity2',file(1:end-4),'.png'));


%% Forward velocity around the jumps

figure,
plot(time,jumpData2.forwardVel)
hold on
plot([0 0],[min(min(jumpData2.forwardVel))-1 max(max(jumpData2.forwardVel))+1],'k')
plot(time,mean(jumpData2.forwardVel,2),'k','lineWidth',2)
ylim([min(min(jumpData2.forwardVel))-1 max(max(jumpData2.forwardVel))+1]);
title('Forward velocity around the jumps')
xlabel('Time (s)'); ylabel('Forward velocity (mm/s)');


%% Response magnitude

%Get the magnitude of the response
%I'll take the response between frames 100 and 140.
%I'll take the baseline between frames 30 and 70.
response2 = mean(jumpData2.angVel(100:140,:))-mean(jumpData2.angVel(30:70,:));
responseMag2 = abs(response2);

%plot a scatter of the speed vs the jump magnitude (using the difference
%between the startPosition and the pre-jump position).
absJumpMag2 = abs(jumpData2.jumpMag);
figure,
plot(absJumpMag2,responseMag2,'ro');
xlabel('Jump magnitude (deg)'); ylabel('Response magnitude (deg/s)');
title('Response magnitude as a function of jump magnitude');


saveas(gcf,strcat(path,'ResponseMagvsJumpMagBackToClosed',file(1:end-4),'.png'));

%% Evolution of the responses in time

figure('Position', [100 100 1200 900]),
% subplot(1,2,1)
% colormap(hot)
% imagesc(time,[1:length(panelsOFF)],abs(jumpData.angVel'))
% colorbar
% hold on
% plot([0, 0], [1,length(panelsOFF)],'k','LineWidth',2);
% xlim([-5,5]);
% title('Angular speed in time');
% xlabel('Time around jump (s)'); ylabel('Trial number');
% 
% subplot(1,2,2)
plot(responseMag2,'ko')
title('Response magnitude in time');
xlabel('Trial number'); ylabel('Response magnitude (deg/s)');
