function [smoothed,posToRadFly,degsFlyDistMoving,probabilitiesDistMoving] = runClosedLoopAnalysis(rawData,path,trialName)

%This function analyses closed-loop data

%OUTPUT




%INPUT
%rawData = the data to be analyzed
%trialName = the name of the trial %I need to add a number to this as well.
%path = the folder where the data is kept



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
data.ficTracIntx = -rawData ( : , xFly); %the negative sign is necessary under my current conditions (z axis facing up)
data.ficTracInty = rawData ( : , yFly); %I think if I wanted to look into this one, I probably need to invert it too


%% Downsample, unwrap and smooth position data, then get velocity and smooth

[smoothed] = posDataDecoding(data,1000);


%% Forward velocity analysis

% The forward velocity is a good indicative of whether the fly is walking
% well in that trial or not

forwardVelocity = smoothed.xVel;
meanVelocity = mean(forwardVelocity);
time = linspace(0,(length(rawData)/1000),length(forwardVelocity));

% figure,
% subplot(2,1,1)
% plot(time,forwardVelocity,'k','HandleVisibility','off')
% xlim([0 time(end)]);
% hold on
% hline = refline([0 meanVelocity]);
% hline.Color = 'r'; hline.LineStyle = '--';
% rline = refline([0 0]);
% rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
% title('Forward velocity of the fly', 'Interpreter', 'none');
% xlabel('Time (s)')
% ylabel('Velocity (mm/s)')
% legend('Mean forward velocity');
% 
% subplot(2,1,2)
% histogram(forwardVelocity,'FaceColor',[.4 .2 .6])
% title('Distribution of forward velocities');
% xlabel('Forward velocity (mm/s)');
% ylabel('Frequency');


%% Angular velocity analysis

% angularVelocity = smoothed.angularVel;
% meanAngVelocity = mean(angularVelocity);
% time = linspace(0,(length(rawData)/1000),length(angularVelocity));
% 
% figure,
% subplot(2,1,1)
% plot(time,angularVelocity,'k','HandleVisibility','off')
% xlim([0 time(end)]);
% hold on
% hline = refline([0 meanAngVelocity]);
% hline.Color = 'r'; hline.LineStyle = '--';
% rline = refline([0 0]);
% rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
% title('Angular velocity of the fly', 'Interpreter', 'none');
% xlabel('Time (s)')
% ylabel('Velocity (deg/s)')
% legend('Mean angular velocity');
% 
% subplot(2,1,2)
% histogram(angularVelocity,'FaceColor',[.2 .8 .6])
% title('Distribution of angular velocities');
% xlabel('Angular velocity (deg/s)');
% ylabel('Frequency');


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


%A more accurate plot
% figure
% set(gcf, 'Position', [500, 500, 1000, 100])
% newMap = flipud(gray);
% xaxis = time;
% trials = [0:1];
% imagesc(xaxis,trials,activity')
% colormap(newMap)
% title(strcat('Activity raster plot, percentage activity:', num2str(percentageActivity), '%'));
% ylabel('Activity');
% xlabel('Time (s)');


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

% figure('Position', [100 100 1200 900]),
% 
% % Stimulus position in time, with every frame
% time = linspace(0,(length(rawData)/1000),length(remapPosToDeg));
% subplot(2,4,[1,5])
% scatter(remapPosToDeg, time, [], forwardVelocity); %add the forward velocity as a color
% hold on
% plot(remapPosToDeg, time,'k','HandleVisibility','off')
% xlabel('Heading angle (deg)'); ylabel('Time (s)');
% title('Angular position of the stimulus');
% xlim([-180 180]); ylim([0 max(time)]);
% ax = gca;
% ax.YDir = 'reverse';
% 
%  % Heading in time, with only moving frames.
% %time = linspace(0,(length(rawData)/1000),length(flyPos180(forwardVelocity>1)));
% subplot(2,4,[4,8])
% scatter(remapPosToDeg(forwardVelocity>1), time(forwardVelocity>1));
% xlabel('Heading angle (deg)'); ylabel('Time (s)');
% title('Angular position of the stimulus');
% xlim([-180 180]); ylim([0 max(time)]);
% ax = gca;
% ax.YDir = 'reverse'; 
% 
% 
% % Plot the histogram and probability density
% %With every frame
% edges = [-180:20:180];
% [counts] = histcounts(remapPosToDeg,edges);
% probabilities = counts./sum(counts);
% degs = linspace(-180,180,length(counts));
% 
% subplot(2,4,2)
% plot(degs,probabilities,'k')
% set(0, 'DefaulttextInterpreter', 'none')
% xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
% title('Probability density of the stimulus position');
% ylabel('Probability density'); xlabel('Stimulus position (deg)');
% 
% %With only frames with velocity>1
% [countsMoving] = histcounts(remapPosToDeg(forwardVelocity>1),edges);
% probabilitiesMoving = countsMoving./sum(countsMoving);
% degsMoving = linspace(-180,180,length(countsMoving));
% 
% subplot(2,4,3)
% plot(degsMoving,probabilitiesMoving,'k')
% set(0, 'DefaulttextInterpreter', 'none')
% xlim([-180 180]); ylim([0 max(probabilitiesMoving)+0.05]);
% title('Probability density of the stimulus position with only moving frames');
% ylabel('Probability density'); xlabel('Stimulus position (deg)');


% Polar coordinates analysis of the stimulus position
posToRad = deg2rad(posToDeg);

%% Plot the histogram in polar coordinates
% circedges = [0:20:360];
% circedges = deg2rad(circedges);
% subplot(2,4,6)
% polarhistogram(posToRad,circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
% title('Probability density of the stimulus position, every frame');
% ax = gca;
% ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
% 
% subplot(2,4,7)
% polarhistogram(posToRad(forwardVelocity>1),circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
% title('Probability density of the stimulus position, moving frames');
% ax = gca;
% ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top';


%% Fly's heading thoughout the experiment

%I think for the fly's heading I don't need to remap anything, cause it
%should be based on fictrac's calibration and have the front be 0?

% figure('Position', [100 100 1200 900]),

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

% time = linspace(0,(length(rawData)/1000),length(flyPos180));
% subplot(2,4,[1,5])
% scatter(flyPos180, time, [], forwardVelocity); %add the forward velocity as a color
% hold on
% plot(flyPos180, time,'k','HandleVisibility','off')
% xlabel('Heading angle (deg)'); ylabel('Time (s)');
% title('Angular position of the Fly');
% xlim([-180 180]); ylim([0 max(time)]);
% ax = gca;
% ax.YDir = 'reverse';
% 
%  % Heading in time, with only moving frames.
% subplot(2,4,[4,8])
% scatter(flyPos180(forwardVelocity>1), time(forwardVelocity>1));
% xlabel('Heading angle (deg)'); ylabel('Time (s)');
% title('Angular position of the Fly');
% xlim([-180 180]); ylim([0 max(time)]);
% ax = gca;
% ax.YDir = 'reverse'; 
%  
% 
% % Plot the histogram and probability density
edges = [-180:20:180];
% [countsFly] = histcounts(flyPos180,edges);
% [countsFlyMoving] = histcounts(flyPos180(forwardVelocity>1),edges);
% probabilitiesFly = countsFly./sum(countsFly);
% probabilitiesFlyMoving = countsFlyMoving./sum(countsFlyMoving);
% degsFly = linspace(-180,180,length(countsFly));
% degsFlyMoving = linspace(-180,180,length(countsFlyMoving));
% 
% subplot(2,4,2) %with every frame
% plot(degsFly,probabilitiesFly,'k')
% xlim([-180 180]); ylim([0 max(probabilitiesFly)+0.05]);
% title('Fly heading with every frame');
% ylabel('Probability density'); xlabel('Fly heading (deg)');
% 
% subplot(2,4,3) %with moving frames only
% plot(degsFlyMoving,probabilitiesFlyMoving,'k')
% xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
% title('Fly heading with moving frames');
% ylabel('Probability density'); xlabel('Fly heading (deg)');

% In polar coordinates...
%Taking every frame into account
posToRadFly = deg2rad(FlyPos360);
% circedges = [0:20:360];
% circedges = deg2rad(circedges);
% subplot(2,4,6)
% polarhistogram(posToRadFly,circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
% ax = gca;
% ax.ThetaDir='clockwise';
% ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
% 
% %with only moving frames
% subplot(2,4,7)
% polarhistogram(posToRadFly(forwardVelocity>1),circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
% ax = gca;
% ax.ThetaDir = 'clockwise';
% ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
% 
% saveas(gcf,strcat(path,'FlyPosition_',trialName,'.png'));


%% Distance to the goal

%Taking all the experiment
% figure('Position', [100 100 1200 900]),
%with every frame
goal = circ_mean(posToRadFly,[],2);
dist2goal2 = circ_dist(posToRadFly,goal);
dist2goal = wrapTo180(rad2deg(dist2goal2));
[countsDist] = histcounts(dist2goal,edges);
probabilitiesDist = countsDist./sum(countsDist);
degsFlyDist = linspace(-180,180,length(countsDist));
% subplot(1,2,1), plot(degsFlyDist,probabilitiesDist,'r')
% title('Distance to the goal with every frame');
% xlim([-180, 180]); xlabel('Distance to the goal (deg)');
% ylabel('Probability');

%taking only 'moving frames'
goalMoving = circ_mean(posToRadFly(forwardVelocity>1),[],2);
dist2goalMoving2 = circ_dist(posToRadFly(forwardVelocity>1),goalMoving);
dist2goalMoving = wrapTo180(rad2deg(dist2goalMoving2));
[countsDistMoving] = histcounts(dist2goalMoving,edges);
probabilitiesDistMoving = countsDistMoving./sum(countsDistMoving);
degsFlyDistMoving = linspace(-180,180,length(countsDistMoving));
% subplot(1,2,2), plot(degsFlyDistMoving,probabilitiesDistMoving,'r')
% title('Distance to the goal with moving frames');
% xlim([-180, 180]); xlabel('Distance to the goal (deg)');
% ylabel('Probability');

%saveas(gcf,strcat(path,'Dist2goal_',trialName,'.png'));

%% Plot 2D trajectory

% [posx,posy]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.angularPosition);
% 
% figure, plot(posx,posy);
% axis equal
% axis tight
% title('2D trajectory of the fly');
% 
% saveas(gcf,strcat(path,'Trajectory_',trialName,'.png'));


end