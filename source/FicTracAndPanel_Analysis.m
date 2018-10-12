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


%% Forward velocity analysis

% The forward velocity is a good indicative of whether the fly is walking
% well in that trial or not

forwardVelocity = smoothed.xVel;
meanVelocity = mean(forwardVelocity);
time = linspace(0,(length(rawData)/1000),length(forwardVelocity));

figure,
subplot(2,1,1)
plot(time,forwardVelocity,'k','HandleVisibility','off')
xlim = ([0 time(end)]);
hold on
hline = refline([0 meanVelocity]);
hline.Color = 'r'; hline.LineStyle = '--';
rline = refline([0 0]);
rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
title({'Forward velocity of the fly', typeOfStim}, 'Interpreter', 'none');
xlabel('Time (s)')
ylabel('Velocity (mm/s)')
legend('Mean forward velocity');

subplot(2,1,2)
histogram(forwardVelocity,'FaceColor',[.4 .2 .6])
title('Distribution of forward velocities');
xlabel('Forward velocity (mm/s)');
ylabel('Frequency');

saveas(gcf,strcat(path,'ForwardVelocity_ExpNum', file(end-4), '.png'))



%%  Keep the frames during which the fly is moving

% We are going to decide whether a fly is moving or not based on the
% forward velocity. If it's above 0.7 mm/s we will consider it is moving
% We will work with the downsampled data

downsampled.xPanelPos = downsample(data.xPanelPos,1000/25);
dataMoving.xPanelPos = downsampled.xPanelPos(forwardVelocity>0.7);
moving = smoothed.angularPosition(forwardVelocity>0.7);

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

% Create a time vector in sec
time = linspace(0,(size(rawData,1)/1000),size(dataMoving.xPanelPos,1)); %1000 is our sampling rate

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
%xlim([-180, 180]); ylim([0 max(probabilities)+0.05]);
title('Histogram of the stimulus position');
ylabel('Probability'); xlabel('Stimulus position (deg)');

subplot(1,2,2),
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')
suptitle(typeOfStim);
%xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
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

%% Polar coordinates analysis of the stimulus position

%Plot the histogram in polar coordinates
posToRad = posToDeg.*pi/180;
figure,
polarhistogram(posToRad,20,'Normalization','pdf','FaceColor',[0.6,0.3,0.3]);
title({'Probability density of the stimulus position';typeOfStim});
% hold on
% polarhistogram(deg2rad(noPanelDeg),20,'FaceColor','k');

saveas(gcf,strcat(path,'PolarHistogramStimPosition_ExpNum', file(end-4), '.png'))

% some statistics...

CircularStats = circ_stats(posToRad);
[pval,z] = circ_rtest(posToRad);

%% Angular position of the stimulus in time

figure,
plot(time,wrapTo360(posToDeg),'k','HandleVisibility','off')
ylabel('Heading angle (deg)'); xlabel('Time (s)');
ylim([0 360]);
title('Angular position of the stimulus as a function of time');;
hline = refline([0 wrapTo360(rad2deg(CircularStats.median))]);
set(hline,'Color',[1,0,0])
legend('Median heading');

saveas(gcf,strcat(path,'AngulaPosStimInTime_ExpNum', file(end-4), '.png'))
%% Probability density of the fly heading

flyPosToDegMoving = moving.*36; %if 10 V = 360 deg, then xV = 36x deg

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
%xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
title('Histogram of the fly heading');
ylabel('Probability'); xlabel('Fly heading (deg)');

subplot(1,2,2),
plot(degsFlyMoving,probabilitiesFlyMoving,'k')
suptitle(typeOfStim)
%xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
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

FlyPosToRad = flyPosToDegMoving.*pi/180;
figure,
polarhistogram(FlyPosToRad,20,'Normalization','pdf','FaceColor',[0.6,0.3,0.3]);
title({'Probability density of the fly heading';typeOfStim});
saveas(gcf,strcat(path,'PolarHistFlyHeading_ExpNum', file(end-4), '.png'))

CircularStatsFly = circ_stats(FlyPosToRad);
%% Angular position of the fly in time

figure,
plot(time,wrapTo360(flyPosToDegMoving),'k','HandleVisibility','off')
ylabel('Heading angle (deg)'); xlabel('Time (s)');
title('Angular position of the Fly as a function of time');
hline = refline([0 wrapTo360(rad2deg(CircularStatsFly.median))]);
set(hline,'Color',[1,0,0])
legend('Median heading');

saveas(gcf,strcat(path,'AngulaPosFlyInTime_ExpNum', file(end-4), '.png'))


%% Plot 2D virtual trajectory of the fly

dataMoving.ficTracIntx = data.ficTracIntx(forwardVelocity>0.7);
dataMoving.ficTracInty = data.ficTracInty(forwardVelocity>0.7);

figure, plot(dataMoving.ficTracIntx,dataMoving.ficTracInty,'o')
% 
% z = zeros(size(dataMoving.ficTracIntx));
% col = dataMoving.ficTracIntx;  % This is the color, vary with x in this case.
% surface([dataMoving.ficTracIntx;dataMoving.ficTracIntx],[dataMoving.ficTracInty;dataMoving.ficTracInty],[z;z],[col;col],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',2);