% Experiment 14 analysis

%%% This code analyses the output of the panels and the FicTrac data
%for experiment 14
close all; clear all;

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment14'
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

    
%% Subset acquisition of x and y pos, as well as FicTrac data

data.xpanelVolts =  rawData (:,xPanels); 
VOLTAGE_RANGE_x = 10;
maxValX =  96 ;

data.yPanelVolts =  rawData (:, yPanels);
VOLTAGE_RANGE_y = 10;
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
subplot(2,1,1)
plot(rawData(:,6))
ylabel('Voltage signal (V)');
title('Panels ON and OFF');
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
%Find the frame when the panels turn on
[M, I] = max(panelON);
%Find the frame when they turn off.
[M2, I2] = min(panelON);
startFrame = I(1);
endFrame = I2(1);
%add text indicating start and end frames
text(startFrame+35000,5,strcat('startFrame',num2str(startFrame)))
text(endFrame-100000,5,strcat('endFrame',num2str(endFrame)))


% get the jumps using the signal from the yPanels channel
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
%j = j(2:end);
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

%% 

%Now that we have the bar jumps, let's bring the y voltage down to be
%contained similarly for all 3 blocks

%first we need to find the jumps that signal a change in block
remove1 = abs(j-900000); %the first change occurs around 980 s
[M,I] = min(remove1); %find the jump index that's closer to that time
changeBlock1 = j(I);

remove2 = abs(j-1800000); %the second change occurs around 1960 s
[M2,I2] = min(remove2); %find the jump index that's closer to that time
changeBlock2 = j(I2);

voltageCorr = data.yPanelVolts(j(I)+1) - data.yPanelVolts(j(I)-1);

figure,
subplot(2,1,1)
plot(data.yPanelVolts)
title('Y panel output before deleting the block change jumps');

data.yPanelVolts(j(I)+1:j(I2)) = data.yPanelVolts(j(I)+1:j(I2)) - voltageCorr;
data.yPanelVolts = 2*data.yPanelVolts;

%check that it now looks fine
subplot(2,1,2)
plot(data.yPanelVolts)
title('Y panel output after deleting the block change jumps');

%remove values from jumps vector
indices = [I,I2];
j(indices) = [];
jsec = j/1000;

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
data.xPanelPos = round ((data.xPanelVolts  * maxValX) /VOLTAGE_RANGE_x); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels
data.yPanelPos = round ((data.yPanelVolts  * maxValY) /VOLTAGE_RANGE_y);

%check the wrapping graphically
figure('Position', [100 100 1200 900]),
subplot(2,1,1),
plot(data.yPanelVolts,'r')
hold on
plot(data.xPanelVoltsUW,'b')
plot(data.xPanelVolts,'g')
ylabel('Voltage (V)');
title('x panels voltage signals before and after wrapping');
ylim([0 max(data.xPanelVoltsUW)]);
legend({'yPanels','Unwrapped xPanels', 'Wrapped xPanels'});

subplot(2,1,2),
plot(data.yPanelVolts,'r')
hold on
plot(data.ficTracAngularPositionUW,'b')
plot(data.ficTracAngularPosition,'g')
xlabel('Time (frame)');ylabel('Voltage (V)');
title('angular position voltage signals before and after wrapping');
ylim([0 max(data.ficTracAngularPositionUW)]);
legend({'yPanels','Unwrapped angular position', 'Wrapped angular position'});


% Getting the degrees that each jump represents and checking they are
% correct

figure('Position', [100 100 1200 900]),
jumpPos = data.yPanelPos(j+1)-data.yPanelPos(j-1);
degJumps = wrapTo180(jumpPos*(360/96));
subplot(1,3,1), plot(degJumps,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from yPanels');
ylabel('deg');xlabel('Trial #');xlim([1 length(j)]);
% compare that to the jump function we have stored
hold on
plot(jumps,'b')
legend({'Jumps from Y data','Jump function used'});

%Check if the jump magnitude appears ok in the x panel data
jumpMag = data.xPanelPos(j+1)-data.xPanelPos(j-1); 
degMag = wrapTo180(jumpMag*(360/96));
subplot(1,3,2), plot(degMag,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from xPanels');
ylabel('deg');xlabel('Trial #'); xlim([1 length(j)]);
% compare that to the jump function we have stored
hold on
plot(jumps,'b')
legend({'Jumps from X data','Jump function used'});

%check if the jump magnitude appears ok in the angular position data
jumpMag2 = data.ficTracAngularPosition(j+1)-data.ficTracAngularPosition(j-1); 
radMag2 = jumpMag2.* 2 .* pi ./ 10; %convert from voltage to radians
degMag2 = wrapTo180(rad2deg(radMag2)); %convert from radians to degrees and wrap 180
subplot(1,3,3), plot(degMag2,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from angular position');
ylabel('deg');xlabel('Trial #');xlim([1 length(j)]);
% compare that to the jump function we have stored
hold on
plot(jumps,'b')
legend({'Jumps from angular position','Jump function used'});


% If the jumps look different because the function got cut off because of
% the frame number, we create a newJumps vector with the real jump values:

 newJumps = zeros(size(j));
 
 for i = 1:length(newJumps)
    if degJumps(i)>-100 & degJumps(i)<-80
        newJumps(i) = -90;
    else
        newJumps(i) = 90;
    end
 end
 
 
figure('Position', [100 100 1200 900]),
jumpPos = data.yPanelPos(j+1)-data.yPanelPos(j-1);
degJumps = wrapTo180(jumpPos*(360/96));
subplot(1,3,1), plot(degJumps,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from yPanels');
ylabel('deg');xlabel('Trial #');xlim([1 length(j)]);
% compare that to the jump function we have stored
hold on
plot(newJumps,'b')
legend({'Jumps from Y data','Jump function used'});

%Check if the jump magnitude appears ok in the x panel data
jumpMag = data.xPanelPos(j+1)-data.xPanelPos(j-1); 
degMag = wrapTo180(jumpMag*(360/96));
subplot(1,3,2), plot(degMag,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from xPanels');
ylabel('deg');xlabel('Trial #'); xlim([1 length(j)]);
% compare that to the jump function we have stored
hold on
plot(newJumps,'b')
legend({'Jumps from X data','Jump function used'});

%check if the jump magnitude appears ok in the angular position data
jumpMag2 = data.ficTracAngularPosition(j+1)-data.ficTracAngularPosition(j-1); 
radMag2 = jumpMag2.* 2 .* pi ./ 10; %convert from voltage to radians
degMag2 = wrapTo180(rad2deg(radMag2)); %convert from radians to degrees and wrap 180
subplot(1,3,3), plot(degMag2,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from angular position');
ylabel('deg');xlabel('Trial #');xlim([1 length(j)]);
% compare that to the jump function we have stored
hold on
plot(newJumps,'b')
legend({'Jumps from angular position','Jump function used'});



%% Downsample, unwrap and smooth position data, then get velocity and smooth

sizeBall = 9;
%[smoothed] = singleTrialVelocityAnalysis9mm(data,1000);
[smoothed] = posDataDecoding(data,1000);

% for most experiments I have used 1000 Hz as a sample rate, and it is what
% I will use from now on, so that's how I'll leave it, but this could be
% changed in case of need


%% 2D trajectories

%Uncorrected trajectory (without correcting the heading). This is to
%visualize the fly's response to the jumps, when they do very well.
[posx,posy]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.AngularPosition);

time = linspace(1,1000,length(posx));

figure, 
subplot(1,2,1)
scatter(posx,posy,0.5,time)
colorbar
axis equal
axis tight
title('Uncorrected 2D trajectory');


%Corrected trajectory: this is the actual trajectory the fly took during
%the block
[posx2,posy2]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.angularPosition);

subplot(1,2,2), scatter(posx2,posy2,0.5,time);
colorbar
axis equal
axis tight
title('Corrected 2D trajectory');

saveas(gcf,strcat(path,'Trajectories',file(1:end-4),'.png'));


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
figure ('Position', [100 100 1200 900]),
subplot(2,1,1)
plot(time,forwardVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
ylim([min(forwardVelocity)-10 max(forwardVelocity)+10]);
hold on
%add the jumps
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
xlim([min(forwardVelocity) max(forwardVelocity)]);

saveas(gcf,strcat(path,'ForwardVelocity',file(1:end-4),'.png'));



%% Angular velocity

%To look at the velocity, I should do so without adding the offsets of the
%jumps, because the offsets will give me weird jumps in velocity that the
%fly might not actually have made. For the angular velocity, I'm using a new 
%function that I made to downsample, unwrap, smooth and get the velocity
%with the angular position uncorrected for the offset

AngularPosition = rawData ( : , headingFly);
%angularVelocity = getAngVel(AngularPosition);
angularVelocity = smoothed.AngularVel;

meanAngVelocity = mean(angularVelocity);
time = linspace(0,(length(rawData)/1000),length(angularVelocity));

figure('Position', [100 100 1200 900]),
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
xlim([min(angularVelocity) max(angularVelocity)]);

saveas(gcf,strcat(path,'SmoothedAngularVelocity',file(1:end-4),'.png'))

%% Plot forward and angular velocities throughout the experiment using a raster plot

newMap = flipud(gray);
figure
set(gcf, 'Position', [300, 500, 1600, 500]),
subplot(2,1,1)
imagesc(time,[],forwardVelocity')
colormap(hot)
colorbar
xlabel('Time (s)');
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('Forward Velocity (mm/s)')

subplot(2,1,2)
imagesc(time,[],angularVelocity')
colormap(hot)
colorbar
xlabel('Time (s)');
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('Angular Velocity (deg/s)')

saveas(gcf,strcat(path,'VelocityRP',file(1:end-4),'.png'))



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
newMap = flipud(gray);
xaxis = time;
trials = [0:1];
imagesc(xaxis,trials,activity')
colormap(newMap)
title(strcat('Activity raster plot, percentage activity:', num2str(percentageActivity), '%'));
ylabel('Activity');
xlabel('Time (s)');

save(strcat(path,'Activity',file(1:end-4),'.mat'),'time','activity','percentageActivity');
saveas(gcf,strcat(path,'ActivityRP',file(1:end-4),'.png'))



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

figure('Position', [100 100 1600 900]),

% Stimulus position in time, with every frame
time = linspace(0,(length(rawData)/1000),length(remapPosToDeg));
subplot(3,4,[1,5,9])
scatter(remapPosToDeg, time, [], forwardVelocity); %add the forward velocity as a color
hold on
plot(remapPosToDeg, time,'k','HandleVisibility','off')
xlabel('Heading angle (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [changeBlock1/1000 changeBlock1/1000],'r','LineWidth',2);
     plot([-180 180], [changeBlock2/1000 changeBlock2/1000],'r','LineWidth',2);
end
ax = gca;
ax.YDir = 'reverse';

 % Heading in time, with only moving frames.
subplot(3,4,[4,8,12])
scatter(remapPosToDeg(forwardVelocity>1), time(forwardVelocity>1));
xlabel('Heading angle (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [changeBlock1/1000 changeBlock1/1000],'r','LineWidth',2);
     plot([-180 180], [changeBlock2/1000 changeBlock2/1000],'r','LineWidth',2);
end
ax = gca;
ax.YDir = 'reverse'; 


% Plot the polar histograms

% Polar coordinates analysis of the stimulus position
posToRad = deg2rad(posToDeg);

%Plot the histogram in polar coordinates
circedges = [0:20:360];
circedges = deg2rad(circedges);
%Block1
subplot(3,4,2)
polarhistogram(posToRad(1:round(length(posToRad)/3)),circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%Block2
subplot(3,4,6)
polarhistogram(posToRad(round(length(posToRad)/3):(round(length(posToRad)/3))*2),circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%Block3
subplot(3,4,10)
polarhistogram(posToRad((round(length(posToRad)/3))*2:(round(length(posToRad)/3))*3),circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top


%with moving frames
movingPos = posToRad(forwardVelocity>1);

%Block1
subplot(3,4,3)
polarhistogram(movingPos(1:round(length(movingPos)/3)),circedges,'Normalization','probability','FaceColor',[0,0.5,0.3],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%Block2
subplot(3,4,7)
polarhistogram(movingPos(round(length(movingPos)/3):(round(length(movingPos)/3))*2),circedges,'Normalization','probability','FaceColor',[0,0.5,0.3],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%Block3
subplot(3,4,11)
polarhistogram(movingPos((round(length(movingPos)/3))*2:end),circedges,'Normalization','probability','FaceColor',[0,0.5,0.3],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

saveas(gcf,strcat(path,'BarPosition',file(1:end-4),'.png'))
saveas(gcf,strcat(path,'BarPosition',file(1:end-4),'.svg'))