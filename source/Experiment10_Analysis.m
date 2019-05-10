% Experiment 10 analysis

%%% This code analyses the output of the panels and the FicTrac data
%for experiment 10, in which the fly gets bouts of gratings in closed-loop
%alternated with optomotor response trials
close all; clear all;

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10'
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
panelsON = find(panelON>1.5);
panelsOFF = find(panelON<-1.5);


%Plot to check for issues
figure,
plot(panelsON,'ko')
hold on
plot(panelsOFF,'ro')
xlabel('Trial number');ylabel('Frame number');
title('Panel signal pattern');
xlim([0 length(panelsOFF)+1]);
%% Separating the data per trial type

%All of the odd trials are closed-loop trials (i.e., between the first
%panelON and the first panelOFF, between the 3rd panelON and the third
%panelOFF, etc,...)
oddJumps = [1:2:length(panelsOFF)-1];

for i = 1:length(oddJumps)
    closedloopTrials{i} = daq_data(:,panelsON(oddJumps(i)):panelsOFF(oddJumps(i)));
end

%First take all of the optomotor response trials
evenJumps = [2:2:length(panelsOFF)];

for i = 1:length(evenJumps)
    optomotorResponseTrials{i} = daq_data(:,panelsON(evenJumps(i)):panelsOFF(evenJumps(i)));
end

%Use the jumpFunctions stored to determine the direction of rotation of the
%stimulus
clockwiseTrials = optomotorResponseTrials(jumpFunction(1:length(optomotorResponseTrials)) == 52);
counterclockwiseTrials = optomotorResponseTrials(jumpFunction(1:length(optomotorResponseTrials)) == 53);

    
%% Subset acquisition of x and y pos, as well as FicTrac data

    VOLTAGE_RANGE_x = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
    maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
    VOLTAGE_RANGE_y = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
    maxValY = 96;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2

    
for i = 1:length(closedloopTrials)
    dataSF{i}.xpanelVolts =  closedloopTrials{i}(xPanels,:); 
    dataSF{i}.yPanelVolts =  closedloopTrials{i}(yPanels,:);
    dataSF{i}.ficTracAngularPosition = closedloopTrials{i}(headingFly,:); 
    dataSF{i}.ficTracIntx = -closedloopTrials{i}(xFly,:); 
    dataSF{i}.ficTracInty = closedloopTrials{i}(yFly,:); 
end

for i = 1:length(clockwiseTrials)
    dataClockwise{i}.xpanelVolts =  clockwiseTrials{i}(xPanels,:); 
    dataClockwise{i}.yPanelVolts =  clockwiseTrials{i}(yPanels,:);
    dataClockwise{i}.ficTracAngularPosition = clockwiseTrials{i}(headingFly,:); 
    dataClockwise{i}.ficTracIntx = -clockwiseTrials{i}(xFly,:); 
    dataClockwise{i}.ficTracInty = clockwiseTrials{i}(yFly,:); 
end
    
for i = 1:length(counterclockwiseTrials)
    dataCounterclockwise{i}.xpanelVolts =  counterclockwiseTrials{i}(xPanels,:); 
    dataCounterclockwise{i}.yPanelVolts =  counterclockwiseTrials{i}(yPanels,:);
    dataCounterclockwise{i}.ficTracAngularPosition = counterclockwiseTrials{i}(headingFly,:); 
    dataCounterclockwise{i}.ficTracIntx = -counterclockwiseTrials{i}(xFly,:); 
    dataCounterclockwise{i}.ficTracInty = counterclockwiseTrials{i}(yFly,:); 
end
    
    
%% Downsample, unwrap and smooth position data, then get velocity and smooth

sizeBall = 9;

for i = 1:length(closedloopTrials)
    [smoothedSF{i}] = posDataDecoding(dataSF{1,i},1000);
end

for i = 1:length(clockwiseTrials)
    [smoothedClockwise{i}] = posDataDecoding(dataClockwise{1,i},1000);
end

for i = 1:length(counterclockwiseTrials)
    [smoothedCounterclockwise{i}] = posDataDecoding(dataCounterclockwise{1,i},1000);
end

% for most experiments I have used 1000 Hz as a sample rate, and it is what
% I will use from now on, so that's how I'll leave it, but this could be
% changed in case of need

%get full data smoothed without sorting
data.xpanelVolts =  rawData (:,xPanels); 
data.yPanelVolts =  rawData (:, yPanels);
data.ficTracAngularPosition = rawData ( : , headingFly); 
data.ficTracIntx = -rawData ( : , xFly); 
data.ficTracInty = rawData ( : , yFly); 

[smoothed] = posDataDecoding(data,1000);
%% Forward velocity analysis

%(1) Stripe fixation trials
for i = 1:length(smoothedSF)
    forwardVelocity{i} = smoothedSF{1,i}.xVel(1:750);
end

forwardVelocity = cell2mat(forwardVelocity);
forwardVelocity = reshape(forwardVelocity,[750,length(forwardVelocity)/750]);

meanVelocity = mean(forwardVelocity,2);
time = linspace(0,30,length(forwardVelocity));

figure('Position', [100 100 1400 900]),
subplot(1,3,1)
plot(time,forwardVelocity)
hold on
plot(time,meanVelocity,'k','lineWidth',2)
title('Forward velocity for closed-loop trials')
xlabel('Time (sec)'); ylabel('Velocity (mm/s)');
xlim([0 8]); ylim([min(min(forwardVelocity))-2, max(max(forwardVelocity))+2]);


%(2) Clockwise optomotor trials
for i = 1:length(smoothedClockwise)
    forwardVelocityClockwise{i} = smoothedClockwise{i}.xVel(1:72);
end

forwardVelocityClockwise = cell2mat(forwardVelocityClockwise);
forwardVelocityClockwise = reshape(forwardVelocityClockwise,[72,length(forwardVelocityClockwise)/72]);

meanVelocityClockwise = mean(forwardVelocityClockwise,2);
time = linspace(0,3,72);

subplot(1,3,2)
plot(time,forwardVelocityClockwise)
hold on
plot(time,meanVelocityClockwise,'k','lineWidth',2)
title('Forward velocity for clockwise optomotor trials')
xlabel('Time (sec)');
xlim([0 3]); ylim([min(min(forwardVelocity))-2, max(max(forwardVelocity))+2]);


%(3) Counterclockwise optomotor trials
for i = 1:length(smoothedCounterclockwise)
    forwardVelocityCounterclockwise{i} = smoothedCounterclockwise{i}.xVel(1:72);
end

forwardVelocityCounterclockwise = cell2mat(forwardVelocityCounterclockwise);
forwardVelocityCounterclockwise = reshape(forwardVelocityCounterclockwise,[72,length(forwardVelocityCounterclockwise)/72]);

meanVelocityCounterclockwise = mean(forwardVelocityCounterclockwise,2);
time = linspace(0,3,72);

subplot(1,3,3)
plot(time,forwardVelocityCounterclockwise)
hold on
plot(time,meanVelocityCounterclockwise,'k','lineWidth',2)
title('Forward velocity for counterclockwise optomotor trials')
xlabel('Time (sec)');
xlim([0 3]); ylim([min(min(forwardVelocity))-2, max(max(forwardVelocity))+2]);

saveas(gcf,strcat(path,'ForwardVelocity',file(1:end-4),'.png'));
%% Angular velocity

%(1) Closed loop trials
for i = 1:length(smoothedSF)
    angularVelocity{i} = smoothedSF{i}.angularVel(1:750);
end

angularVelocity = cell2mat(angularVelocity);
angularVelocity = reshape(angularVelocity,[750,length(angularVelocity)/750]);

meanVelocity = mean(angularVelocity,2);
time = linspace(0,30,length(angularVelocity));

figure('Position', [100 100 1400 900]),
subplot(1,3,1)
plot(time,angularVelocity)
hold on
plot(time,meanVelocity,'k','lineWidth',2)
title('Angular velocity for closed loop trials')
xlabel('Time (sec)'); ylabel('Velocity (deg/s)');
xlim([0 8]); ylim([min(min(angularVelocity))-2, max(max(angularVelocity))+2]);


%(2) Clockwise optomotor trials
for i = 1:length(smoothedClockwise)
    angularVelocityClockwise{i} = smoothedClockwise{i}.angularVel(1:72);
end

angularVelocityClockwise = cell2mat(angularVelocityClockwise);
angularVelocityClockwise = reshape(angularVelocityClockwise,[72,length(angularVelocityClockwise)/72]);

meanVelocityClockwise = mean(angularVelocityClockwise,2);

%(3) Counterclockwise optomotor trials
for i = 1:length(smoothedCounterclockwise)
    angularVelocityCounterclockwise{i} = smoothedCounterclockwise{i}.angularVel(1:72);
end

angularVelocityCounterclockwise = cell2mat(angularVelocityCounterclockwise);
angularVelocityCounterclockwise = reshape(angularVelocityCounterclockwise,[72,length(angularVelocityCounterclockwise)/72]);

meanVelocityCounterclockwise = mean(angularVelocityCounterclockwise,2);
time = linspace(0,3,72);

subplot(1,3,2)
plot(time,angularVelocityClockwise)
hold on
plot(time,meanVelocityClockwise,'k','lineWidth',2)
title('Angular velocity for clockwise optomotor trials')
xlabel('Time (sec)');
xlim([0 3]); ylim([min(min(angularVelocityClockwise))-2, max(max(angularVelocityCounterclockwise))+2]);

subplot(1,3,3)
plot(time,angularVelocityCounterclockwise)
hold on
plot(time,meanVelocityCounterclockwise,'k','lineWidth',2)
title('Angular velocity for counterclockwise optomotor trials')
xlabel('Time (sec)');
xlim([0 3]); ylim([min(min(angularVelocityClockwise))-2, max(max(angularVelocityCounterclockwise))+2]);

saveas(gcf,strcat(path,'AngularVelocity',file(1:end-4),'.png'));


%Look at their evolution with a heatmap.
figure('Position', [100 100 1000 900]),
colormap(hot)
subplot(1,2,1),
imagesc(time,[1:length(clockwiseTrials)],angularVelocityClockwise')
colorbar
title('Angular Velocity in clockwise trials');

subplot(1,2,2),
imagesc(time,[1:length(counterclockwiseTrials)],angularVelocityCounterclockwise')
colorbar
title('Angular Velocity in counterclockwise trials');

save([path,'AngularVelocityClockwise',file(1:end-4),'.mat'],'angularVelocityClockwise');
save([path,'AngularVelocityCounterclockwise',file(1:end-4),'.mat'],'angularVelocityCounterclockwise');
saveas(gcf,strcat(path,'AngularVelocityEvolution',file(1:end-4),'.png'));


%% 2D trajectories

%closed-loop bouts
figure('Position', [100 100 1400 900]),
for i = 1:length(smoothedSF)
subplot(10,5,i)
[posx{i},posy{i}]=FlyTrajectory(smoothedSF{1,i}.Intx',smoothedSF{1,i}.Inty',smoothedSF{1,i}.angularPosition');
timeTraj{i}= linspace(1,1000,length(posx{i}));
plot(posx{i},posy{i},'k')
axis equal
axis tight
end
suptitle('2D trajectories in individual closed-loop bouts');

saveas(gcf,strcat(path,'ClosedLoopTrajectories',file(1:end-4),'.png'));
%% Heading during closed-loop trials

pxToDeg = 360/96;

figure('Position', [100 100 1600 900]),

for i = 1:length(smoothedSF)
 flyPosPx{i} = rad2deg(smoothedSF{1,i}.angularPosition)/pxToDeg;
 FlyPosDeg{i} = zeros(1,length(smoothedSF{1,i}.angularPosition));

% Convert from xpos to degrees, knowing that xpos 70 = 0 deg
    for j = 1:length(smoothedSF{1,i}.angularPosition)
        if flyPosPx{i}(j) == 70
            FlyPosDeg{i}(j) = 0;
        elseif flyPosPx{i}(j) >70 
            FlyPosDeg{i}(j) = (flyPosPx{i}(j)-70)*pxToDeg; % Correct the offset and multiply by factor to get deg
        else
            FlyPosDeg{i}(j) = (flyPosPx{i}(j)+27)*pxToDeg; % Correct the offset and multiply by factor to get deg
        end
    end

    FlyPos360{i} = wrapTo360(FlyPosDeg{i});
    flyPos180{i} = wrapTo180(FlyPos360{i});

% In polar coordinates...
% Taking every frame into account
    subplot(10,5,i)
    posToRadFly{i} = deg2rad(FlyPos360{i});
    circedges = [0:20:360];
    circedges = deg2rad(circedges);
    polarhistogram(posToRadFly{i},circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
    ax = gca;
    ax.ThetaDir='clockwise';
    ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

end

suptitle('Polar distribution of heading for closed-loop trials');
saveas(gcf,strcat(path,'FlyPositionClosedLoop',file(1:end-4),'.png'))


figure('Position', [100 100 1600 900]),

for i = 1:length(smoothedSF)

%with only moving frames
    subplot(10,5,i)
    posToRadFlyMoving{i} = posToRadFly{i}(smoothedSF{1,i}.xVel>1);
    polarhistogram(posToRadFlyMoving{i},circedges,'Normalization','probability','FaceColor',[0,0,0],'HandleVisibility','off');
    ax = gca;
    ax.ThetaDir='clockwise';
    ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

end

suptitle('Polar distribution of heading for closed-loop trials; moving frames only');
saveas(gcf,strcat(path,'FlyPositionClosedLoopMovingFrames',file(1:end-4),'.png'))

%% Optomotor response analysis


%Import the empty trial data, take random 3 sec snippets of it and smooth
darknessData = load(strcat(path,'EmptyTrialExpNum1'));
darknessData = darknessData.rawData;

%take random 2932 snippets of time (length of opto trials)
sample = round(rand(1,60)*(length(darknessData)-2931));
for i = 1:length(sample)
    darkness{i} = darknessData(sample(i):sample(i)+2931,:);
    dataDarkness{i}.xpanelVolts =  darkness{1,i}(:,xPanels); 
    dataDarkness{i}.yPanelVolts =  darkness{1,i}(:,yPanels);
    dataDarkness{i}.ficTracAngularPosition = darkness{1,i}(:,headingFly); 
    dataDarkness{i}.ficTracIntx = -darkness{1,i}(:,xFly); 
    dataDarkness{i}.ficTracInty = darkness{1,i}(:,yFly); 
    [smoothedDark{i}] = posDataDecoding(dataDarkness{1,i},1000);
    angularVelocityDark{i} = smoothedDark{i}.angularVel(1:72);
end

angularVelocityDark = cell2mat(angularVelocityDark);


%(1) Take different windows of time tau
%each opto response trial has 72 frames after smoothing. I can take windows
%of about ten frames from 0 to 72, or take increasingly long windows?

%For clockwise trials 

%I'll start testing 7 windows of 1 frames each
figure('Position', [100 100 1600 900]),
for tau = 1:6

%(2) Detect the maximum slope of each curve during that window
    for i = 1:size(angularVelocityClockwise,2)
        %get 6 11 frame windows for each trial
        PartialCurve{tau,i} = angularVelocityClockwise(tau*10:tau*10+10,i);
        %take the derivative
        diffPartialCurve{tau,i} = diff(PartialCurve{tau,i});
        %get the maximum of the derivative and that will be the maximum
        %slope
        maxDiffPartialCurve(tau,i) = max(diffPartialCurve{tau,i});
        
%(4) Repeat 1-3 for randomly sampled 3s snipped of the empty trial
        PartialCurveDark{tau,i} = angularVelocityDark(tau*10:tau*10+10,i);
        diffPartialCurveDark{tau,i} = diff(PartialCurveDark{tau,i});
        maxDiffPartialCurveDark(tau,i) = max(diffPartialCurveDark{tau,i});

    end
%(5) Compute and plot 2 histograms for the max slopes for the exp vs
%darkness
    subplot(2,3,tau)
    histogram(maxDiffPartialCurve(tau,:),10)
    hold on
    histogram(maxDiffPartialCurveDark(tau,:),10,'FaceColor','red')
    %d(tau) = pdist2(maxDiffPartialCurve(tau,:),maxDiffPartialCurveDark(tau,:)); %compute the distance between the 2 types of data
    %[mC,iC] = max(d); 
    legend('Trials','Darkness');
%(6) Measure the d'
    discrim(tau) =((mean(maxDiffPartialCurve(tau,:)))-(mean(maxDiffPartialCurveDark(tau,:))))/(sqrt((var(maxDiffPartialCurve(tau,:))+var(maxDiffPartialCurveDark(tau,:)))/2));
end
suptitle('Distribution of maximum slopes for clockwise trials');
saveas(gcf,strcat(path,'ClockwiseTaus',file(1:end-4),'.png'))

[mC,iC] = max(discrim);
%% 
%(7) Find the best window (i.e., the one with the biggest d')

%say that based on this, the perfect window is the one that yields the
%maximum distance d.
%I will substract the response of each trial to an average baseline
%activity in the first 10 frames
baseline = mean(mean(angularVelocityClockwise(1:10,:)));
response = mean(angularVelocityClockwise(10*iC+1:10*iC+10,:))-baseline;

%plot evolution of the response in time
figure, plot(response,'ko')
xlabel('Trial number');ylabel('Response magnitude')
xlim([0, size(angularVelocityClockwise,2)]);
hold on
plot([0, size(angularVelocityClockwise,2)],[0, 0],'-.r');
title('Evolution of clockwise optomotor trials');

saveas(gcf,strcat(path,'ClockwiseResponseEv',file(1:end-4),'.png'))


%% 

%counterclockwise trials
figure('Position', [100 100 1600 900]),
for tau = 1:6

%(2) Detect the maximum slope of each curve during that window
    for i = 1:size(angularVelocityCounterclockwise,2)
        PartialCurveCC{tau,i} = angularVelocityCounterclockwise(tau*10:tau*10+10,i);
        diffPartialCurveCC{tau,i} = diff(PartialCurveCC{tau,i});
        minDiffPartialCurveCC(tau,i) = min(diffPartialCurveCC{tau,i});
        PartialCurveDarkCC{tau,i} = angularVelocityDark(tau*10:tau*10+10,i);
        diffPartialCurveDarkCC{tau,i} = diff(PartialCurveDarkCC{tau,i});
        minDiffPartialCurveDark(tau,i) = min(diffPartialCurveDarkCC{tau,i});

    end
%(5) Compute and plot 2 histograms for the max slopes for the exp vs
%darkness
    subplot(2,3,tau)
    histogram(minDiffPartialCurveCC(tau,:),10)
    hold on
    histogram(minDiffPartialCurveDark(tau,:),10,'FaceColor','red')
    %dCC(tau) = pdist2(minDiffPartialCurveCC(tau,:),minDiffPartialCurveDark(tau,:)); %compute the distance between the 2 types of data
    %[mCC,iCC] = max(dCC); %likewise here, I want to make sure that the exp curves are maximally shifted to the LEFT
    legend('Trials','Darkness');
%(6) Measure the d'

    discrimCC(tau) =((mean(minDiffPartialCurveCC(tau,:)))-(mean(minDiffPartialCurveDark(tau,:))))/(sqrt((var(minDiffPartialCurveCC(tau,:))+var(minDiffPartialCurveDark(tau,:)))/2));

end
suptitle('Distribution of minimum slopes for counterclockwise trials');
saveas(gcf,strcat(path,'CounterclockwiseTaus',file(1:end-4),'.png'))

%% 

[mCC,iCC] = min(discrimCC);

%say that based on this, the perfect window is the one that yields the
%maximum distance dCC.
%I will substract the response of each trial to an average baseline
%activity in the first 10 frames
baselineCC = mean(mean(angularVelocityCounterclockwise(1:10,:)));
responseCC = baselineCC-mean(angularVelocityCounterclockwise(10*iCC+1:10*iCC+10,:));

%plot evolution of the response in time
figure, plot(responseCC,'ko')
xlabel('Trial number');ylabel('Response magnitude')
xlim([0, size(angularVelocityCounterclockwise,2)]);
hold on
plot([0, size(angularVelocityCounterclockwise,2)],[0, 0],'-.r');
title('Evolution of counterclockwise optomotor trials');
saveas(gcf,strcat(path,'CounterclockwiseResponseEv',file(1:end-4),'.png'))


%save both responses 
save([path,'Responses',file(1:end-4),'.mat'],'response','responseCC');


%Plotting them both together taking into account the actual trial number in
%the full number of trials:

clockTrials = jumpFunction==52;
counterclockTrials = jumpFunction==53;

allTrials(clockTrials) = response;
allTrials(counterclockTrials) = NaN;

allTrialsCC(clockTrials) = NaN;
allTrialsCC(counterclockTrials) = responseCC;

figure,
plot(allTrials,'bo')
hold on
plot(allTrialsCC,'ko')
title('Response evolution for all optomotor trials')
xlabel('Trial number');ylabel('Response magnitude')
xlim([0, size(allTrials,2)]);
hold on
plot([0, size(allTrials,2)],[0, 0],'-.r');
legend('Clockwise trials','Counterclockwise trials');
saveas(gcf,strcat(path,'AllResponsesEv',file(1:end-4),'.png'))


%% PCA with the opto trials

%Compute the PCAs
[coeffClock,scoresClock] = pca(angularVelocityClockwise');
[coeffCounterclock,scoresCounterclock] = pca(angularVelocityCounterclockwise');

%Plot the biplots
figure,
subplot(1,2,1)
biplot(coeffClock(:,1:2),'Scores',scoresClock(:,1:2));
title('Biplot for clockwise trials');
subplot(1,2,2)
biplot(coeffCounterclock(:,1:2),'Scores',scoresCounterclock(:,1:2));
title('Biplot for counterclockwise trials');
saveas(gcf,strcat(path,'PCA',file(1:end-4),'.svg'))


%cross correlation between PC1 and the responses
for i = 1:size(angularVelocityClockwise,2)
    [crossCorrClockwise{i},lagsClockwise{i}] = xcorr(coeffClock(:,1),angularVelocityClockwise(:,i));
    figure,
    stem(lagsClockwise{i},crossCorrClockwise{i})
end

%  
%  
% 
% %% add trajectory for the full experiment, and plot in 3 different colors the 3 different trial types
% 
% 
% [posx,posy]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.angularPosition);
% 
% 
% figure, 
% %plot the trajectories for each trial type in a different color
% plot(posx,posy,'r')
% axis equal
% axis tight
% title('2D trajectory of the fly');
% 
