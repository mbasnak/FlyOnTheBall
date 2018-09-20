%% Analysis of the panel/FicTrac data

%%% This code analyses the output of the panels and the FicTrac data

%I'm going to add something that loads every dataExpNum.mat in the folder
%and then at the bottom saves the figures

% Finding every file we're interested to load
files = dir('dataExpNum*.mat'); %(the ones of the form datExpNum
% load every file's rawData into a cell array
allData = {};
for i = 1:size(files,1)
 allData{i} = load(files(i).name);
end

% Define Ni-Daq channels ID
headingFly = 1;
xFly = 2;
yFly = 3;
xPanels = 4;
yPanels = 5;

%% Subset acquisition of x and y pos, as well as FicTrac data

%Panels data
VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
maxValY = 1;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2


for i = 1:size(files,1)
%Panel data
data.xPanelVolts{i} =  allData{1,i}.rawData (:,xPanels);
data.xPanelPos{i} = round ((data.xPanelVolts{i}  * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels
data.yPanelVolts{i} =  allData{1,i}.rawData (:,yPanels);
data.yPanelPos{i} = round ((data.yPanelVolts{i}  * maxValY) /VOLTAGE_RANGE);

%FicTrac data
data.ficTracAngularPosition{i} = allData{1,i}.rawData (:,headingFly); 
data.ficTracIntx{i} = allData{1,i}.rawData (:,xFly);
data.ficTracInty{i} = allData{1,i}.rawData (:,yFly);
end


%% 
% 
% data.xPanelVolts =  rawData (:,xPanels); 
% VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
% maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
% data.xPanelPos = round ((data.xPanelVolts  * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels
% 
% data.yPanelVolts =  rawData (:, yPanels);
% VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
% maxValY = 1;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2
% data.yPanelPos = round ((data.yPanelVolts  * maxValY) /VOLTAGE_RANGE);
% 
% %FicTrac data
% data.ficTracAngularPosition = rawData ( : , headingFly); 
% data.ficTracIntx = rawData ( : , xFly); 
% data.ficTracInty = rawData ( : , yFly); 
%% 

% Pos x=5 is 0 deg (ie facing the fly), I measured this empirically

pxToDeg = 360/97; % There are 97 possible positions (the last one = first one) and this represents 360 deg
posToDeg = {};

% Convert from xpos to degrees, knowing that xpos 5 = 0 deg
for j = 1:size(data.xPanelVolts,2)
    for i=1:length(data.xPanelPos{1,j})
        if data.xPanelPos{1,j}(i) ==1 | data.xPanelPos{1,j}(i) ==2 | data.xPanelPos{1,j}(i) ==3 | data.xPanelPos{1,j}(i) ==4
        posToDeg{j,i} = (data.xPanelPos{1,j}(i)+93)*pxToDeg; % Correct the offset and multiply by factor to get deg
        else
        posToDeg{j,i} = (data.xPanelPos{1,j}(i)-5)*pxToDeg; % Correct the offset and multiply by factor to get deg
        end       
    end
            % Create a time vector in sec
        time{j} = linspace(0,(size(allData{1,j}.rawData,1)/1000),size(allData{1,j}.rawData,1)); %1000 is our sampling rate
        posTo = posToDeg{j,:};
        posTo = cell2mat(posTo);
        
        % Plot the position of the stimulus in degrees for the trial
        figure,
        plot(time{j},posTo)
        ylim([0 365]);
        ylabel('Position of the stimulus (deg)');
        xlabel('Time (s)');
        title('Position of the stimulus in deg as a function of time');
end

%this is currently not working.

%% Probability density function of the stimulus position

% Remapping the positions to span -180 to 180 deg

%for e = 1:length(ExpNum)
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
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Probability density of the stimulus position');
ylabel('Probability density'); xlabel('Stimulus position (deg)');

%saveas(gcf,strc('ProbabilityOfStimulusPositionExpNum', e, '.png'))
%end

%% Using Yvette's function to filtrate and unwrap the data

LOWPASS_FILTER_CUTOFF= 100; % Hz
THRESHOLD_ANGULAR_VELOCITY = 900; % degrees / s  this is the max velocity that can be allowed into analysis
% decode angular velocity and accumulated position
[ angularVelocity , filteredFicTracPos ] = ficTracSignalDecoding(data.ficTracAngularPosition, 1000, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);

figure,
subplot(3,1,1)
plot(time,posToDeg)
ylim([0 365]);
ylabel('Position of the stimulus in degrees');
xlabel('Time (s)');
title('Position of the stimulus in degrees');

subplot(3,1,2)
plot(time,angularVelocity)
ylabel('Angular velocity of the fly');
xlabel('Time (s)');
title('Angular velocity of the fly');

subplot(3,1,3)
plot(time,filteredFicTracPos)
ylabel('Heading of the fly? (deg)');
xlabel('Time (s)');
ylim([0 365]);
title('Heading of the fly?');
%For some reason this last plot has different time units.