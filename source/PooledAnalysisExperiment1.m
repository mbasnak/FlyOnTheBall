% Analysis for pooled flies, experiment 1

close all; clear all;

%% Load all the data for an individual fly

% prompt the user to select the file to open and load it.
% choose directory of cell to be analyzed
FlyPath = uigetdir()
% set the current path to that directory
path = cd(FlyPath);
files = dir(path);
fileNames = {};

for i = 1:size(files,1)
    fileNames{i} = files(i).name;
end

names = find(contains(fileNames,'data'));

flyData = {};
for i=1:size(names,2)  
    name(i) = strcat(path,'\',fileNames(names(i)));
    flyData{i} = load(name{i});
end

for i = 1:length(flyData)
    Data{i} = flyData{1,i}.rawData;
end

%% % Calculate the mean heading angle of each trial

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;

% Subset acquisition of FicTrac data
for i = 1:length(Data)
    ficTracAngularPosition{i} = Data{i}( : , headingFly);
end

% Process the data
for i=1:length(ficTracAngularPosition)
    downsampledAngularPosition{i} = downsample(ficTracAngularPosition{i},1000/25);
    downsRadAngularPosition{i} = downsampledAngularPosition{i} .* 2 .* pi ./ 10;
    unwrappedAngularPosition{i} = unwrap(downsRadAngularPosition{i});
    smoothedAngularPosition{i} = smoothdata(unwrappedAngularPosition{i},10);
    degAngularPosition{i} = (smoothedAngularPosition{i} / (2*pi)) * 360;
    degAngularPosition{i} = wrapTo360(degAngularPosition{i});
    angPos{i} = deg2rad(degAngularPosition{i});
end

% Calculate the mean of means and name it the "goal"

angPos = cell2mat(degAngularPosition);
for i = 1:size(angPos,2)
    CircularStatsFly = circ_stats(angPos(:,i));
    circMean(i) = CircularStatsFly.mean;
    circStd(i) = CircularStatsFly.std;
end


% Plot the mean heading with error for every trial
figure,
bar(rad2deg(circMean))
hold on
errorbar(rad2deg(circMean),rad2deg(circStd),'.')
ylim([-180 180]);
ylabel('Angle (deg)'); xlabel('Trial #');
title('Mean heading across trials');

%% Establish distance to goal

goal = mean(circMean);
goalDeg = rad2deg(goal);

dist2goal = wrapTo180(angPos - goal);
save('distanceToGoal.mat','dist2goal','goalDeg');

figure, histogram(dist2goal,'Normalization','probability');
