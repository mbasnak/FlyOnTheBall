%Experiment 7, pooled blocks analysis

%This code pools the data from the various blocks of an individual fly and
%analyses it.

clear all; close all;

order = [1,3,2];
%% Dist2goal

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment7'
[path] = uigetdir();

files = dir(strcat(path,'\goals*.mat'));

for i = 1:length(files)
data{i} = load(strcat(files(i).folder,'\',files(i).name));
blockName{i} = files(i).name(6:end-4);
end

data = data(order);
blockName = blockName(order);

%plot distance to the goal for all three blocks on the same plot
figure,
for i = 1:length(data)
plot(data{1,i}.degsFlyDistMoving,data{1,i}.probabilitiesDistMoving)
hold on
end
title('Distance to the goal with moving frames');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');
legend(blockName{1},blockName{2},blockName{3});

saveas(gcf,strcat(files(1).folder,'\Dist2GoalAllBlocks.png'))
%% Activity

files = dir(strcat(path,'\Activity*.mat'));

for i = 1:length(files)
data{i} = load(strcat(files(i).folder,'\',files(i).name));
blockName{i} = files(i).name(9:end-4);
end
data = data(order);
blockName = blockName(order);

figure,
newMap = flipud(gray);

trials = [0:1];
colormap(newMap)

for i = 1:length(data)
    subplot(length(data),1,i)
    xaxis{i} = data{1,i}.time;
    imagesc(xaxis{i},trials,data{1,i}.activity')
    ylabel('Activity');
    xlabel('Time (s)');
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    title(strcat('ActivityRP',blockName{i},'percentage activity:', num2str(data{1,i}.percentageActivity), '%'));
end

saveas(gcf,strcat(files(1).folder,'\ActivityRPAllBlocks.png'))

%% Fly's heading

files = dir(strcat(path,'\FlyPosition*.mat'));

for i = 1:length(files)
data{i} = load(strcat(files(i).folder,'\',files(i).name));
blockName{i} = files(i).name(12:end-4);
end
data = data(order);
blockName = blockName(order);

figure('Position', [100 100 1600 900]),
for i = 1:length(data)
    subplot(1,length(data),i)
    polarhistogram(data{1,i}.posToRadFlyMoving,data{1,i}.circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
    ax = gca;
    ax.ThetaDir='clockwise';
    ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
    hold on
    title(strcat('FlyHeading',blockName{i}));
end

saveas(gcf,strcat(files(1).folder,'\FlyHeadingAllBlocks.png'))

%% Around jump velocities

files = dir(strcat(path,'\perTrialData*.mat'));

for i = 1:length(files)
data{i} = load(strcat(files(i).folder,'\',files(i).name));
blockName{i} = files(i).name(13:end-4);
end
data = data(order);
blockName = blockName(order);

%classify data into '90' and '-90' trials
for i = 1:length(files)   
Data90(i).forwardVel = data{1,i}.perTrialData.forwardVel(:,data{1,i}.perTrialData.jumpMag == 90);
Data90(i).angVel = data{1,i}.perTrialData.angVel(:,data{1,i}.perTrialData.jumpMag == 90);
DataNeg90(i).forwardVel = data{1,i}.perTrialData.forwardVel(:,data{1,i}.perTrialData.jumpMag == -90);
DataNeg90(i).angVel = data{1,i}.perTrialData.angVel(:,data{1,i}.perTrialData.jumpMag == -90);
end

%get mean forward and angular velocity per group
for i = 1:length(files)
meanForwardVel90{i} = mean(Data90(i).forwardVel,2);
meanForwardVelNeg90{i} = mean(DataNeg90(i).forwardVel,2);
meanAngVel90{i} = mean(Data90(i).angVel,2);
meanAngVelNeg90{i} = mean(DataNeg90(i).angVel,2);
end

%Plot
time = linspace(-10,10,length(meanAngVel90{1,1}));

figure('Position', [100 100 1600 900]),
for i = 1:length(files)
subplot(1,length(data),i)
plot(time,meanAngVel90{i},'r')
hold on
plot(time,meanAngVelNeg90{i},'k')
xlim([-5,5]); ylim([-80,80]);
title(strcat('MeanAngVel',blockName{i}));
legend({'90','-90'});
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump');
plot([0,0],[-80, 80],'k','HandleVisibility','off');
plot([-5,5],[0,0],'-.k','HandleVisibility','off');
end

saveas(gcf,strcat(files(1).folder,'\AJangvelAllBlocks.png'))

%Use heatmaps to see the evolution
figure,
for i = 1:length(files)
    subplot(1,length(files),i)
    trials = size(DataNeg90(i).angVel,2);
    imagesc(time,trials,DataNeg90(i).angVel')
    colorbar
    title(strcat('AngVel ',blockName{i}))
end

%% Dist2Goal 10 sec around jumps

files = dir(strcat(path,'\shortData10sec*.mat'));

for i = 1:length(files)
    data{i} = load(strcat(files(i).folder,'\',files(i).name));
    blockName{i} = files(i).name(15:end-4);
end
data = data(order);
blockName = blockName(order);

figure('Position', [100 100 1600 900]),
for i = 1:length(files)
    subplot(1,length(files),i)
    newMap = flipud(gray);
    xaxis = [-180:360/17:180];
    trials = [1:48];
    imagesc(xaxis,trials,data{1,i}.probaDist10secMoving')
    colormap(newMap)
    colorbar
    title(blockName{i})
    xlabel('Dist2goal 10 sec around jump (deg)')
end

saveas(gcf,strcat(files(1).folder,'\Dist2Goal10secAllBlocks.png'))
