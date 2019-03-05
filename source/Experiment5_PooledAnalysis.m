%Experiment 5, pooled analysis

%This code pools the data from all the flies in experiment 5
%and does some pooled analyses.

clear all; close all;


%% Around jump velocities

FirstLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\perTrialDataLowErrorBlockExp3.mat');
HighErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\perTrialDataHighErrorBlockExp4.mat');
SecondLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\perTrialDataLowErrorBlockExp5.mat');

%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from.
for i = 1:length(FirstLowErrorBlock)
    dataFirstLE{i} = load(strcat(FirstLowErrorBlock(i).folder,'\',FirstLowErrorBlock(i).name));
    dataHE{i} = load(strcat(HighErrorBlock(i).folder,'\',HighErrorBlock(i).name));
    dataSecondLE{i} = load(strcat(SecondLowErrorBlock(i).folder,'\',SecondLowErrorBlock(i).name));
end

% First low error block
%concatenate as columns the arrays for angular velocity, forward
%velocity and jump magnitude for every fly
angVelAll = zeros(0,0);
jumpMagAll = zeros(0,0);

for i = 1:length(dataFirstLE)
    angVelAll = [angVelAll,dataFirstLE{1,i}.perTrialData.angVel];
    jumpMagAll = [jumpMagAll,dataFirstLE{1,i}.perTrialData.jumpMag];
end

%subset into the 2 magnitude groups
Data90.angVel = angVelAll(:,jumpMagAll == 90);
DataNeg90.angVel = angVelAll(:,jumpMagAll == -90);

%calculate mean and std ang vel
meanAngVel90 = mean(Data90.angVel,2);
stdAngVel90 = std(Data90.angVel,[],2);
meanAngVelNeg90 = mean(DataNeg90.angVel,2);
stdAngVelNeg90 = std(DataNeg90.angVel,[],2);

%plot mean and error for each
time = linspace(-15,15,length(angVelAll));

figure('Position', [100 100 1600 900]),
subplot(1,3,1)
p1 = boundedline(time,meanAngVel90,stdAngVel90/sqrt(i),'r','alpha')
hold on
p2 = boundedline(time,meanAngVelNeg90,stdAngVelNeg90/sqrt(i),'k','alpha')
title('First low error block');
ylim([-80, 80]);xlim([-8, 8]);
plot([0,0],[-80, 80],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p1,p2], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');


% High error block
%concatenate as columns the arrays for angular velocity, forward
%velocity and jump magnitude for every fly
angVelAllHE = zeros(0,0);
jumpMagAllHE = zeros(0,0);

for i = 1:length(dataFirstLE)
    angVelAllHE = [angVelAllHE,dataHE{1,i}.perTrialData.angVel];
    jumpMagAllHE = [jumpMagAllHE,dataHE{1,i}.perTrialData.jumpMag];
end

%subset into the 2 magnitude groups
Data90HE.angVel = angVelAllHE(:,jumpMagAllHE == 90);
DataNeg90HE.angVel = angVelAllHE(:,jumpMagAllHE == -90);

%calculate mean and std ang vel
meanAngVel90HE = mean(Data90HE.angVel,2);
stdAngVel90HE = std(Data90HE.angVel,[],2);
meanAngVelNeg90HE = mean(DataNeg90HE.angVel,2);
stdAngVelNeg90HE = std(DataNeg90HE.angVel,[],2);

%plot mean and error for each
subplot(1,3,2)
p21 = boundedline(time,meanAngVel90HE,stdAngVel90HE/sqrt(i),'r','alpha')
hold on
p22 = boundedline(time,meanAngVelNeg90HE,stdAngVelNeg90HE/sqrt(i),'k','alpha')
title('High error block');
ylim([-80, 80]);xlim([-8, 8]);
plot([0,0],[-80, 80],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p21,p22], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');


% Second low block
%concatenate as columns the arrays for angular velocity, forward
%velocity and jump magnitude for every fly
angVelAllSecondLE = zeros(0,0);
jumpMagAllSecondLE = zeros(0,0);

for i = 1:length(dataFirstLE)
    angVelAllSecondLE = [angVelAllSecondLE,dataSecondLE{1,i}.perTrialData.angVel];
    jumpMagAllSecondLE = [jumpMagAllSecondLE,dataSecondLE{1,i}.perTrialData.jumpMag];
end

%subset into the 2 magnitude groups
Data90SecondLE.angVel = angVelAllSecondLE(:,jumpMagAllSecondLE == 90);
DataNeg90SecondLE.angVel = angVelAllSecondLE(:,jumpMagAllSecondLE == -90);

%calculate mean and std ang vel
meanAngVel90SecondLE = mean(Data90SecondLE.angVel,2);
stdAngVel90SecondLE = std(Data90SecondLE.angVel,[],2);
meanAngVelNeg90SecondLE = mean(DataNeg90SecondLE.angVel,2);
stdAngVelNeg90SecondLE = std(DataNeg90SecondLE.angVel,[],2);

%plot mean and error for each
subplot(1,3,3)
p31 = boundedline(time,meanAngVel90SecondLE,stdAngVel90SecondLE/sqrt(i),'r','alpha')
hold on
p32 = boundedline(time,meanAngVelNeg90SecondLE,stdAngVelNeg90SecondLE/sqrt(i),'k','alpha')
title('Second low error block');
ylim([-80, 80]);xlim([-8, 8]);
plot([0,0],[-80, 80],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p31,p32], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\AverageVelocityChange.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\AverageVelocityChange.svg'))


%% distance to the goal

FilesFirstLE = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\goalsLowErrorBlockExp3.mat');
FilesHE = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\goalsHighErrorBlockExp4.mat');
FilesSecondLE = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\goalsLowErrorBlockExp5.mat');

%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(FilesHE)
    goalsFirstLE{i} = load(strcat(FilesFirstLE(i).folder,'\',FilesFirstLE(i).name));
    goalsHE{i} = load(strcat(FilesHE(i).folder,'\',FilesHE(i).name));
    goalsSecondLE{i} = load(strcat(FilesSecondLE(i).folder,'\',FilesSecondLE(i).name));
end

%concatenate as columns the arrays for all the variables
goalMovingAll = zeros(0,0);
dist2goalMovingAll = {};
goalMovingAllHE = zeros(0,0);
dist2goalMovingAllHE = {};
goalMovingAllSecondLE = zeros(0,0);
dist2goalMovingAllSecondLE = {};

for i = 1:length(FilesHE)
    goalMovingAll = [goalMovingAll;goalsFirstLE{1,i}.goalMoving];
    dist2goalMovingAll{i} = goalsFirstLE{1,i}.dist2goalMoving;
    goalMovingAllHE = [goalMovingAllHE;goalsHE{1,i}.goalMoving];
    dist2goalMovingAllHE{i} = goalsHE{1,i}.dist2goalMoving;
    goalMovingAllSecondLE = [goalMovingAllSecondLE;goalsSecondLE{1,i}.goalMoving];
    dist2goalMovingAllSecondLE{i} = goalsSecondLE{1,i}.dist2goalMoving;
end

%Look at goal distribution per block
figure('Position', [100 100 1600 900]),
subplot(2,3,1)
plot(rad2deg(goalMovingAll),'ro')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('First low error block');
subplot(2,3,4)
polarhistogram(goalMovingAll,10,'FaceColor','red')
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top';

subplot(2,3,2)
plot(rad2deg(goalMovingAllHE),'ro')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('High error block');
subplot(2,3,5)
polarhistogram(goalMovingAllHE,10,'FaceColor','red')
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top';

subplot(2,3,3)
plot(rad2deg(goalMovingAllSecondLE),'ro')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('Second low error block');
subplot(2,3,6)
polarhistogram(goalMovingAllSecondLE,10,'FaceColor','red')
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top';


%Looks at the distribution of the variable 'Distance to the goal' per block
figure('Position', [100 100 1600 900]),
edges = [-180:20:180];
%For moving frames only
for i = 1:length(dist2goalMovingAll)
    [countsDistMoving{i}] = histcounts(dist2goalMovingAll{i},edges);
    probabilitiesDistMoving{i} = countsDistMoving{i}./sum(countsDistMoving{i});
    degsFlyDistMoving{i} = linspace(-180,180,length(countsDistMoving{i}));
    subplot(1,3,1), plot(degsFlyDistMoving{i},probabilitiesDistMoving{i})
    hold on
    ylim([0, 0.25]);
    xlim([-180, 180]); xlabel('Distance to the goal (deg)');
    ylabel('Probability');
    title('First low error block');
end
for i = 1:length(dist2goalMovingAllHE)
    [countsDistMovingHE{i}] = histcounts(dist2goalMovingAllHE{i},edges);
    probabilitiesDistMovingHE{i} = countsDistMovingHE{i}./sum(countsDistMovingHE{i});
    degsFlyDistMovingHE{i} = linspace(-180,180,length(countsDistMovingHE{i}));
    subplot(1,3,2), plot(degsFlyDistMovingHE{i},probabilitiesDistMovingHE{i})
    hold on
    ylim([0, 0.25]);
    xlim([-180, 180]); xlabel('Distance to the goal (deg)');
    ylabel('Probability');
    title('High error block');
end
for i = 1:length(dist2goalMovingAll)
    [countsDistMovingSecondLE{i}] = histcounts(dist2goalMovingAllSecondLE{i},edges);
    probabilitiesDistMovingSecondLE{i} = countsDistMovingSecondLE{i}./sum(countsDistMovingSecondLE{i});
    degsFlyDistMovingSecondLE{i} = linspace(-180,180,length(countsDistMovingSecondLE{i}));
    subplot(1,3,3), plot(degsFlyDistMovingSecondLE{i},probabilitiesDistMovingSecondLE{i})
    hold on
    ylim([0, 0.25]);
    xlim([-180, 180]); xlabel('Distance to the goal (deg)');
    ylabel('Probability');
    title('Second low error block');
end

%Plot mean and error for each one
DistanceMoving = zeros(0,0);
allFliesprobMoving = zeros(0,0);
DistanceMovingHE = zeros(0,0);
allFliesprobMovingHE = zeros(0,0);
DistanceMovingSecondLE = zeros(0,0);
allFliesprobMovingSecondLE = zeros(0,0);

for i = 1:length(dist2goalMovingAll)
    DistanceMoving = [DistanceMoving,dist2goalMovingAll{i}];
    allFliesprobMoving = [allFliesprobMoving;probabilitiesDistMoving{1,i}];
    DistanceMovingHE = [DistanceMovingHE,dist2goalMovingAllHE{i}];
    allFliesprobMovingHE = [allFliesprobMovingHE;probabilitiesDistMovingHE{1,i}];
    DistanceMovingSecondLE = [DistanceMovingSecondLE,dist2goalMovingAllSecondLE{i}];
    allFliesprobMovingSecondLE = [allFliesprobMovingSecondLE;probabilitiesDistMovingSecondLE{1,i}];
end

[allcountsMoving] = histcounts(DistanceMoving,edges);
AllprobabilitiesMoving = allcountsMoving./sum(allcountsMoving);
alldegsMoving = linspace(-180,180,length(allcountsMoving));
errorMoving = std(allFliesprobMoving);

[allcountsMovingHE] = histcounts(DistanceMovingHE,edges);
AllprobabilitiesMovingHE = allcountsMovingHE./sum(allcountsMovingHE);
alldegsMovingHE = linspace(-180,180,length(allcountsMovingHE));
errorMovingHE = std(allFliesprobMovingHE);

[allcountsMovingSecondLE] = histcounts(DistanceMovingSecondLE,edges);
AllprobabilitiesMovingSecondLE = allcountsMovingSecondLE./sum(allcountsMovingSecondLE);
alldegsMovingSecondLE = linspace(-180,180,length(allcountsMovingSecondLE));
errorMovingSecondLE = std(allFliesprobMovingSecondLE);

figure,
h1 = boundedline(alldegsMoving,AllprobabilitiesMoving,errorMoving,'k','alpha')
set(h1, 'linewidth', 3);
hold on
h2 = boundedline(alldegsMovingHE,AllprobabilitiesMovingHE,errorMovingHE,'r','alpha')
set(h2, 'linewidth', 3);
h3 = boundedline(alldegsMovingSecondLE,AllprobabilitiesMovingSecondLE,errorMovingSecondLE,'b','alpha')
set(h3, 'linewidth', 3);
%legend('First low error block','High error block','Second low error block');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\Dist2goal.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\Dist2goal.svg'))


%% Distance to the goal in the 10 sec around the jumps

FilesLE = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\shortData10secLowErrorBlockExp*.mat');
FilesHE = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\shortData10secHighErrorBlockExp4.mat');


%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(FilesLE)
    shortData{i} = load(strcat(FilesLE(i).folder,'\',FilesLE(i).name));
    end
for i = 1:length(FilesHE)
    shortDataHE{i} = load(strcat(FilesHE(i).folder,'\',FilesHE(i).name));
end

proba10secAll = zeros(0,0);
proba10secAllHE = zeros(0,0);
for i = 1:length(FilesLE)
    proba10secAll = [proba10secAll,shortData{1,i}.probaDist10secMoving];
end
for i = 1:length(FilesHE)
    proba10secAllHE = [proba10secAllHE,shortDataHE{1,i}.probaDist10secMoving];
end

figure('Position', [300 100 1000 900]),
subplot(1,2,1)
newMap = flipud(gray);
xaxis = [-180:360/17:180];
trials = [1:12*length(shortData)];
imagesc(xaxis,trials,proba10secAll')
colormap(newMap)
colorbar
title('Low error trials')
xlabel('Distance to the goal (deg)'); ylabel('Trial #');
subplot(1,2,2)
imagesc(xaxis,trials,proba10secAllHE')
colormap(newMap)
colorbar
title('High error trials')
xlabel('Distance to the goal (deg)'); ylabel('Trial #');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\Dist2Goal10secHeatMap.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\Dist2Goal10secHeatMap.svg'))

