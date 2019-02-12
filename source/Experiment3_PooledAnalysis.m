%Experiment 3, pooled analysis

%This code pools the data from all the flies in experiment 3
%and does some pooled analyses.

clear all; close all;

%% Around jump velocities

files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\*\*\perTrialData.mat');

%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(files)
data{i} = load(strcat(files(i).folder,'\',files(i).name));
end

%concatenate as columns the arrays for angular velocity, forward
%velocity and jump magnitude for every fly

angVelAll = zeros(0,0);
forwardVelAll = zeros(0,0);
jumpMagAll = zeros(0,0);

for i = 1:length(files)
    angVelAll = [angVelAll,data{1,i}.perTrialData.angVel];
    forwardVelAll = [forwardVelAll,data{1,i}.perTrialData.forwardVel];
    jumpMagAll = [jumpMagAll,data{1,i}.perTrialData.jumpMag];
end


%subset into the 4 magnitude groups
Data45.forwardVel = forwardVelAll(:,jumpMagAll == 45);
Data45.angVel = angVelAll(:,jumpMagAll == 45);
DataNeg45.forwardVel = forwardVelAll(:,jumpMagAll == -45);
DataNeg45.angVel = angVelAll(:,jumpMagAll == -45);
Data90.forwardVel = forwardVelAll(:,jumpMagAll == 90);
Data90.angVel = angVelAll(:,jumpMagAll == 90);
DataNeg90.forwardVel = forwardVelAll(:,jumpMagAll == -90);
DataNeg90.angVel = angVelAll(:,jumpMagAll == -90);


%calculate mean and std ang and forw vel
meanForwardVel45 = mean(Data45.forwardVel,2);
stdForwardVel45 = std(Data45.forwardVel,[],2);
meanForwardVelNeg45 = mean(DataNeg45.forwardVel,2);
stdForwardVelNeg45 = std(DataNeg45.forwardVel,[],2);
meanForwardVel90 = mean(Data90.forwardVel,2);
stdForwardVel90 = std(Data90.forwardVel,[],2);
meanForwardVelNeg90 = mean(DataNeg90.forwardVel,2);
stdForwardVelNeg90 = std(DataNeg90.forwardVel,[],2);

meanAngVel45 = mean(Data45.angVel,2);
stdAngVel45 = std(Data45.angVel,[],2);
meanAngVelNeg45 = mean(DataNeg45.angVel,2);
stdAngVelNeg45 = std(DataNeg45.angVel,[],2);
meanAngVel90 = mean(Data90.angVel,2);
stdAngVel90 = std(Data90.angVel,[],2);
meanAngVelNeg90 = mean(DataNeg90.angVel,2);
stdAngVelNeg90 = std(DataNeg90.angVel,[],2);


%plot mean and error for each
time = linspace(-10,10,length(forwardVelAll));

% figure,
% subplot(1,2,1)
% plot(time,meanForwardVel45)
% hold on
% plot(time,meanForwardVelNeg45)
% plot(time,meanForwardVel90)
% plot(time,meanForwardVelNeg90)
% title('Mean forward velocity');
% legend({'45','-45','90','-90'});
% ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
% subplot(1,2,2)
% plot(time,meanAngVel45)
% hold on
% plot(time,meanAngVelNeg45)
% plot(time,meanAngVel90)
% plot(time,meanAngVelNeg90)
% title('Mean angular velocity');
% legend({'45','-45','90','-90'});
% ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');


figure,
subplot(1,2,1)
p1 = boundedline(time,meanForwardVel45,stdForwardVel45/sqrt(i),'r','alpha')
hold on
p2 = boundedline(time,meanForwardVelNeg45,stdForwardVelNeg45/sqrt(i),'b','alpha')
p3 = boundedline(time,meanForwardVel90,stdForwardVel90/sqrt(i),'y','alpha')
p4 = boundedline(time,meanForwardVelNeg90,stdForwardVelNeg90/sqrt(i),'k','alpha')
title('Mean forward velocity');
legend([p1,p2,p3,p4], '45','-45','90','-90');
ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
ylim([-5,15]);
plot([0,0],[-5,15],'--k','HandleVisibility','off');

subplot(1,2,2)
p12 = boundedline(time,meanAngVel45,stdAngVel45/sqrt(i),'r','alpha')
hold on
p22 = boundedline(time,meanAngVelNeg45,stdAngVelNeg45/sqrt(i),'b','alpha')
p32 = boundedline(time,meanAngVel90,stdAngVel90/sqrt(i),'y','alpha')
p42 = boundedline(time,meanAngVelNeg90,stdAngVelNeg90/sqrt(i),'k','alpha')
title('Mean angular velocity');
ylim([-100;100]);
plot([0,0],[-100,100],'--k','HandleVisibility','off');
legend([p12,p22,p32,p42], '45','-45','90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');
%% distance to the goal

Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\*\*\goals.mat');

%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(Files)
data{i} = load(strcat(Files(i).folder,'\',Files(i).name));
end

%concatenate as columns the arrays for all the variables

goalAll = zeros(0,0);
goalMovingAll = zeros(0,0);
dist2goalAll = {};
dist2goalMovingAll = {};

for i = 1:length(Files)
    goalAll = [goalAll;data{1,i}.goal];
    goalMovingAll = [goalMovingAll;data{1,i}.goalMoving];
    dist2goalAll{i} = data{1,i}.dist2goal;
    dist2goalMovingAll{i} = data{1,i}.dist2goalMoving;
end

%Look at goal distribution
figure,
subplot(1,2,1)
plot(rad2deg(goalAll),'ko')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('Heading goal with every frame');
subplot(1,2,2)
plot(rad2deg(goalMovingAll),'ro')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('Heading goal with only moving frames');

%Looks at the distribution of the variable 'Distance to the goal'
figure,
edges = [-180:20:180];
for i = 1:length(dist2goalAll)
    [countsDist{i}] = histcounts(dist2goalAll{i},edges);
    probabilitiesDist{i} = countsDist{i}./sum(countsDist{i});
    degsFlyDist{i} = linspace(-180,180,length(countsDist{i}));
    subplot(1,2,1), plot(degsFlyDist{i},probabilitiesDist{i})
    hold on
end
title('Distance to the goal with every frame');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

%add mean and error
Distance = zeros(0,0);
allFliesprob = zeros(0,0);
for i = 1:length(dist2goalAll)
    Distance = [Distance,dist2goalAll{i}];
    allFliesprob = [allFliesprob;probabilitiesDist{1,i}];
end
[allcounts] = histcounts(Distance,edges);
Allprobabilities = allcounts./sum(allcounts);
alldegs = linspace(-180,180,length(allcounts));
error = std(allFliesprob);
boundedline(alldegs,Allprobabilities,error,'k','alpha')

%For moving frames only
for i = 1:length(dist2goalMovingAll)
    [countsDistMoving{i}] = histcounts(dist2goalMovingAll{i},edges);
    probabilitiesDistMoving{i} = countsDistMoving{i}./sum(countsDistMoving{i});
    degsFlyDistMoving{i} = linspace(-180,180,length(countsDistMoving{i}));
    subplot(1,2,2), plot(degsFlyDistMoving{i},probabilitiesDistMoving{i})
    hold on
end
title('Distance to the goal with moving frames');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

%add mean and error
DistanceMoving = zeros(0,0);
allFliesprobMoving = zeros(0,0);
for i = 1:length(dist2goalMovingAll)
    DistanceMoving = [DistanceMoving,dist2goalMovingAll{i}];
    allFliesprobMoving = [allFliesprobMoving;probabilitiesDistMoving{1,i}];
end
[allcountsMoving] = histcounts(DistanceMoving,edges);
AllprobabilitiesMoving = allcountsMoving./sum(allcountsMoving);
alldegsMoving = linspace(-180,180,length(allcountsMoving));
errorMoving = std(allFliesprobMoving);
boundedline(alldegsMoving,AllprobabilitiesMoving,errorMoving,'k','alpha')
