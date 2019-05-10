%Experiment 4, pooled analysis

%This code pools the data from all the flies in experiment 4
%and does some pooled analyses.

clear all; close all;

%% Around jump velocities

%files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\*\experimental flies\*\perTrialData.mat');

files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\*\experimental flies\*\Smoothed50perTrialData.mat');

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
angPosAll = zeros(0,0);

for i = 1:length(files)
    angVelAll = [angVelAll,data{1,i}.perTrialData.angVel];
    forwardVelAll = [forwardVelAll,data{1,i}.perTrialData.forwardVel];
    jumpMagAll = [jumpMagAll,data{1,i}.perTrialData.jumpMag];
    angPosAll = [angPosAll,data{1,i}.perTrialData.angPos];
end


%subset into the 2 magnitude groups
Data90.forwardVel = forwardVelAll(:,jumpMagAll == 90);
Data90.angVel = angVelAll(:,jumpMagAll == 90);
DataNeg90.forwardVel = forwardVelAll(:,jumpMagAll == -90);
DataNeg90.angVel = angVelAll(:,jumpMagAll == -90);


%calculate mean and std ang and forw vel
meanForwardVel90 = mean(Data90.forwardVel,2);
stdForwardVel90 = std(Data90.forwardVel,[],2);
meanForwardVelNeg90 = mean(DataNeg90.forwardVel,2);
stdForwardVelNeg90 = std(DataNeg90.forwardVel,[],2);

meanAngVel90 = mean(Data90.angVel,2);
stdAngVel90 = std(Data90.angVel,[],2);
meanAngVelNeg90 = mean(DataNeg90.angVel,2);
stdAngVelNeg90 = std(DataNeg90.angVel,[],2);


%plot mean and error for each
time = linspace(-15,15,length(forwardVelAll));

figure('Position', [100 100 1600 900]),
subplot(1,2,1)
p1 = boundedline(time,meanForwardVel90,stdForwardVel90/sqrt(i),'r','alpha')
hold on
p2 = boundedline(time,meanForwardVelNeg90,stdForwardVelNeg90/sqrt(i),'k','alpha')
title('Mean forward velocity');
legend([p1,p2], '90','-90');
ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
ylim([0, 18]);xlim([-15, 15]);
plot([0,0],[0, 18],'--k','HandleVisibility','off');

subplot(1,2,2)
p32 = boundedline(time,meanAngVel90,stdAngVel90/sqrt(i),'r','alpha')
hold on
p42 = boundedline(time,meanAngVelNeg90,stdAngVelNeg90/sqrt(i),'k','alpha')
title('Mean angular velocity');
ylim([-40, 40]);xlim([-8, 8]);
plot([0,0],[-40, 40],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([,p32,p42], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\AverageVelocityChange.png'))
%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\AverageVelocityChange.svg'))

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Smoothed50AverageVelocityChange.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Smoothed50AverageVelocityChange.svg'))

%% Add Jonathan's quantitative analysis:
%he takes the mean velocity from -2 to 0 and compares it with the one from
%0 to 2, for each fly, using a wilcoxon signed-rank test

jumpFrame = round(size(angVelAll,1)/2+1);

%2 sec = 50 frames, 5 sec = 150 frames

for i = 1:length(data)
BJ90degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,data{1,i}.perTrialData.jumpMag==90));
BJNeg90degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,data{1,i}.perTrialData.jumpMag==-90));

AJ90degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,data{1,i}.perTrialData.jumpMag==90));
AJNeg90degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,data{1,i}.perTrialData.jumpMag==-90));

BJ90degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,data{1,i}.perTrialData.jumpMag==90));
BJNeg90degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,data{1,i}.perTrialData.jumpMag==-90));

AJ90degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,data{1,i}.perTrialData.jumpMag==90));
AJNeg90degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,data{1,i}.perTrialData.jumpMag==-90));

end

%convert to doubles
BJ90degMeanForwardVel  = cell2mat(BJ90degMeanForwardVel);
BJNeg90degMeanForwardVel  = cell2mat(BJNeg90degMeanForwardVel);

AJ90degMeanForwardVel = cell2mat(AJ90degMeanForwardVel);
AJNeg90degMeanForwardVel = cell2mat(AJNeg90degMeanForwardVel);

BJ90degMeanAngVel = cell2mat(BJ90degMeanAngVel);
BJNeg90degMeanAngVel = cell2mat(BJNeg90degMeanAngVel);

AJ90degMeanAngVel = cell2mat(AJ90degMeanAngVel);
AJNeg90degMeanAngVel = cell2mat(AJNeg90degMeanAngVel);


%group BJ and AJ for every angle
around90 = [BJ90degMeanForwardVel;AJ90degMeanForwardVel];
aroundNeg90 = [BJNeg90degMeanForwardVel;AJNeg90degMeanForwardVel];

around90Ang = [BJ90degMeanAngVel;AJ90degMeanAngVel];
aroundNeg90Ang = [BJNeg90degMeanAngVel;AJNeg90degMeanAngVel];


%compute wilcoxon tests
p3 = ranksum(BJ90degMeanForwardVel,AJ90degMeanForwardVel);
p4 = ranksum(BJNeg90degMeanForwardVel,AJNeg90degMeanForwardVel);

p7 = ranksum(BJ90degMeanAngVel,AJ90degMeanAngVel);
p8 = ranksum(BJNeg90degMeanAngVel,AJNeg90degMeanAngVel);


% Plotting them
x = ones(1, length(around90));
figure('Position', [100 100 1600 900]),
%Forward velocity
subplot(2,2,1)
plot(x, around90(1,:),'o')
hold on
plot(2*x,around90(2,:),'o')
xlim([0,3]);
for i = 1:length(around90)
    plot([1,2],[around90(1,i) around90(2,i)],'k')
end
text(1,max(max(around90)),strcat('p = ',num2str(p3)))
title('90 deg jumps')
ylabel('Forward velocity (mm/s)');
set(gca,'xticklabel',{[]})
subplot(2,2,2)
plot(x, aroundNeg90(1,:),'o')
hold on
plot(2*x,aroundNeg90(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg90)
    plot([1,2],[aroundNeg90(1,i) aroundNeg90(2,i)],'k')
end
text(1,max(max(aroundNeg90)),strcat('p = ',num2str(p4)))
title('-90 deg jumps');
set(gca,'xticklabel',{[]})
%Angular velocity
subplot(2,2,3)
plot(x, around90Ang(1,:),'o')
hold on
plot(2*x,around90Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(around90Ang)
    plot([1,2],[around90Ang(1,i) around90Ang(2,i)],'k')
end
ylabel('Angular velocity (deg/s)');
text(1,max(max(around90Ang)),strcat('p = ',num2str(p7)))
set(gca,'xticklabel',{[]})
subplot(2,2,4)
plot(x, aroundNeg90Ang(1,:),'o')
hold on
plot(2*x,aroundNeg90Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg90Ang)
    plot([1,2],[aroundNeg90Ang(1,i) aroundNeg90Ang(2,i)],'k')
end
text(1,max(max(aroundNeg90Ang)),strcat('p = ',num2str(p8)))
set(gca,'xticklabel',{[]})
%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\QuantVelChange.png'))

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Smoothed50QuantVelChange.png'))

%% 

%Using mixed models statistics
%random factor: fly, fixed factors: time (BJ/AJ), jump magnitude

FlyID = zeros(0,0);
BJVelocity = zeros(0,0);
AJVelocity = zeros(0,0);
JumpMagnitude = zeros(0,0);

for i = 1:length(data)
    FlyID = [FlyID, repelem(i,size(data{1,i}.perTrialData.forwardVel,2))]; %we repeat the fly's idea for double the number of trials, since we will compute the velocity on a given trial before and after the jump
    BJVelocity = [BJVelocity,mean(data{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,:))];
    AJVelocity = [AJVelocity,mean(data{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,:))];
    JumpMagnitude = [JumpMagnitude,data{1,i}.perTrialData.jumpMag];
end

FlyID = [FlyID,FlyID];
FlyID = sprintfc('%d',FlyID);
Velocity = [BJVelocity,AJVelocity];
JumpMagnitude = [JumpMagnitude,JumpMagnitude];
Time = {};
Time(1:length(BJVelocity)) = {'BJ'};
Time(length(BJVelocity)+1:length(BJVelocity)+length(AJVelocity)) = {'AJ'};

factorsVel = table(Velocity',FlyID',Time',JumpMagnitude');
factorsVel.Properties.VariableNames = {'Velocity','FlyID','Time','JumpMagnitude'};

%run a mixed-effects model
glme0 = fitglme(factorsVel,'Velocity ~ 1 + Time*JumpMagnitude +(1|FlyID)');
disp(glme0)
stats = anova(glme0)

%% distance to the goal

%Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\*\experimental flies\*\goals.mat');
Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\*\experimental flies\*\Smoothedgoals.mat');

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
figure('Position', [100 100 1600 900]),
subplot(2,2,1)
plot(rad2deg(goalAll),'ko')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('Heading goal with every frame');
subplot(2,2,2)
plot(rad2deg(goalMovingAll),'ro')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('Heading goal with only moving frames');
subplot(2,2,3)
polarhistogram(goalAll,10,'FaceColor','black')
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top';
subplot(2,2,4)
polarhistogram(goalMovingAll,10,'FaceColor','red')
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top';

%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\GoalDistribution.png'))
%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\GoalDistribution.svg'))

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\SmoothedGoalDistribution.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\SmoothedGoalDistribution.svg'))


%Looks at the distribution of the variable 'Distance to the goal'
figure('Position', [100 100 1600 900]),
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
h1 = boundedline(alldegs,Allprobabilities,error,'k','alpha')
set(h1, 'linewidth', 3);
ylim([0 0.27]);

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
h2 = boundedline(alldegsMoving,AllprobabilitiesMoving,errorMoving,'k','alpha')
set(h2, 'linewidth', 3);
ylim([0 0.27]);
%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Dist2goal.png'))
%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Dist2goal.svg'))

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\SmoothedDist2goal.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\SmoothedDist2goal.svg'))


%% Regress distance to goal using several variables

% I need to get the FlyData from each folder
Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\*\experimental flies\*\flyData.mat');

%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(Files)
flyData{i} = load(strcat(Files(i).folder,'\',Files(i).name));
TimeStamp{i} = cellstr(Files(i).date);
Hours(i) = hour(TimeStamp{i});
Day(i) = day(TimeStamp{i});
Genotype{i} = flyData{1,i}.FlyData.line;
end

%I can take the half-width and my y variable and try to account for that
%with my predictors

% 1) Get max height of every probability distribution

for i = 1:length(Files)
height(i) = max(probabilitiesDistMoving{i});
end

%assess circadian time
circPeakCS = 11;
circPeakDL = 16;

circTime = zeros(1,length(Files));

for i = 1:length(Files)
    if Genotype{i} == 'CS'
        circTime(i) = Hours(i)-circPeakCS;
    else
        circTime(i) = Hours(i)-circPeakDL;
    end    
    
    if circTime(i) < 0
        circTime(i) = 12 + circTime(i);
    end
end


%create a table with my response variable and my predictors
factors = table(height',Hours',Day',Genotype',circTime');
factors.Properties.VariableNames = {'Height' 'Hours' 'Day' 'Genotype' 'circTime'};

% %run a mixed-effects model
% glme = fitglme(factors,'Height ~ 1 + circTime + Hours + Day + Genotype');
% disp(glme)



%boxplots of the height according to these different factors

figure('Position', [100 100 1600 900]), 
subplot(2,2,1), boxplot(height,circTime); title('Contribution of the circadian time');
subplot(2,2,2), boxplot(height,Hours); title('Contribution of the hour of the day');
subplot(2,2,3), boxplot(height,Day); title('Contribution of the day');
subplot(2,2,4), boxplot(height,Genotype); title('Contribution of the genotype');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\FactorEffects.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\FactorEffects.svg'))



%% Distance to the goal in the 100 sec around the jumps

%Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\*\experimental flies\*\shortData.mat');
Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\*\experimental flies\*\SmoothedshortData.mat');


%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(Files)
shortData{i} = load(strcat(Files(i).folder,'\',Files(i).name));
end

proba100secAll = zeros(0,0);

for i = 1:length(files)
    proba100secAll = [proba100secAll,shortData{1,i}.probaDist100secMoving];
end

newMap = flipud(gray);
xaxis = [-180:360/17:180];
trials = [1:110];
figure('Position', [300 100 1000 900]), imagesc(xaxis,trials,proba100secAll')
colormap(newMap)
colorbar
title('Distribution of the distance to the goal 100 sec before to 100 sec after jumps')
xlabel('Distance to the goal (deg)'); ylabel('Trial #');

%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Dist2Goal100sec.png'))
%saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Dist2Goal100sec.svg'))

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\SmoothedDist2Goal100sec.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\SmoothedDist2Goal100sec.svg'))

%% Distance to the goal distribution for empty trials


Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\*\experimental flies\*\EmptyTrialGoals.mat');

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
figure('Position', [100 100 1600 900]),
subplot(2,2,1)
plot(rad2deg(goalAll),'ko')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('Heading goal with every frame');
subplot(2,2,2)
plot(rad2deg(goalMovingAll),'ro')
ylim([-180,180]);
xlabel('Fly #');ylabel('Heading goadl (deg)');
title('Heading goal with only moving frames');
subplot(2,2,3)
polarhistogram(goalAll,10,'FaceColor','black')
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top';
subplot(2,2,4)
polarhistogram(goalMovingAll,10,'FaceColor','red')
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top';

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\GoalDistributionEmptyTrial.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\GoalDistributionEmptyTrial.svg'))


%Looks at the distribution of the variable 'Distance to the goal'
figure('Position', [100 100 1600 900]),
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
h1 = boundedline(alldegs,Allprobabilities,error,'k','alpha')
set(h1, 'linewidth', 3); ylim([0,0.25]); 

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
h2 = boundedline(alldegsMoving,AllprobabilitiesMoving,errorMoving,'k','alpha')
set(h2, 'linewidth', 3);
ylim([0,0.25]);

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Dist2goalEmmptyTrial.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment4\Dist2goalEmptyTrial.svg'))


