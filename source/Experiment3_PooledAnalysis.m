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
angPosAll = zeros(0,0);

for i = 1:length(files)
    angVelAll = [angVelAll,data{1,i}.perTrialData.angVel];
    forwardVelAll = [forwardVelAll,data{1,i}.perTrialData.forwardVel];
    jumpMagAll = [jumpMagAll,data{1,i}.perTrialData.jumpMag];
    angPosAll = [angPosAll,data{1,i}.perTrialData.angPos];
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
time = linspace(-15,15,length(forwardVelAll));

figure('Position', [100 100 1600 900]),
subplot(1,2,1)
p1 = boundedline(time,meanForwardVel45,stdForwardVel45/sqrt(i),'r','alpha')
hold on
p2 = boundedline(time,meanForwardVelNeg45,stdForwardVelNeg45/sqrt(i),'b','alpha')
p3 = boundedline(time,meanForwardVel90,stdForwardVel90/sqrt(i),'y','alpha')
p4 = boundedline(time,meanForwardVelNeg90,stdForwardVelNeg90/sqrt(i),'k','alpha')
title('Mean forward velocity');
legend([p1,p2,p3,p4], '45','-45','90','-90');
ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
ylim([0, 8]);xlim([-15, 15]);
plot([0,0],[0, 8],'--k','HandleVisibility','off');

subplot(1,2,2)
p12 = boundedline(time,meanAngVel45,stdAngVel45/sqrt(i),'r','alpha')
hold on
p22 = boundedline(time,meanAngVelNeg45,stdAngVelNeg45/sqrt(i),'b','alpha')
p32 = boundedline(time,meanAngVel90,stdAngVel90/sqrt(i),'y','alpha')
p42 = boundedline(time,meanAngVelNeg90,stdAngVelNeg90/sqrt(i),'k','alpha')
title('Mean angular velocity');
ylim([-80, 80]);xlim([-8, 8]);
plot([0,0],[-80, 80],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p12,p22,p32,p42], '45','-45','90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\AverageVelocityChange.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\AverageVelocityChange.svg'))

%% Add Jonathan's quantitative analysis:
%he takes the mean velocity from -2 to 0 and compares it with the one from
%0 to 2, for each fly, using a wilcoxon signed-rank test

jumpFrame = round(size(angVelAll,1)/2+1);

%2 sec = 50 frames, 5 sec = 150 frames

for i = 1:length(data)
BJ45MeandegForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,data{1,i}.perTrialData.jumpMag==45));
BJNeg45MeandegForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,data{1,i}.perTrialData.jumpMag==-45));
BJ90degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,data{1,i}.perTrialData.jumpMag==90));
BJNeg90degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,data{1,i}.perTrialData.jumpMag==-90));

AJ45degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,data{1,i}.perTrialData.jumpMag==45));
AJNeg45degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,data{1,i}.perTrialData.jumpMag==-45));
AJ90degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,data{1,i}.perTrialData.jumpMag==90));
AJNeg90degMeanForwardVel{i}= mean(data{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,data{1,i}.perTrialData.jumpMag==-90));


BJ45MeandegAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,data{1,i}.perTrialData.jumpMag==45));
BJNeg45MeandegAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,data{1,i}.perTrialData.jumpMag==-45));
BJ90degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,data{1,i}.perTrialData.jumpMag==90));
BJNeg90degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,data{1,i}.perTrialData.jumpMag==-90));

AJ45degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,data{1,i}.perTrialData.jumpMag==45));
AJNeg45degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,data{1,i}.perTrialData.jumpMag==-45));
AJ90degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,data{1,i}.perTrialData.jumpMag==90));
AJNeg90degMeanAngVel{i}= mean(data{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,data{1,i}.perTrialData.jumpMag==-90));

end

%convert to doubles
BJ45MeandegForwardVel = cell2mat(BJ45MeandegForwardVel);
BJNeg45MeandegForwardVel  = cell2mat(BJNeg45MeandegForwardVel);
BJ90degMeanForwardVel  = cell2mat(BJ90degMeanForwardVel);
BJNeg90degMeanForwardVel  = cell2mat(BJNeg90degMeanForwardVel);

AJ45degMeanForwardVel = cell2mat(AJ45degMeanForwardVel);
AJNeg45degMeanForwardVel = cell2mat(AJNeg45degMeanForwardVel);
AJ90degMeanForwardVel = cell2mat(AJ90degMeanForwardVel);
AJNeg90degMeanForwardVel = cell2mat(AJNeg90degMeanForwardVel);

BJ45MeandegAngVel = cell2mat(BJ45MeandegAngVel);
BJNeg45MeandegAngVel = cell2mat(BJNeg45MeandegAngVel);
BJ90degMeanAngVel = cell2mat(BJ90degMeanAngVel);
BJNeg90degMeanAngVel = cell2mat(BJNeg90degMeanAngVel);

AJ45degMeanAngVel = cell2mat(AJ45degMeanAngVel);
AJNeg45degMeanAngVel = cell2mat(AJNeg45degMeanAngVel);
AJ90degMeanAngVel = cell2mat(AJ90degMeanAngVel);
AJNeg90degMeanAngVel = cell2mat(AJNeg90degMeanAngVel);


%group BJ and AJ for every angle
around45 = [BJ45MeandegForwardVel;AJ45degMeanForwardVel];
aroundNeg45 = [BJNeg45MeandegForwardVel;AJNeg45degMeanForwardVel];
around90 = [BJ90degMeanForwardVel;AJ90degMeanForwardVel];
aroundNeg90 = [BJNeg90degMeanForwardVel;AJNeg90degMeanForwardVel];

around45Ang = [BJ45MeandegAngVel;AJ45degMeanAngVel];
aroundNeg45Ang = [BJNeg45MeandegAngVel;AJNeg45degMeanAngVel];
around90Ang = [BJ90degMeanAngVel;AJ90degMeanAngVel];
aroundNeg90Ang = [BJNeg90degMeanAngVel;AJNeg90degMeanAngVel];


%compute wilcoxon tests
p1 = ranksum(BJ45MeandegForwardVel,AJ45degMeanForwardVel);
p2 = ranksum(BJNeg45MeandegForwardVel,AJNeg45degMeanForwardVel);
p3 = ranksum(BJ90degMeanForwardVel,AJ90degMeanForwardVel);
p4 = ranksum(BJNeg90degMeanForwardVel,AJNeg90degMeanForwardVel);

p5 = ranksum(BJ45MeandegAngVel,AJ45degMeanAngVel);
p6 = ranksum(BJNeg45MeandegAngVel,AJNeg45degMeanAngVel);
p7 = ranksum(BJ90degMeanAngVel,AJ90degMeanAngVel);
p8 = ranksum(BJNeg90degMeanAngVel,AJNeg90degMeanAngVel);


% Plotting them
x = ones(1, length(around45));
figure('Position', [100 100 1600 900]),
%Forward velocity
subplot(2,4,1)
plot(x, around45(1,:),'o') %before jump
hold on
plot(2*x,around45(2,:),'o') %after jump
xlim([0,3]);
for i = 1:length(around45)
    plot([1,2],[around45(1,i) around45(2,i)],'k') %add lines to join them
end
ylabel('Forward velocity (mm/s)');
title('45 deg jumps')
text(1,max(max(around45)),strcat('p = ',num2str(p1))) %add p value
set(gca,'xticklabel',{[]})
subplot(2,4,2)
plot(x, aroundNeg45(1,:),'o')
hold on
plot(2*x,aroundNeg45(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg45)
    plot([1,2],[aroundNeg45(1,i) aroundNeg45(2,i)],'k')
end
text(1,max(max(aroundNeg45)),strcat('p = ',num2str(p2)))
title('-45 deg jumps')
set(gca,'xticklabel',{[]})
subplot(2,4,3)
plot(x, around90(1,:),'o')
hold on
plot(2*x,around90(2,:),'o')
xlim([0,3]);
for i = 1:length(around90)
    plot([1,2],[around90(1,i) around90(2,i)],'k')
end
text(1,max(max(around90)),strcat('p = ',num2str(p3)))
title('90 deg jumps')
set(gca,'xticklabel',{[]})
subplot(2,4,4)
plot(x(2:end), aroundNeg90(1,:),'o')
hold on
plot(2*x(2:end),aroundNeg90(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg90)
    plot([1,2],[aroundNeg90(1,i) aroundNeg90(2,i)],'k')
end
text(1,max(max(aroundNeg90)),strcat('p = ',num2str(p4)))
title('-90 deg jumps');
set(gca,'xticklabel',{[]})
%Angular velocity
subplot(2,4,5)
plot(x, around45Ang(1,:),'o')
hold on
plot(2*x,around45Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(around45Ang)
    plot([1,2],[around45Ang(1,i) around45Ang(2,i)],'k')
end
ylabel('Angular velocity (deg/s)');
text(1,max(max(around45Ang)),strcat('p = ',num2str(p5)))
set(gca,'xticklabel',{[]})
subplot(2,4,6)
plot(x, aroundNeg45Ang(1,:),'o')
hold on
plot(2*x,aroundNeg45Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg45Ang)
    plot([1,2],[aroundNeg45Ang(1,i) aroundNeg45Ang(2,i)],'k')
end
text(1,max(max(aroundNeg45Ang)),strcat('p = ',num2str(p6)))
set(gca,'xticklabel',{[]})
subplot(2,4,7)
plot(x, around90Ang(1,:),'o')
hold on
plot(2*x,around90Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(around90Ang)
    plot([1,2],[around90Ang(1,i) around90Ang(2,i)],'k')
end
text(1,max(max(around90Ang)),strcat('p = ',num2str(p7)))
set(gca,'xticklabel',{[]})
subplot(2,4,8)
plot(x(2:end), aroundNeg90Ang(1,:),'o')
hold on
plot(2*x(2:end),aroundNeg90Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg90Ang)
    plot([1,2],[aroundNeg90Ang(1,i) aroundNeg90Ang(2,i)],'k')
end
text(1,max(max(aroundNeg90Ang)),strcat('p = ',num2str(p8)))
set(gca,'xticklabel',{[]})

%I see an effect for the angular velocity in 90 and -90 if I compare every
%individual trial BJ and AJ, but not if I take the mean for each fly...


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

Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\*\*\goals.mat');

%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(Files)
data{i} = load(strcat(Files(i).folder,'\',Files(i).name));
end

%concatenate as columns the arrays for all the variables

goalAll = zeros(0,0);
goalMovingAll = zeros(0,0);
goalMovingAll300sec = zeros(0,0);
dist2goalAll = {};
dist2goalMovingAll = {};
dist2goalMovingAll300sec = {};

for i = 1:length(Files)
    goalAll = [goalAll;data{1,i}.goal];
    goalMovingAll = [goalMovingAll;data{1,i}.goalMoving];
    goalMovingAll300sec = [goalMovingAll300sec;data{1,i}.goalMoving300sec];
    dist2goalAll{i} = data{1,i}.dist2goal;
    dist2goalMovingAll{i} = data{1,i}.dist2goalMoving;
    dist2goalMovingAll300sec{i} = data{1,i}.dist2goalMoving300sec;
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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\GoalDistribution.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\GoalDistribution.svg'))


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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\Dist2goal.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\Dist2goal.svg'))



%Look at dist2goal for the first 300 sec
Distance300sec = zeros(0,0);
allFliesprob300sec = zeros(0,0);
for i = 1:length(dist2goalMovingAll300sec)
    [countsDistMoving300sec{i}] = histcounts(dist2goalMovingAll300sec{i},edges);
    probabilitiesDistMoving300sec{i} = countsDistMoving300sec{i}./sum(countsDistMoving300sec{i});
    Distance300sec = [Distance300sec,dist2goalMovingAll300sec{i}];
    allFliesprob300sec = [allFliesprob300sec;probabilitiesDistMoving300sec{1,i}];
end
[allcounts300sec] = histcounts(Distance300sec,edges);
Allprobabilities300sec = allcounts300sec./sum(allcounts300sec);
alldegs300sec = linspace(-180,180,length(allcounts300sec));
error300sec = std(allFliesprob300sec);
figure,
h3 = boundedline(alldegs300sec,Allprobabilities300sec,error300sec,'k','alpha')
set(h3, 'linewidth', 3);
title('Distance to the goal with moving frames');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

%% Regress distance to goal using several variables

% I need to get the FlyData from each folder
Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\*\*\flyData.mat');

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
height300sec(i) = max(probabilitiesDistMoving300sec{i});
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
factors = table(height',Hours',Day',Genotype',circTime',height300sec');
factors.Properties.VariableNames = {'Height' 'Hours' 'Day' 'Genotype' 'circTime' 'height300sec'};

%run a mixed-effects model
glme = fitglme(factors,'Height ~ 1 + circTime + Hours + Day + Genotype + height300sec');
disp(glme)

glme2 = fitglme(factors,'Height ~ 1 + circTime + Hours + Day + Genotype');
disp(glme2)
%the day seems to be slightly significant, and maybe the genotype slightly.
% I don't think there are clear effects

%There is a clear relationship with the menotaxis during the first 300 sec

%boxplots of the height according to these different factors

figure('Position', [100 100 1600 900]), 
subplot(2,2,1), boxplot(height,circTime); title('Contribution of the circadian time');
subplot(2,2,2), boxplot(height,Hours); title('Contribution of the hour of the day');
subplot(2,2,3), boxplot(height,Day); title('Contribution of the day');
subplot(2,2,4), boxplot(height,Genotype); title('Contribution of the genotype');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\FactorEffects.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\FactorEffects.svg'))


figure('Position', [100 100 1200 900]),
scatter(height300sec,height,'k')
title('Correlation between the menotactic activity during the full experiment and the first 5 min')
xlabel('Menotaxis during first 5 min');
ylabel('Menotaxis during all experiment');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\menotacticCorr.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\menotacticCorr.svg'))

%% Distance to the goal in the 100 sec around the jumps

Files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\*\*\shortData.mat');

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
figure, imagesc(proba100secAll')
colormap(newMap)
colorbar
title('Distribution of the distance to the goal 100 sec before to 100 sec after jumps')
xlabel('Distance to the goal (deg)'); ylabel('Trial #');



%% Dist to the goal 10 sec before jump in time

%He keeps only those flies that had a circular sd <45 deg before the jump,
%so I should start by assessing that

%Define the goal per jump as the average heading in the 10 s preceding it.
jumpFrame = round(size(angVelAll,1)/2+1);
angPosBJ = angPosAll(jumpFrame-250:jumpFrame,:); %get the position in the time before the jump
deviationHeadingBJ = rad2deg(circ_std(deg2rad(angPosBJ))); %get the std of that position

%exclude trials with a deviation above 45 deg
narrowFlies = find(deviationHeadingBJ<45);
angPosBJselect = angPosBJ(:,narrowFlies);

narrowGoal = circ_mean(deg2rad(angPosBJselect)); %calculate the goal as the mean heading

%Calulate the distance to the goal at every time point with respect to that
%goal
angPosAllselect = angPosAll(:,deviationHeadingBJ<45);
for i = 1:length(narrowGoal)
narrowDist2Goal{i} = circ_dist(deg2rad(angPosAllselect(:,i)),narrowGoal(i));
end

narrowDist2Goal = cell2mat(narrowDist2Goal);

% % I have to didive into the 4 different groups now to see the goal shift to
% % those directions
% narrowDist2Goal45deg = narrowDist2Goal(:,jumpMagAll == 45);
% narrowDist2GoalNeg45deg = narrowDist2Goal(:,jumpMagAll == -45);
% narrowDist2Goal90deg = narrowDist2Goal(:,jumpMagAll == 90);
% narrowDist2GoalNeg90deg = narrowDist2Goal(:,jumpMagAll == -90);
% 
% % %Check that after the jumps they are offset....
% % figure,
% % subplot(2,2,1), polarhistogram(narrowDist2Goal45deg(254,:))
% % subplot(2,2,2), polarhistogram(narrowDist2GoalNeg45deg(254,:))
% % subplot(2,2,3), polarhistogram(narrowDist2Goal90deg(254,:))
% % subplot(2,2,4), polarhistogram(narrowDist2GoalNeg90deg(254,:))
% 
% %Plot like he does it: distributions of distance from goal at 0-1 s (red)
% %and at 10-20 s (black) after the barjump
% 
% %0-1 s = frames 252 to 277
% %10-20 s = frames 501 to 751
% afterjumpFrames = [round(size(angVelAll,1)/2)+1:round(size(angVelAll,1)/2)+26];
% LaterFrames = [jumpFrame+250:size(angVelAll,1)];
% 
% narrowDist2Goal45degAJ = wrapTo180(rad2deg(circ_mean(narrowDist2Goal45deg(afterjumpFrames,:))));
% narrowDist2GoalNeg45degAJ = wrapTo180(rad2deg(circ_mean(narrowDist2GoalNeg45deg(afterjumpFrames,:))));
% narrowDist2Goal90degAJ = wrapTo180(rad2deg(circ_mean(narrowDist2Goal90deg(afterjumpFrames,:))));
% narrowDist2GoalNeg90degAJ = wrapTo180(rad2deg(circ_mean(narrowDist2GoalNeg90deg(afterjumpFrames,:))));
% 
% narrowDist2Goal45degLater = wrapTo180(rad2deg(circ_mean(narrowDist2Goal45deg(LaterFrames,:))));
% narrowDist2GoalNeg45degLater = wrapTo180(rad2deg(circ_mean(narrowDist2GoalNeg45deg(LaterFrames,:))));
% narrowDist2Goal90degLater = wrapTo180(rad2deg(circ_mean(narrowDist2Goal90deg(LaterFrames,:))));
% narrowDist2GoalNeg90degLater = wrapTo180(rad2deg(circ_mean(narrowDist2GoalNeg90deg(LaterFrames,:))));
% 
% 
% figure,
% for i = 1:size(narrowDist2Goal45degAJ,1)
%     [countsDist45deg{i}] = histcounts(narrowDist2Goal45degAJ(i),edges);
%     probabilitiesDist45deg{i} = countsDist45deg{i}./sum(countsDist45deg{i});
%     degsFlyDist45deg{i} = linspace(-180,180,length(countsDist45deg{i}));
%     subplot(2,2,1), plot(degsFlyDist45deg{i},probabilitiesDist45deg{i},'r')    
%     
%     [countsDistNeg45deg{i}] = histcounts(narrowDist2GoalNeg45degAJ(i),edges);
%     probabilitiesDistNeg45deg{i} = countsDistNeg45deg{i}./sum(countsDistNeg45deg{i});
%     degsFlyDistNeg45deg{i} = linspace(-180,180,length(countsDistNeg45deg{i}));
%     subplot(2,2,2), plot(degsFlyDistNeg45deg{i},probabilitiesDistNeg45deg{i},'r')
%     
%     [countsDist90deg{i}] = histcounts(narrowDist2Goal90degAJ(i),edges);
%     probabilitiesDist90deg{i} = countsDist90deg{i}./sum(countsDist90deg{i});
%     degsFlyDist90deg{i} = linspace(-180,180,length(countsDist90deg{i}));
%     subplot(2,2,3), plot(degsFlyDist90deg{i},probabilitiesDist90deg{i},'r')
%     
%     [countsDistNeg90deg{i}] = histcounts(narrowDist2GoalNeg90degAJ(i),edges);
%     probabilitiesDistNeg90deg{i} = countsDistNeg90deg{i}./sum(countsDistNeg90deg{i});
%     degsFlyDistNeg90deg{i} = linspace(-180,180,length(countsDistNeg90deg{i}));
%     subplot(2,2,4), plot(degsFlyDistNeg90deg{i},probabilitiesDistNeg90deg{i},'r')
%         
% end
% 
% 
% figure,
% for i = 1:size(narrowDist2Goal45degLater,1)
%     [countsDist45degLater{i}] = histcounts(narrowDist2Goal45degLater(i),edges);
%     probabilitiesDist45degLater{i} = countsDist45degLater{i}./sum(countsDist45degLater{i});
%     degsFlyDist45degLater{i} = linspace(-180,180,length(countsDist45degLater{i}));
%     subplot(2,2,1), plot(degsFlyDist45degLater{i},probabilitiesDist45degLater{i},'k')    
%     
%     [countsDistNeg45degLater{i}] = histcounts(narrowDist2GoalNeg45degLater(i),edges);
%     probabilitiesDistNeg45degLater{i} = countsDistNeg45degLater{i}./sum(countsDistNeg45degLater{i});
%     degsFlyDistNeg45degLater{i} = linspace(-180,180,length(countsDistNeg45degLater{i}));
%     subplot(2,2,2), plot(degsFlyDistNeg45degLater{i},probabilitiesDistNeg45degLater{i},'k')
%     
%     [countsDist90degLater{i}] = histcounts(narrowDist2Goal90degLater(i),edges);
%     probabilitiesDist90degLater{i} = countsDist90degLater{i}./sum(countsDist90degLater{i});
%     degsFlyDist90degLater{i} = linspace(-180,180,length(countsDist90degLater{i}));
%     subplot(2,2,3), plot(degsFlyDist90degLater{i},probabilitiesDist90degLater{i},'k')
%     
%     [countsDistNeg90degLater{i}] = histcounts(narrowDist2GoalNeg90degLater(i),edges);
%     probabilitiesDistNeg90degLater{i} = countsDistNeg90degLater{i}./sum(countsDistNeg90degLater{i});
%     degsFlyDistNeg90degLater{i} = linspace(-180,180,length(countsDistNeg90degLater{i}));
%     subplot(2,2,4), plot(degsFlyDistNeg90degLater{i},probabilitiesDistNeg90degLater{i},'k')
%         
% end

%Either I did something very wrong, or they're not bringing the bar back
%like they should...

%% Choose the flies that menotaxed well based on the height of their dist2goal,
%and repeat the velocity plots for them

%Keep flies with height of the 'Dist2goal' variable above average

goodFlies = find(height>mean(height));

files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\*\*\perTrialData.mat');

%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(files)
data{i} = load(strcat(files(i).folder,'\',files(i).name));
end

%Keep the data for the 'good flies'

for i = 1:length(goodFlies)
    goodData{i} = data{goodFlies(i)};
end


%concatenate as columns the arrays for angular velocity, forward
%velocity and jump magnitude for every fly

angVelAll = zeros(0,0);
forwardVelAll = zeros(0,0);
jumpMagAll = zeros(0,0);
angPosAll = zeros(0,0);

for i = 1:length(goodFlies)
    angVelAll = [angVelAll,goodData{1,i}.perTrialData.angVel];
    forwardVelAll = [forwardVelAll,goodData{1,i}.perTrialData.forwardVel];
    jumpMagAll = [jumpMagAll,goodData{1,i}.perTrialData.jumpMag];
    angPosAll = [angPosAll,goodData{1,i}.perTrialData.angPos];
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
time = linspace(-15,15,length(forwardVelAll));

figure('Position', [100 100 1600 900]),
subplot(1,2,1)
p1 = boundedline(time,meanForwardVel45,stdForwardVel45/sqrt(i),'r','alpha')
hold on
p2 = boundedline(time,meanForwardVelNeg45,stdForwardVelNeg45/sqrt(i),'b','alpha')
p3 = boundedline(time,meanForwardVel90,stdForwardVel90/sqrt(i),'y','alpha')
p4 = boundedline(time,meanForwardVelNeg90,stdForwardVelNeg90/sqrt(i),'k','alpha')
title('Mean forward velocity for good flies');
legend([p1,p2,p3,p4], '45','-45','90','-90');
ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
ylim([-2, 14]);xlim([-15, 15]);
plot([0,0],[-2, 14],'--k','HandleVisibility','off');

subplot(1,2,2)
p12 = boundedline(time,meanAngVel45,stdAngVel45/sqrt(i),'r','alpha')
hold on
p22 = boundedline(time,meanAngVelNeg45,stdAngVelNeg45/sqrt(i),'b','alpha')
p32 = boundedline(time,meanAngVel90,stdAngVel90/sqrt(i),'y','alpha')
p42 = boundedline(time,meanAngVelNeg90,stdAngVelNeg90/sqrt(i),'k','alpha')
title('Mean angular velocity for good flies');
ylim([-80, 80]);xlim([-8, 8]);
plot([0,0],[-80, 80],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p12,p22,p32,p42], '45','-45','90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

%with no errors, for a clearer view
figure,
subplot(1,2,1)
plot(time,meanForwardVel45,'r')
hold on
plot(time,meanForwardVelNeg45,'m')
plot(time,meanForwardVel90,'b')
plot(time,meanForwardVelNeg90,'c')
title('Mean forward velocity for good flies');
legend('45','-45','90','-90');
ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
ylim([-2, 14]);xlim([-15, 15]);
plot([0,0],[-2, 14],'--k','HandleVisibility','off');

subplot(1,2,2)
plot(time,meanAngVel45,'r')
hold on
plot(time,meanAngVelNeg45,'m')
plot(time,meanAngVel90,'b')
plot(time,meanAngVelNeg90,'c')
title('Mean angular velocity for good flies');
ylim([-80, 80]);xlim([-8, 8]);
plot([0,0],[-80, 80],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p12,p22,p32,p42], '45','-45','90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

%% Calculate the difference in velocities for these good flies

for i = 1:length(goodData)
GBJ45MeandegForwardVel{i}= mean(goodData{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,goodData{1,i}.perTrialData.jumpMag==45));
GBJNeg45MeandegForwardVel{i}= mean(goodData{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,goodData{1,i}.perTrialData.jumpMag==-45));
GBJ90degMeanForwardVel{i}= mean(goodData{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,goodData{1,i}.perTrialData.jumpMag==90));
GBJNeg90degMeanForwardVel{i}= mean(goodData{1,i}.perTrialData.forwardVel(jumpFrame-250:jumpFrame,goodData{1,i}.perTrialData.jumpMag==-90));

GAJ45degMeanForwardVel{i}= mean(goodData{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,goodData{1,i}.perTrialData.jumpMag==45));
GAJNeg45degMeanForwardVel{i}= mean(goodData{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,goodData{1,i}.perTrialData.jumpMag==-45));
GAJ90degMeanForwardVel{i}= mean(goodData{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,goodData{1,i}.perTrialData.jumpMag==90));
GAJNeg90degMeanForwardVel{i}= mean(goodData{1,i}.perTrialData.forwardVel(jumpFrame:jumpFrame+250,goodData{1,i}.perTrialData.jumpMag==-90));


GBJ45MeandegAngVel{i}= mean(goodData{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,goodData{1,i}.perTrialData.jumpMag==45));
GBJNeg45MeandegAngVel{i}= mean(goodData{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,goodData{1,i}.perTrialData.jumpMag==-45));
GBJ90degMeanAngVel{i}= mean(goodData{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,goodData{1,i}.perTrialData.jumpMag==90));
GBJNeg90degMeanAngVel{i}= mean(goodData{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,goodData{1,i}.perTrialData.jumpMag==-90));

GAJ45degMeanAngVel{i}= mean(goodData{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,goodData{1,i}.perTrialData.jumpMag==45));
GAJNeg45degMeanAngVel{i}= mean(goodData{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,goodData{1,i}.perTrialData.jumpMag==-45));
GAJ90degMeanAngVel{i}= mean(goodData{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,goodData{1,i}.perTrialData.jumpMag==90));
GAJNeg90degMeanAngVel{i}= mean(goodData{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,goodData{1,i}.perTrialData.jumpMag==-90));

end

%convert to doubles
GBJ45MeandegForwardVel = cell2mat(GBJ45MeandegForwardVel);
GBJNeg45MeandegForwardVel  = cell2mat(GBJNeg45MeandegForwardVel);
GBJ90degMeanForwardVel  = cell2mat(GBJ90degMeanForwardVel);
GBJNeg90degMeanForwardVel  = cell2mat(GBJNeg90degMeanForwardVel);

GAJ45degMeanForwardVel = cell2mat(GAJ45degMeanForwardVel);
GAJNeg45degMeanForwardVel = cell2mat(GAJNeg45degMeanForwardVel);
GAJ90degMeanForwardVel = cell2mat(GAJ90degMeanForwardVel);
GAJNeg90degMeanForwardVel = cell2mat(GAJNeg90degMeanForwardVel);

GBJ45MeandegAngVel = cell2mat(GBJ45MeandegAngVel);
GBJNeg45MeandegAngVel = cell2mat(GBJNeg45MeandegAngVel);
GBJ90degMeanAngVel = cell2mat(GBJ90degMeanAngVel);
GBJNeg90degMeanAngVel = cell2mat(GBJNeg90degMeanAngVel);

GAJ45degMeanAngVel = cell2mat(GAJ45degMeanAngVel);
GAJNeg45degMeanAngVel = cell2mat(GAJNeg45degMeanAngVel);
GAJ90degMeanAngVel = cell2mat(GAJ90degMeanAngVel);
GAJNeg90degMeanAngVel = cell2mat(GAJNeg90degMeanAngVel);


%group BJ and AJ for every angle
Garound45 = [GBJ45MeandegForwardVel;GAJ45degMeanForwardVel];
GaroundNeg45 = [GBJNeg45MeandegForwardVel;GAJNeg45degMeanForwardVel];
Garound90 = [GBJ90degMeanForwardVel;GAJ90degMeanForwardVel];
GaroundNeg90 = [GBJNeg90degMeanForwardVel;GAJNeg90degMeanForwardVel];

Garound45Ang = [GBJ45MeandegAngVel;GAJ45degMeanAngVel];
GaroundNeg45Ang = [GBJNeg45MeandegAngVel;GAJNeg45degMeanAngVel];
Garound90Ang = [GBJ90degMeanAngVel;GAJ90degMeanAngVel];
GaroundNeg90Ang = [GBJNeg90degMeanAngVel;GAJNeg90degMeanAngVel];


%compute wilcoxon tests
p1 = ranksum(GBJ45MeandegForwardVel,GAJ45degMeanForwardVel);
p2 = ranksum(GBJNeg45MeandegForwardVel,GAJNeg45degMeanForwardVel);
p3 = ranksum(GBJ90degMeanForwardVel,GAJ90degMeanForwardVel);
p4 = ranksum(GBJNeg90degMeanForwardVel,GAJNeg90degMeanForwardVel);

p5 = ranksum(GBJ45MeandegAngVel,GAJ45degMeanAngVel);
p6 = ranksum(GBJNeg45MeandegAngVel,GAJNeg45degMeanAngVel);
p7 = ranksum(GBJ90degMeanAngVel,GAJ90degMeanAngVel);
p8 = ranksum(GBJNeg90degMeanAngVel,GAJNeg90degMeanAngVel);


% Plotting them
x = ones(1, length(Garound45));
figure,
%Forward velocity
subplot(2,4,1)
plot(x, Garound45(1,:),'o') %before jump
hold on
plot(2*x,Garound45(2,:),'o') %after jump
xlim([0,3]);
for i = 1:length(Garound45)
    plot([1,2],[Garound45(1,i) Garound45(2,i)],'k') %add lines to join them
end
ylabel('Forward velocity (mm/s)');
title('45 deg jumps')
text(1,max(max(Garound45)),strcat('p = ',num2str(p1))) %add p value
set(gca,'xticklabel',{[]})
subplot(2,4,2)
plot(x, GaroundNeg45(1,:),'o')
hold on
plot(2*x,GaroundNeg45(2,:),'o')
xlim([0,3]);
for i = 1:length(GaroundNeg45)
    plot([1,2],[GaroundNeg45(1,i) GaroundNeg45(2,i)],'k')
end
text(1,max(max(GaroundNeg45)),strcat('p = ',num2str(p2)))
title('-45 deg jumps')
set(gca,'xticklabel',{[]})
subplot(2,4,3)
plot(x, Garound90(1,:),'o')
hold on
plot(2*x,Garound90(2,:),'o')
xlim([0,3]);
for i = 1:length(Garound90)
    plot([1,2],[Garound90(1,i) Garound90(2,i)],'k')
end
text(1,max(max(Garound90)),strcat('p = ',num2str(p3)))
title('90 deg jumps')
set(gca,'xticklabel',{[]})
subplot(2,4,4)
plot(x, GaroundNeg90(1,:),'o')
hold on
plot(2*x,GaroundNeg90(2,:),'o')
xlim([0,3]);
for i = 1:length(GaroundNeg90)
    plot([1,2],[GaroundNeg90(1,i) GaroundNeg90(2,i)],'k')
end
text(1,max(max(GaroundNeg90)),strcat('p = ',num2str(p4)))
title('-90 deg jumps');
set(gca,'xticklabel',{[]})
%Angular velocity
subplot(2,4,5)
plot(x, Garound45Ang(1,:),'o')
hold on
plot(2*x,Garound45Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(Garound45Ang)
    plot([1,2],[Garound45Ang(1,i) Garound45Ang(2,i)],'k')
end
ylabel('Angular velocity (deg/s)');
text(1,max(max(Garound45Ang)),strcat('p = ',num2str(p5)))
set(gca,'xticklabel',{[]})
subplot(2,4,6)
plot(x, GaroundNeg45Ang(1,:),'o')
hold on
plot(2*x,GaroundNeg45Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(GaroundNeg45Ang)
    plot([1,2],[GaroundNeg45Ang(1,i) GaroundNeg45Ang(2,i)],'k')
end
text(1,max(max(GaroundNeg45Ang)),strcat('p = ',num2str(p6)))
set(gca,'xticklabel',{[]})
subplot(2,4,7)
plot(x, Garound90Ang(1,:),'o')
hold on
plot(2*x,Garound90Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(Garound90Ang)
    plot([1,2],[Garound90Ang(1,i) Garound90Ang(2,i)],'k')
end
text(1,max(max(Garound90Ang)),strcat('p = ',num2str(p7)))
set(gca,'xticklabel',{[]})
subplot(2,4,8)
plot(x, GaroundNeg90Ang(1,:),'o')
hold on
plot(2*x,GaroundNeg90Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(GaroundNeg90Ang)
    plot([1,2],[GaroundNeg90Ang(1,i) GaroundNeg90Ang(2,i)],'k')
end
text(1,max(max(GaroundNeg90Ang)),strcat('p = ',num2str(p8)))
set(gca,'xticklabel',{[]})