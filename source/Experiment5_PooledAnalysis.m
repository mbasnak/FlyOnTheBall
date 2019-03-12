%Experiment 5, pooled analysis

%This code pools the data from all the flies in experiment 5
%and does some pooled analyses.

clear all; close all;


%% Activity level

AFirstLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\ActivityLowErrorBlockExp3.mat');
AHighErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\ActivityHighErrorBlockExp4.mat');
ASecondLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\*\experimental flies\*\ActivityLowErrorBlockExp5.mat');

for i = 1:length(AFirstLowErrorBlock)
    ActFirstLE{i} = load(strcat(AFirstLowErrorBlock(i).folder,'\',AFirstLowErrorBlock(i).name));
    ActHE{i} = load(strcat(AHighErrorBlock(i).folder,'\',AHighErrorBlock(i).name));
    ActSecondLE{i} = load(strcat(ASecondLowErrorBlock(i).folder,'\',ASecondLowErrorBlock(i).name));
end

for i = 1:length(ActFirstLE)
    ActivityPercentageFirstLE(i) = ActFirstLE{1,i}.percentageActivity;
    ActivityPercentageHE(i) = ActHE{1,i}.percentageActivity;
    ActivityPercentageSecondLE(i) = ActSecondLE{1,i}.percentageActivity;
end

AllActivityPercentage = [ActivityPercentageFirstLE;ActivityPercentageHE;ActivityPercentageSecondLE];

%Plot the activity percentage per block
figure, 
subplot(1,2,1)
plot(AllActivityPercentage,'o')
hold on, plot(AllActivityPercentage)
ylabel('Activity (%)');
xlim([0,4]); ylim([0,100]);
xticks([1 2 3])
xticklabels({'Block 1','Block 2','Block 3'})
subplot(1,2,2)
boxplot((AllActivityPercentage'))
ylim([0,100]);
xticklabels({'Block 1','Block 2','Block 3'})
suptitle('Activity percentage across blocks');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\ActivityAcrossBlocks.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\ActivityAcrossBlocks.svg'))

%% State transitions
%I'm defining a state transition as a changing from 'quiescence' to
%'movement' of from 'movement' to 'quiescence'
activityAll = zeros(0,0);
activityAllHE = zeros(0,0);
activityAllSecondLE = zeros(0,0);

for i = 1:length(ActFirstLE)
    activityAll = [activityAll,ActFirstLE{1,i}.activity];
    activityAllHE = [activityAllHE,ActHE{1,i}.activity];
    activityAllSecondLE = [activityAllSecondLE,ActSecondLE{1,i}.activity];
end

%I will determine the state transitions by generating vectors that take the
%difference between two consecutive points. If the state remains the same,
%the cell will read 0, otherwise it will read 1 or -1.
stateTransFirstLE = diff(activityAll);
stateTransHE = diff(activityAllHE);
stateTransSecondLE = diff(activityAllSecondLE);

%we are now going to set every non-zero cell to one, and those are the
%cells where the state transitioned
stateTransFirstLE(stateTransFirstLE~=0)=1;
stateTransHE(stateTransHE~=0)=1;
stateTransSecondLE(stateTransSecondLE~=0)=1;

%The probability of state transitions for each block should be the mean
allStateTrans = [mean(stateTransFirstLE);mean(stateTransHE);mean(stateTransSecondLE)];

figure, 
subplot(1,2,1)
plot(allStateTrans,'o')
hold on, plot(allStateTrans)
ylabel('Probability of state transition');
xlim([0,4]); ylim([0,0.03]);
xticks([1 2 3])
xticklabels({'Block 1','Block 2','Block 3'})
subplot(1,2,2)
boxplot(allStateTrans')
ylim([0,0.03]);
xticklabels({'Block 1','Block 2','Block 3'})
suptitle('State transitions across blocks');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\StateTransitions.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\StateTransitions.svg'))

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
ylim([-60, 60]);xlim([-8, 8]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
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
ylim([-60, 60]);xlim([-8, 8]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
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
ylim([-60, 60]);xlim([-8, 8]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p31,p32], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\AverageVelocityChange.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\AverageVelocityChange.svg'))

%% Look at the progression in angular velocity around jumps within block

firstHalfAngVel = zeros(0,0);
firstHalfjumpMag = zeros(0,0);
secondHalfAngVel = zeros(0,0);
secondHalfjumpMag = zeros(0,0);

for i = 1:length(dataFirstLE)
    firstHalfAngVel = [firstHalfAngVel,dataHE{1,i}.perTrialData.angVel(:,1:12)];
    secondHalfAngVel = [secondHalfAngVel,dataHE{1,i}.perTrialData.angVel(:,13:24)];
    firstHalfjumpMag = [firstHalfjumpMag,dataHE{1,i}.perTrialData.jumpMag(1,1:12)];
    secondHalfjumpMag = [secondHalfjumpMag,dataHE{1,i}.perTrialData.jumpMag(1,13:24)];
end

%subset into the 2 magnitude groups
Data90HEFirstHalf.angVel = firstHalfAngVel(:,firstHalfjumpMag == 90);
DataNeg90HEFirstHalf.angVel = firstHalfAngVel(:,firstHalfjumpMag == -90);
Data90HESecondHalf.angVel = secondHalfAngVel(:,secondHalfjumpMag == 90);
DataNeg90HESecondHalf.angVel = secondHalfAngVel(:,secondHalfjumpMag == -90);

%calculate mean and std ang vel
meanAngVel90HEFirstHalf = mean(Data90HEFirstHalf.angVel,2);
stdAngVel90HEFirstHalf = std(Data90HEFirstHalf.angVel,[],2);
meanAngVelNeg90HEFirstHalf = mean(DataNeg90HEFirstHalf.angVel,2);
stdAngVelNeg90HEFirstHalf = std(DataNeg90HEFirstHalf.angVel,[],2);

meanAngVel90HESecondHalf = mean(Data90HESecondHalf.angVel,2);
stdAngVel90HESecondHalf = std(Data90HESecondHalf.angVel,[],2);
meanAngVelNeg90HESecondHalf = mean(DataNeg90HESecondHalf.angVel,2);
stdAngVelNeg90HESecondHalf = std(DataNeg90HESecondHalf.angVel,[],2);

%plot mean and error for each
figure('Position', [200 200 1200 900]),
subplot(1,2,1)
p11 = boundedline(time,meanAngVel90HEFirstHalf,stdAngVel90HEFirstHalf/sqrt(i),'r','alpha')
hold on
p12 = boundedline(time,meanAngVelNeg90HEFirstHalf,stdAngVelNeg90HEFirstHalf/sqrt(i),'k','alpha')
title('High error block, first half');
ylim([-60, 60]);xlim([-8, 8]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p11,p12], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

subplot(1,2,2)
p21 = boundedline(time,meanAngVel90HESecondHalf,stdAngVel90HESecondHalf/sqrt(i),'r','alpha')
hold on
p22 = boundedline(time,meanAngVelNeg90HESecondHalf,stdAngVelNeg90HESecondHalf/sqrt(i),'k','alpha')
title('High error block, second half');
ylim([-60, 60]);xlim([-8, 8]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p21,p22], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\AverageVelocityProgressionHE.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\AverageVelocityProgressionHE.svg'))

%% Add Jonathan's quantitative analysis:
%he takes the mean velocity from -2 to 0 and compares it with the one from
%0 to 2, for each fly, using a wilcoxon signed-rank test

jumpFrame = round(size(angVelAll,1)/2+1);

%2 sec = 50 frames, 5 sec = 150 frames

for i = 1:length(dataFirstLE)
    BJ90degMeanAngVel{i}= mean(dataFirstLE{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,dataFirstLE{1,i}.perTrialData.jumpMag==90));
    BJNeg90degMeanAngVel{i}= mean(dataFirstLE{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,dataFirstLE{1,i}.perTrialData.jumpMag==-90));
    AJ90degMeanAngVel{i}= mean(dataFirstLE{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,dataFirstLE{1,i}.perTrialData.jumpMag==90));
    AJNeg90degMeanAngVel{i}= mean(dataFirstLE{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,dataFirstLE{1,i}.perTrialData.jumpMag==-90));

    BJ90degMeanAngVelHE{i}= mean(dataHE{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,dataHE{1,i}.perTrialData.jumpMag==90));
    BJNeg90degMeanAngVelHE{i}= mean(dataHE{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,dataHE{1,i}.perTrialData.jumpMag==-90));
    AJ90degMeanAngVelHE{i}= mean(dataHE{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,dataHE{1,i}.perTrialData.jumpMag==90));
    AJNeg90degMeanAngVelHE{i}= mean(dataHE{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,dataHE{1,i}.perTrialData.jumpMag==-90));

    BJ90degMeanAngVelSecondLE{i}= mean(dataSecondLE{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,dataSecondLE{1,i}.perTrialData.jumpMag==90));
    BJNeg90degMeanAngVelSecondLE{i}= mean(dataSecondLE{1,i}.perTrialData.angVel(jumpFrame-50:jumpFrame,dataSecondLE{1,i}.perTrialData.jumpMag==-90));
    AJ90degMeanAngVelSecondLE{i}= mean(dataSecondLE{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,dataSecondLE{1,i}.perTrialData.jumpMag==90));
    AJNeg90degMeanAngVelSecondLE{i}= mean(dataSecondLE{1,i}.perTrialData.angVel(jumpFrame:jumpFrame+50,dataSecondLE{1,i}.perTrialData.jumpMag==-90));
end

%convert to doubles
BJ90degMeanAngVel = cell2mat(BJ90degMeanAngVel);
BJNeg90degMeanAngVel = cell2mat(BJNeg90degMeanAngVel);
AJ90degMeanAngVel = cell2mat(AJ90degMeanAngVel);
AJNeg90degMeanAngVel = cell2mat(AJNeg90degMeanAngVel);

BJ90degMeanAngVelHE = cell2mat(BJ90degMeanAngVelHE);
BJNeg90degMeanAngVelHE = cell2mat(BJNeg90degMeanAngVelHE);
AJ90degMeanAngVelHE = cell2mat(AJ90degMeanAngVelHE);
AJNeg90degMeanAngVelHE = cell2mat(AJNeg90degMeanAngVelHE);

BJ90degMeanAngVelSecondLE = cell2mat(BJ90degMeanAngVelSecondLE);
BJNeg90degMeanAngVelSecondLE = cell2mat(BJNeg90degMeanAngVelSecondLE);
AJ90degMeanAngVelSecondLE = cell2mat(AJ90degMeanAngVelSecondLE);
AJNeg90degMeanAngVelSecondLE = cell2mat(AJNeg90degMeanAngVelSecondLE);

%group BJ and AJ for every angle
around90Ang = [BJ90degMeanAngVel;AJ90degMeanAngVel];
aroundNeg90Ang = [BJNeg90degMeanAngVel;AJNeg90degMeanAngVel];

around90AngHE = [BJ90degMeanAngVelHE;AJ90degMeanAngVelHE];
aroundNeg90AngHE = [BJNeg90degMeanAngVelHE;AJNeg90degMeanAngVelHE];

around90AngSecondLE = [BJ90degMeanAngVelSecondLE;AJ90degMeanAngVelSecondLE];
aroundNeg90AngSecondLE = [BJNeg90degMeanAngVelSecondLE;AJNeg90degMeanAngVelSecondLE];

%compute wilcoxon tests
p1 = ranksum(BJ90degMeanAngVel,AJ90degMeanAngVel);
p2 = ranksum(BJNeg90degMeanAngVel,AJNeg90degMeanAngVel);
p3 = ranksum(BJ90degMeanAngVelHE,AJ90degMeanAngVelHE);
p4 = ranksum(BJNeg90degMeanAngVelHE,AJNeg90degMeanAngVelHE);
p5 = ranksum(BJ90degMeanAngVelSecondLE,AJ90degMeanAngVelSecondLE);
p6 = ranksum(BJNeg90degMeanAngVelSecondLE,AJNeg90degMeanAngVelSecondLE);

% Plotting them
x = ones(1, length(around90Ang));
figure('Position', [100 100 1600 900]),
%Angular velocity
subplot(2,3,1)
plot(x, around90Ang(1,:),'o')
hold on
plot(2*x,around90Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(around90Ang)
    plot([1,2],[around90Ang(1,i) around90Ang(2,i)],'k')
end
ylabel('Angular velocity (deg/s)');
text(1,max(max(around90Ang)),strcat('p = ',num2str(p1)))
set(gca,'xticklabel',{[]})
title('First low error block');
subplot(2,3,4)
plot(x, aroundNeg90Ang(1,:),'o')
hold on
plot(2*x,aroundNeg90Ang(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg90Ang)
    plot([1,2],[aroundNeg90Ang(1,i) aroundNeg90Ang(2,i)],'k')
end
text(1,max(max(aroundNeg90Ang)),strcat('p = ',num2str(p2)))
set(gca,'xticklabel',{[]})
x = ones(1, length(around90AngHE));
subplot(2,3,2)
plot(x, around90AngHE(1,:),'o')
hold on
plot(2*x,around90AngHE(2,:),'o')
xlim([0,3]);
for i = 1:length(around90AngHE)
    plot([1,2],[around90AngHE(1,i) around90AngHE(2,i)],'k')
end
ylabel('Angular velocity (deg/s)');
text(1,max(max(around90AngHE)),strcat('p = ',num2str(p3)))
set(gca,'xticklabel',{[]})
title('High error block');
subplot(2,3,5)
x = ones(1, length(aroundNeg90AngHE));
plot(x, aroundNeg90AngHE(1,:),'o')
hold on
plot(2*x,aroundNeg90AngHE(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg90AngHE)
    plot([1,2],[aroundNeg90AngHE(1,i) aroundNeg90AngHE(2,i)],'k')
end
text(1,max(max(aroundNeg90Ang)),strcat('p = ',num2str(p4)))
set(gca,'xticklabel',{[]})
subplot(2,3,3)
x = ones(1, length(around90AngSecondLE));
plot(x, around90AngSecondLE(1,:),'o')
hold on
plot(2*x,around90AngSecondLE(2,:),'o')
xlim([0,3]);
for i = 1:length(around90AngSecondLE)
    plot([1,2],[around90AngSecondLE(1,i) around90AngSecondLE(2,i)],'k')
end
text(1,max(max(around90AngSecondLE)),strcat('p = ',num2str(p5)))
set(gca,'xticklabel',{[]})
title('Second low error block');
subplot(2,3,6)
x = ones(1, length(aroundNeg90AngSecondLE));
plot(x, aroundNeg90AngSecondLE(1,:),'o')
hold on
plot(2*x,aroundNeg90AngSecondLE(2,:),'o')
xlim([0,3]);
for i = 1:length(aroundNeg90AngSecondLE)
    plot([1,2],[aroundNeg90AngSecondLE(1,i) aroundNeg90AngSecondLE(2,i)],'k')
end
text(1,max(max(aroundNeg90AngSecondLE)),strcat('p = ',num2str(p6)))
set(gca,'xticklabel',{[]})

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

suptitle('Goal distribution');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\GoalDist.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\GoalDist.svg'))


%Look at the goal per flie across blocks
allGoals = [goalMovingAll,goalMovingAllHE,goalMovingAllSecondLE];
figure, plot(rad2deg(allGoals'),'o')
hold on, plot(rad2deg(allGoals'))
ylabel('Goal (deg)');
xlim([0,4]);
xticks([1 2 3])
xticklabels({'Block 1','Block 2','Block 3'})

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\GoalDistAcrossBlocks.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\GoalDistAcrossBlocks.svg'))



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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\Dist2GoalAllTrials.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\Dist2GoalAllTrials.svg'))


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
[h1,h1bis] = boundedline(alldegsMoving,AllprobabilitiesMoving,errorMoving,'k','alpha')
set(h1, 'linewidth', 3);
set(h1bis, 'HandleVisibility','off');
hold on
[h2,h2bis] = boundedline(alldegsMovingHE,AllprobabilitiesMovingHE,errorMovingHE,'r','alpha')
set(h2, 'linewidth', 3);
set(h2bis, 'HandleVisibility','off');
[h3,h3bis] = boundedline(alldegsMovingSecondLE,AllprobabilitiesMovingSecondLE,errorMovingSecondLE,'b','alpha')
set(h3, 'linewidth', 3);
set(h3bis, 'HandleVisibility','off');
legend('First low error block','High error block','Second low error block');
xlabel('Distance to goal (deg)');
ylabel('Probability');

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

suptitle('Distance to goal 10 sec around jumps');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\Dist2Goal10secHeatMap.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment5\Dist2Goal10secHeatMap.svg'))

