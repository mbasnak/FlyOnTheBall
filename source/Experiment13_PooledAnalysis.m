%Experiment 5, pooled analysis

%This code pools the data from all the flies in experiment 5
%and does some pooled analyses.

clear all; close all;


%% Activity level

AFirstLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\ActivityHighErrorBlockExp3.mat');
AHighErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\ActivityHighErrorBlockExp4.mat');
ASecondLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\ActivityHighErrorBlockExp5.mat');

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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\ActivityAcrossBlocks.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\ActivityAcrossBlocks.svg'))

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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\StateTransitions.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\StateTransitions.svg'))

%% Around jump velocities

FirstLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\perTrialDataHighErrorBlockExp3.mat');
HighErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\perTrialDataHighErrorBlockExp4.mat');
SecondLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\perTrialDataHighErrorBlockExp5.mat');

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
time = linspace(-15,15,size(angVelAll,1));

figure('Position', [100 100 1600 900]),
subplot(1,3,1)
p1 = boundedline(time,meanAngVel90,stdAngVel90/sqrt(i),'r','alpha')
hold on
p2 = boundedline(time,meanAngVelNeg90,stdAngVelNeg90/sqrt(i),'k','alpha')
title('First low error block');
ylim([-30, 30]);xlim([-8, 8]);
plot([0,0],[-30, 30],'--k','HandleVisibility','off');
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
ylim([-30, 30]);xlim([-8, 8]);
plot([0,0],[-30, 30],'--k','HandleVisibility','off');
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
ylim([-30, 30]);xlim([-8, 8]);
plot([0,0],[-30, 30],'--k','HandleVisibility','off');
plot([-8,8],[0, 0],'k','HandleVisibility','off');
legend([p31,p32], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\AverageVelocityChange.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\AverageVelocityChange.svg'))

%% distance to the goal

FilesFirstLE = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\goalsHighErrorBlockExp3.mat');
FilesHE = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\goalsHighErrorBlockExp4.mat');
FilesSecondLE = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\goalsHighErrorBlockExp5.mat');

%load every file with name %perTrialData, adding to the name the date and
%fly, reading the folder it came from

for i = 1:length(FilesHE)
    goalsFirstLE{i} = load(strcat(FilesFirstLE(i).folder,'\',FilesFirstLE(i).name));
    goalsHE{i} = load(strcat(FilesHE(i).folder,'\',FilesHE(i).name));
    goalsSecondLE{i} = load(strcat(FilesSecondLE(i).folder,'\',FilesSecondLE(i).name));
end

%concatenate as columns the arrays for all the variables
goalMovingAll = zeros(0,0);
devGoalMovingAll = zeros(0,0);
dist2goalMovingAll = {};
goalMovingAllHE = zeros(0,0);
devGoalMovingAllHE = zeros(0,0);
dist2goalMovingAllHE = {};
goalMovingAllSecondLE = zeros(0,0);
devGoalMovingAllSecondLE = zeros(0,0);
dist2goalMovingAllSecondLE = {};

for i = 1:length(FilesHE)
    goalMovingAll = [goalMovingAll;goalsFirstLE{1,i}.goalMoving];
    devGoalMovingAll = [devGoalMovingAll;goalsFirstLE{1,i}.devGoalMoving];
    dist2goalMovingAll{i} = goalsFirstLE{1,i}.dist2goalMoving;
    goalMovingAllHE = [goalMovingAllHE;goalsHE{1,i}.goalMoving];
    devGoalMovingAllHE = [devGoalMovingAllHE;goalsHE{1,i}.devGoalMoving];
    dist2goalMovingAllHE{i} = goalsHE{1,i}.dist2goalMoving;
    goalMovingAllSecondLE = [goalMovingAllSecondLE;goalsSecondLE{1,i}.goalMoving];
    devGoalMovingAllSecondLE = [devGoalMovingAllSecondLE;goalsSecondLE{1,i}.devGoalMoving];
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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\GoalDist.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\GoalDist.svg'))


%Look at the goal per flie across blocks
allGoals = [goalMovingAll,goalMovingAllHE,goalMovingAllSecondLE];
allDevGoals = [devGoalMovingAll,devGoalMovingAllHE,devGoalMovingAllSecondLE];
colors=colormap(jet(length(allGoals)));
figure,
for i = 1:8
    plot(rad2deg(allGoals(i,:)),'-o','color',colors(i,:),'MarkerFaceColor',colors(i,:))
    hold on
end
ylabel('Goal (deg)');
xlim([0,4]); ylim([-180,180]);
xticks([1 2 3])
xticklabels({'Block 1','Block 2','Block 3'})

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\GoalDistAcrossBlocks.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\GoalDistAcrossBlocks.svg'))

%Make individual within fly across blocks goal plots
for i = 1:8
figure,
plot(rad2deg(allGoals(i,:)),'-ko')
% hold on
% plot(rad2deg(allGoals(i,:)),'k')
ylabel('Goal (deg)');
xlim([0,4]);
ylim([-180 180]);
xticks([1 2 3])
xticklabels({'Block 1','Block 2','Block 3'})
title('Goal evolution')
saveas(gcf,['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\GoalEv',num2str(i),'.png'])
end

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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\Dist2GoalAllTrials.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\Dist2GoalAllTrials.svg'))


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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\Dist2goal.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\Dist2goal.svg'))


%% Look at the goal distribution whithin each block and across blocks

GoalEvFirstLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\goalTrialLowErrorBlockExp3.mat');
GoalEvHighErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\goalTrialHighErrorBlockExp4.mat');
GoalEvSecondLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\*\experimental flies\*\goalTrialLowErrorBlockExp5.mat');


for i = 1:length(GoalEvHighErrorBlock)
    
    shortDataFirstLE{i} = load(strcat(GoalEvFirstLowErrorBlock(i).folder,'\',GoalEvFirstLowErrorBlock(i).name));
    TrialGoalFirstBlock{i} = shortDataFirstLE{1,i}.goalMovingTrial;    
    
    shortDataHE{i} = load(strcat(GoalEvHighErrorBlock(i).folder,'\',GoalEvHighErrorBlock(i).name));
    TrialGoalSecondBlock{i} = shortDataHE{1,i}.goalMovingTrial;
      
    shortDataSecondLE{i} = load(strcat(GoalEvFirstLowErrorBlock(i).folder,'\',GoalEvSecondLowErrorBlock(i).name));
    TrialGoalThirdBlock{i} = shortDataSecondLE{1,i}.goalMovingTrial;
    
end

 
for i = 1:length(GoalEvHighErrorBlock)
    
    figure ('Position', [200 100 800 1000]),
    subplot(3,1,1)
    plot(TrialGoalFirstBlock{1,i},'k')
    hold on
    plot(TrialGoalFirstBlock{1,i},'ko')
    title('First low error block');
    ylim([-180 180]); xlim([1 10]);
    ylabel({'Heading goal', '[circ_mean of angular position] (deg)'});
    xlabel('Time');
    
    subplot(3,1,2)
    plot(TrialGoalSecondBlock{1,i},'r')
    hold on
    plot(TrialGoalSecondBlock{1,i},'ro')
    title('High error block');
    ylim([-180 180]); xlim([1 10]);
    ylabel({'Heading goal', '[circ_mean of angular position] (deg)'});
    xlabel('Time');  
    
    subplot(3,1,3)
    plot(TrialGoalThirdBlock{1,i},'k')
    hold on
    plot(TrialGoalThirdBlock{1,i},'ko')
    title('Second low error block');
    ylim([-180 180]); xlim([1 10]);
    ylabel({'Heading goal', '[circ_mean of angular position] (deg)'});
    xlabel('Time');
    
    saveas(gcf,['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\TrialGoal',num2str(i),'.png'])
end



