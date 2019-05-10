%% Activity level

AFirstLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\*\experimental flies\*\ActivityHighErrorBlockExp3.mat');
AHighErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\*\experimental flies\*\ActivityHigherErrorBlockExp4.mat');
ASecondLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\*\experimental flies\*\ActivityHighErrorBlockExp5.mat');

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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\ActivityAcrossBlocks.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\ActivityAcrossBlocks.svg'))

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

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\StateTransitions.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\StateTransitions.svg'))

%% Around jump velocities

FirstLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\*\experimental flies\*\perTrialDataHighErrorBlockExp3.mat');
HighErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\*\experimental flies\*\perTrialDataHigherErrorBlockExp4.mat');
SecondLowErrorBlock = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\*\experimental flies\*\perTrialDataHighErrorBlockExp5.mat');

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
time = linspace(-3,3,length(angVelAll));

figure('Position', [100 100 1600 900]),
subplot(1,3,1)
p1 = boundedline(time,meanAngVel90,stdAngVel90/sqrt(i),'r','alpha')
hold on
p2 = boundedline(time,meanAngVelNeg90,stdAngVelNeg90/sqrt(i),'k','alpha')
title('First low error block');
ylim([-60, 60]);xlim([-3, 3]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
plot([-3,3],[0, 0],'k','HandleVisibility','off');
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
ylim([-60, 60]);xlim([-3, 3]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
plot([-3,3],[0, 0],'k','HandleVisibility','off');
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
ylim([-60, 60]);xlim([-3, 3]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
plot([-3,3],[0, 0],'k','HandleVisibility','off');
legend([p31,p32], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\AverageVelocityChange.png'))

%% Look at the progression in angular velocity around jumps within block

firstHalfAngVel = zeros(0,0);
firstHalfjumpMag = zeros(0,0);
secondHalfAngVel = zeros(0,0);
secondHalfjumpMag = zeros(0,0);

for i = 1:length(dataFirstLE)
    firstHalfAngVel = [firstHalfAngVel,dataHE{1,i}.perTrialData.angVel(:,1:48)];
    secondHalfAngVel = [secondHalfAngVel,dataHE{1,i}.perTrialData.angVel(:,49:96)];
    firstHalfjumpMag = [firstHalfjumpMag,dataHE{1,i}.perTrialData.jumpMag(1,1:48)];
    secondHalfjumpMag = [secondHalfjumpMag,dataHE{1,i}.perTrialData.jumpMag(1,49:96)];
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
ylim([-60, 60]);xlim([-3, 3]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
plot([-3,3],[0, 0],'k','HandleVisibility','off');
legend([p11,p12], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

subplot(1,2,2)
p21 = boundedline(time,meanAngVel90HESecondHalf,stdAngVel90HESecondHalf/sqrt(i),'r','alpha')
hold on
p22 = boundedline(time,meanAngVelNeg90HESecondHalf,stdAngVelNeg90HESecondHalf/sqrt(i),'k','alpha')
title('High error block, second half');
ylim([-60, 60]);xlim([-3, 3]);
plot([0,0],[-60, 60],'--k','HandleVisibility','off');
plot([-3,3],[0, 0],'k','HandleVisibility','off');
legend([p21,p22], '90','-90');
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment8\AverageVelocityProgressionHE.png'))
