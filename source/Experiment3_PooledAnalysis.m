%Experiment 3, pooled analysis

%This code pools the data from all the flies in experiment 3
%and does some pooled analyses.

clear all; close all;
%list every file inside Experiment 3 folder

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


%calculate mean ang and forw vel
meanForwardVel45 = mean(Data45.forwardVel,2);
meanForwardVelNeg45 = mean(DataNeg45.forwardVel,2);
meanForwardVel90 = mean(Data90.forwardVel,2);
meanForwardVelNeg90 = mean(DataNeg90.forwardVel,2);

meanAngVel45 = mean(Data45.angVel,2);
meanAngVelNeg45 = mean(DataNeg45.angVel,2);
meanAngVel90 = mean(Data90.angVel,2);
meanAngVelNeg90 = mean(DataNeg90.angVel,2);


%plot mean and error for each
time = linspace(-10,10,length(forwardVelAll));

figure,
subplot(1,2,1)
plot(time,meanForwardVel45)
hold on
plot(time,meanForwardVelNeg45)
plot(time,meanForwardVel90)
plot(time,meanForwardVelNeg90)
title('Mean forward velocity');
legend({'45','-45','90','-90'});
ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
subplot(1,2,2)
plot(time,meanAngVel45)
hold on
plot(time,meanAngVelNeg45)
plot(time,meanAngVel90)
plot(time,meanAngVelNeg90)
title('Mean angular velocity');
legend({'45','-45','90','-90'});
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');