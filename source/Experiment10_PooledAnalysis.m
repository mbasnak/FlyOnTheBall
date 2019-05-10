%Experiment 10 pooled analysis

clear all; close all;


%% Angular Velocity

AngVelFirstMB = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\AngularVelocityClockwisedataExpNum3.mat');
AngVelYB = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\AngularVelocityClockwiseyokedDataExpNum4.mat');
AngVelSecondMB = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\AngularVelocityClockwisedataExpNum5.mat');

for i = 1:length(AngVelFirstMB)
    AngularVelFirstMB{i} = load(strcat(AngVelFirstMB(i).folder,'\',AngVelFirstMB(i).name));
    AngularVelYB{i} = load(strcat(AngVelYB(i).folder,'\',AngVelYB(i).name));
    AngularVelSecondMB{i} = load(strcat(AngVelSecondMB(i).folder,'\',AngVelSecondMB(i).name));
end


AngVelFirstMBCC = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\AngularVelocityCounterclockwisedataExpNum3.mat');
AngVelYBCC = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\AngularVelocityCounterclockwiseyokedDataExpNum4.mat');
AngVelSecondMBCC = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\AngularVelocityCounterclockwisedataExpNum5.mat');

for i = 1:length(AngVelFirstMBCC)
    AngularVelFirstMBCC{i} = load(strcat(AngVelFirstMBCC(i).folder,'\',AngVelFirstMBCC(i).name));
    AngularVelYBCC{i} = load(strcat(AngVelYBCC(i).folder,'\',AngVelYBCC(i).name));
    AngularVelSecondMBCC{i} = load(strcat(AngVelSecondMBCC(i).folder,'\',AngVelSecondMBCC(i).name));
end

angVelAll = zeros(0,0);
angVelAllCC = zeros(0,0);
angVelAllYB = zeros(0,0);
angVelAllYBCC = zeros(0,0);
angVelAllSecondMB = zeros(0,0);
angVelAllSecondMBCC = zeros(0,0);

for i = 1:length(AngVelFirstMBCC)
    angVelAll = [angVelAll, AngularVelFirstMB{1,i}.angularVelocityClockwise];
    angVelAllCC = [angVelAllCC, AngularVelFirstMBCC{1,i}.angularVelocityCounterclockwise];
    
    angVelAllYB = [angVelAllYB, AngularVelYB{1,i}.angularVelocityClockwise];
    angVelAllYBCC = [angVelAllYBCC, AngularVelYBCC{1,i}.angularVelocityCounterclockwise];
    
    angVelAllSecondMB = [angVelAllSecondMB, AngularVelSecondMB{1,i}.angularVelocityClockwise];
    angVelAllSecondMBCC = [angVelAllSecondMBCC, AngularVelSecondMBCC{1,i}.angularVelocityCounterclockwise];       
end

%calculate mean and std ang vel
%(1) first block
meanAngVelC = mean(angVelAll,2);
stdAngVelC = std(angVelAll,[],2);
meanAngVelCC = mean(angVelAllCC,2);
stdAngVelCC = std(angVelAllCC,[],2);

meanAngVelYB = mean(angVelAllYB,2);
stdAngVelYB = std(angVelAllYB,[],2);
meanAngVelYBCC = mean(angVelAllYBCC,2);
stdAngVelYBCC = std(angVelAllYBCC,[],2);

meanAngVelSecondMB = mean(angVelAllSecondMB,2);
stdAngVelSecondMB = std(angVelAllSecondMB,[],2);
meanAngVelSecondMBCC = mean(angVelAllSecondMBCC,2);
stdAngVelSecondMBCC = std(angVelAllSecondMBCC,[],2);


%plot mean and error for each
time = linspace(0,3,size(angVelAll,1));

figure('Position', [100 100 1600 900]),
subplot(1,3,1)
p1 = boundedline(time,meanAngVelC,stdAngVelC/sqrt(i),'r','alpha')
hold on
p2 = boundedline(time,meanAngVelCC,stdAngVelCC/sqrt(i),'k','alpha')
title('First master block');
ylim([-150, 150]);xlim([0, 3]);
plot([0,3],[0, 0],'--k','HandleVisibility','off');
legend([p1,p2], 'clockwise','counterclockwise');
ylabel('Angular velocity (deg/s)'); xlabel('Time (s)');

subplot(1,3,2)
p3 = boundedline(time,meanAngVelYB,stdAngVelYB/sqrt(i),'r','alpha')
hold on
p4 = boundedline(time,meanAngVelYBCC,stdAngVelYBCC/sqrt(i),'k','alpha')
title('Yoked block');
ylim([-150, 150]);xlim([0, 3]);
plot([0,3],[0, 0],'--k','HandleVisibility','off');
legend([p3,p4], 'clockwise','counterclockwise');
xlabel('Time (s)');

subplot(1,3,3)
p5 = boundedline(time,meanAngVelSecondMB,stdAngVelSecondMB/sqrt(i),'r','alpha')
hold on
p6 = boundedline(time,meanAngVelSecondMBCC,stdAngVelSecondMBCC/sqrt(i),'k','alpha')
title('Second master block');
ylim([-150, 150]);xlim([0, 3]);
plot([0,3],[0, 0],'--k','HandleVisibility','off');
legend([p5,p6], 'clockwise','counterclockwise');
xlabel('Time (s)');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\AverageVelocityChange.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\AverageVelocityChange.svg'))

%% Evolution of the response

ResponseMB = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\ResponsesdataExpNum3.mat');
ResponseYB = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\ResponsesyokedDataExpNum4.mat');
ResponseSecondMB = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\*\*\ResponsesdataExpNum5.mat');

for i = 1:length(ResponseMB)
    ResponseFirstMB{i} = load(strcat(ResponseMB(i).folder,'\',ResponseMB(i).name));
    ResponsesYB{i} = load(strcat(ResponseYB(i).folder,'\',ResponseYB(i).name));
    ResponsesSecondMB{i} = load(strcat(ResponseSecondMB(i).folder,'\',ResponseSecondMB(i).name));
end

figure('Position', [100 100 1600 900]),
subplot(1,3,1)
for i = 1:length(ResponseFirstMB)    
    plot(ResponseFirstMB{1,i}.response,'ro')
    hold on
    plot(ResponseFirstMB{1,i}.responseCC,'ko')
    hold on    
end
plot([0,30],[0, 0],'--k','HandleVisibility','off');
ylim([-350, 350]); xlim([0 30]);
ylabel('Magnitude of the optomotor response (deg/s)'); xlabel('Trial number');
title('First master block');
legend('clockwise','counterclockwise');

subplot(1,3,2)
for i = 1:length(ResponseFirstMB)    
    plot(ResponsesYB{1,i}.response,'ro')
    hold on
    plot(ResponsesYB{1,i}.responseCC,'ko')
    hold on    
end
plot([0,30],[0, 0],'--k','HandleVisibility','off');
ylim([-350, 350]); xlim([0 30]);
title('Yoked block');
xlabel('Trial number');
legend('clockwise','counterclockwise');

subplot(1,3,3)
for i = 1:length(ResponseFirstMB)    
    plot(ResponsesSecondMB{1,i}.response,'ro')
    hold on
    plot(ResponsesSecondMB{1,i}.responseCC,'ko')
    hold on    
end
plot([0,30],[0, 0],'--k','HandleVisibility','off');
ylim([-350, 350]); xlim([0, 30]);
xlabel('Trial number');
title('Second master block');
legend('clockwise','counterclockwise');
legend('clockwise','counterclockwise');

saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\ResponseEvolution.png'))
saveas(gcf,strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment10\ResponseEvolution.svg'))
