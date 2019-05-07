%This code analyses the multi-trial experiment (for now experiment 11)

close all; clear all;

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment11'
[file,path] = uigetfile();
load([path,file]);

%% Analyze closed-loop trials

LightClosedLoopIndex = find(contains(trials,'lightClosedLoop'));
DarkClosedLoopIndex = find(contains(trials,'darkClosedLoop'));

for i = 1:length(LightClosedLoopIndex)    
   [smoothedL{i},posToRadFlyL{i},degsFlyDistMovingL{i},probabilitiesDistMovingL{i}] = runClosedLoopAnalysis(rawData{LightClosedLoopIndex(i)},path,strcat(trials{LightClosedLoopIndex(i)},num2str(LightClosedLoopIndex(i))));  
   [smoothedD{i},posToRadFlyD{i},degsFlyDistMovingD{i},probabilitiesDistMovingD{i}] = runClosedLoopAnalysis(rawData{DarkClosedLoopIndex(i)},path,strcat(trials{DarkClosedLoopIndex(i)},num2str(DarkClosedLoopIndex(i))));      
end


% Compare light stripe and dark stripe effects

%(1) Distance to the goal for light vs dark stripes
Max = [max([probabilitiesDistMovingL{:}]),max([probabilitiesDistMovingD{:}])];
Max = max(Max);

figure,
subplot(1,2,1)
for i = 1:length(LightClosedLoopIndex)
    plot(degsFlyDistMovingL{i},probabilitiesDistMovingL{i})
    hold on
end
ylim([0 Max+0.05]);
title('Dist2goal distribution for light stripes');
subplot(1,2,2)
for i = 1:length(DarkClosedLoopIndex)
    plot(degsFlyDistMovingD{i},probabilitiesDistMovingD{i})
    hold on
end
title('Dist2goal distribution for dark stripes');
ylim([0 Max+0.05]);
saveas(gcf,strcat(path,'Dist2goal',file(1:end-4),'.png'));

%(2) Polar histograms

figure('Position', [100 100 1400 400]),
circedges = [0:20:360];
circedges = deg2rad(circedges);
for i = 1:length(LightClosedLoopIndex)
subplot(1,length(LightClosedLoopIndex),i)
polarhistogram(posToRadFlyL{i}(smoothedL{1,i}.xVel>1),circedges,'Normalization','probability','FaceColor',[0,0,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
end
suptitle('Polar distribution for light stripes');
saveas(gcf,strcat(path,'PolardistL',file(1:end-4),'.png'));


figure('Position', [100 100 1400 400]),
for i = 1:length(DarkClosedLoopIndex)
subplot(1,length(DarkClosedLoopIndex),i)
polarhistogram(posToRadFlyD{i}(smoothedD{1,i}.xVel>1),circedges,'Normalization','probability','FaceColor',[0,0,0],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
end
suptitle('Polar distribution for dark stripes');
saveas(gcf,strcat(path,'PolardistD',file(1:end-4),'.png'));
%% Analyze velocity open-loop trials

%Divide in the different trial types
FastClockOpenLoopIndex = find(contains(trials,'fastClockwise'));
FastCounterclockOpenLoopIndex = find(contains(trials,'fastCounter'));
SlowClockOpenLoopIndex = find(contains(trials,'slowClockwise'));
SlowCounterclockOpenLoopIndex = find(contains(trials,'slowCounter'));

for i = 1:length(FastClockOpenLoopIndex)    
    [smoothedFC{i},remapPosToDegFC{i},flyPos180FC{i}] = runOpenLoopAnalysis(rawData{FastClockOpenLoopIndex(i)},path,strcat(trials{FastClockOpenLoopIndex(i)},num2str(FastClockOpenLoopIndex(i))));
    [smoothedFCC{i},remapPosToDegFCC{i},flyPos180FCC{i}] = runOpenLoopAnalysis(rawData{FastCounterclockOpenLoopIndex(i)},path,strcat(trials{FastCounterclockOpenLoopIndex(i)},num2str(FastCounterclockOpenLoopIndex(i))));
    [smoothedSC{i},remapPosToDegSC{i},flyPos180SC{i}] = runOpenLoopAnalysis(rawData{SlowClockOpenLoopIndex(i)},path,strcat(trials{SlowClockOpenLoopIndex(i)},num2str(SlowClockOpenLoopIndex(i))));
    [smoothedSCC{i},remapPosToDegSCC{i},flyPos180SCC{i}] = runOpenLoopAnalysis(rawData{SlowCounterclockOpenLoopIndex(i)},path,strcat(trials{SlowCounterclockOpenLoopIndex(i)},num2str(SlowCounterclockOpenLoopIndex(i))));
end

figure('Position', [100 100 1400 900]),
subplot(2,2,1)
for i = 1:length(FastClockOpenLoopIndex)
    plot(remapPosToDegFC{1,i}(5:147),smoothedFC{1,i}.angularVel(5:147))
    angVelFC(:,i) = smoothedFC{1,i}.angularVel(5:147);
    hold on
end
title('Angular Velocity for fast clockwise bars');
ylim([-150, 150]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-150 150],'Color','k','LineStyle','--');
%add mean
meanAngVelFC = mean(angVelFC,2);
plot(remapPosToDegFC{1,1}(5:147),meanAngVelFC,'k','LineWidth',2)

subplot(2,2,2)
for i = 1:length(FastCounterclockOpenLoopIndex)
    plot(remapPosToDegFCC{1,i}(1:142),smoothedFCC{1,i}.angularVel(1:142))
    angVelFCC(:,i) = smoothedFCC{1,i}.angularVel(1:142);
    hold on
end
title('Angular Velocity for fast counterclockwise bars');
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-150 150],'Color','k','LineStyle','--');
ylim([-150, 150]);
%add mean
meanAngVelFCC = mean(angVelFCC,2);
plot(remapPosToDegFCC{1,1}(1:142),meanAngVelFCC,'k','LineWidth',2)

subplot(2,2,3)
for i = 1:length(SlowClockOpenLoopIndex)
    plot(remapPosToDegSC{1,i}(11:442),smoothedSC{1,i}.angularVel(11:442))
    angVelSC(:,i) = smoothedSC{1,i}.angularVel(11:442);
    hold on
end
title('Angular Velocity for slow clockwise bars');
ylim([-150, 150]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-150 150],'Color','k','LineStyle','--');
%add mean
meanAngVelSC = mean(angVelSC,2);
plot(remapPosToDegSC{1,1}(11:442),meanAngVelSC,'k','LineWidth',2)

subplot(2,2,4)
for i = 1:length(SlowCounterclockOpenLoopIndex)
    plot(remapPosToDegSCC{1,i}(1:424),smoothedSCC{1,i}.angularVel(1:424))
    angVelSCC(:,i) = smoothedSCC{1,i}.angularVel(1:424);
    hold on
end
title('Angular Velocity for slow counterclockwise bars');
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-150 150],'Color','k','LineStyle','--');
ylim([-150, 150]);
%add mean
meanAngVelSCC = mean(angVelSCC,2);
plot(remapPosToDegSCC{1,1}(1:424),meanAngVelSCC,'k','LineWidth',2)

saveas(gcf,strcat(path,'AngVelOpenLoop',file(1:end-4),'.png'));


%% 

%Adding Bahl's analyses

%(1) P = (Rcw+Rccw)/2
figure
subplot(1,2,1)
for i = 1:length(FastClockOpenLoopIndex)
    Pfast(:,i) = (smoothedFC{1,i}.angularVel(5:142)+smoothedFCC{1,i}.angularVel(5:142))/2;
    plot(remapPosToDegFC{1,i}(5:142),Pfast(:,i))
    hold on
end
title('P for fast stimuli');
ylim([-150, 150]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-150 150],'Color','k','LineStyle','--');
%add mean
meanPfast = mean(Pfast,2);
plot(remapPosToDegFC{1,1}(5:142),meanPfast,'k','LineWidth',2)

subplot(1,2,2)
for i = 1:length(FastClockOpenLoopIndex)
    Pslow(:,i) = (smoothedSC{1,i}.angularVel(11:424)+smoothedSCC{1,i}.angularVel(11:424))/2;
    plot(remapPosToDegSC{1,i}(11:424),Pslow(:,i))
    hold on
end
title('P for slow stimuli');
ylim([-150, 150]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-150 150],'Color','k','LineStyle','--');
%add mean
meanPslow = mean(Pslow,2);
plot(remapPosToDegSC{1,1}(11:424),meanPslow,'k','LineWidth',2)


%(2) M = (Rcw-Rccw)/2
figure
subplot(1,2,1)
for i = 1:length(FastClockOpenLoopIndex)
    Mfast(:,i) = (smoothedFC{1,i}.angularVel(5:142)-smoothedFCC{1,i}.angularVel(5:142))/2;
    plot(remapPosToDegFC{1,i}(5:142),Mfast(:,i))
    hold on
end
title('M for fast stimuli');
ylim([-150, 150]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-150 150],'Color','k','LineStyle','--');
%add mean
meanMfast = mean(Mfast,2);
plot(remapPosToDegFC{1,1}(5:142),meanMfast,'k','LineWidth',2)

subplot(1,2,2)
for i = 1:length(FastClockOpenLoopIndex)
    Mslow(:,i) = (smoothedSC{1,i}.angularVel(11:424)-smoothedSCC{1,i}.angularVel(11:424))/2;
    plot(remapPosToDegSC{1,i}(11:424),Mslow(:,i))
    hold on
end
title('M for slow stimuli');
ylim([-150, 150]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-150 150],'Color','k','LineStyle','--');
%add mean
meanMslow = mean(Mslow,2);
plot(remapPosToDegSC{1,1}(11:424),meanMslow,'k','LineWidth',2)