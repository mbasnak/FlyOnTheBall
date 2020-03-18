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
probL = cell2mat(probabilitiesDistMovingL);
probL2 = reshape(probL,  18, length(probL)/18);
meanProbL = mean(probL2,2);
degs = degsFlyDistMovingL{1,1};
plot(degs,meanProbL,'k','LineWidth',2)

subplot(1,2,2)
for i = 1:length(DarkClosedLoopIndex)
    plot(degsFlyDistMovingD{i},probabilitiesDistMovingD{i})
    hold on
end
title('Dist2goal distribution for dark stripes');
ylim([0 Max+0.05]);
probD = cell2mat(probabilitiesDistMovingD);
probD2 = reshape(probD,  18, length(probD)/18);
meanProbD = mean(probD2,2);
plot(degs,meanProbD,'k','LineWidth',2)

saveas(gcf,strcat(path,'Dist2goal',file(1:end-4),'.png'));

%(2) Polar histograms

figure('Position', [100 100 1400 400]),
circedges = [0:20:360];
circedges = deg2rad(circedges);
for i = 1:length(LightClosedLoopIndex)
subplot(1,length(LightClosedLoopIndex),i)
polarhistogram(posToRadFlyL{i}(smoothedL{1,i}.xVel>0.5),circedges,'Normalization','probability','FaceColor',[0,0,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
end
suptitle('Polar distribution for light stripes');
saveas(gcf,strcat(path,'PolardistL',file(1:end-4),'.png'));
saveas(gcf,strcat(path,'PolardistL',file(1:end-4),'.svg'));

figure('Position', [100 100 1400 400]),
for i = 1:length(DarkClosedLoopIndex)
subplot(1,length(DarkClosedLoopIndex),i)
polarhistogram(posToRadFlyD{i}(smoothedD{1,i}.xVel>0.5),circedges,'Normalization','probability','FaceColor',[0,0,0],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
end
suptitle('Polar distribution for dark stripes');

saveas(gcf,strcat(path,'PolardistD',file(1:end-4),'.png'));
saveas(gcf,strcat(path,'PolardistD',file(1:end-4),'.svg'));
%% Analyze velocity open-loop trials

%Divide in the different trial types
FastClockOpenLoopIndex = find(contains(trials,'fastClockwise'));
FastCounterclockOpenLoopIndex = find(contains(trials,'fastCounter'));
FClockOpenLoopIndex = find(contains(trials,'clockwiseF'));
FCounterclockOpenLoopIndex = find(contains(trials,'counterclockwiseF'));

for i = 1:length(FastClockOpenLoopIndex)    
    [smoothedFC{i},remapPosToDegFC{i},flyPos180FC{i}] = runOpenLoopAnalysis(rawData{FastClockOpenLoopIndex(i)},path,strcat(trials{FastClockOpenLoopIndex(i)},num2str(FastClockOpenLoopIndex(i))));
    [smoothedFCC{i},remapPosToDegFCC{i},flyPos180FCC{i}] = runOpenLoopAnalysis(rawData{FastCounterclockOpenLoopIndex(i)},path,strcat(trials{FastCounterclockOpenLoopIndex(i)},num2str(FastCounterclockOpenLoopIndex(i))));
    [smoothedSC{i},remapPosToDegSC{i},flyPos180SC{i}] = runOpenLoopAnalysis(rawData{FClockOpenLoopIndex(i)},path,strcat(trials{FClockOpenLoopIndex(i)},num2str(FClockOpenLoopIndex(i))));
    [smoothedSCC{i},remapPosToDegSCC{i},flyPos180SCC{i}] = runOpenLoopAnalysis(rawData{FCounterclockOpenLoopIndex(i)},path,strcat(trials{FCounterclockOpenLoopIndex(i)},num2str(FCounterclockOpenLoopIndex(i))));
end


for i = 1:length(FastClockOpenLoopIndex)
    angVelFC(:,i) = smoothedFC{1,i}.angularVel;
    angVelFCC(:,i) = smoothedFCC{1,i}.angularVel;
end
meanAngVelFC = mean(angVelFC,2);
meanAngVelFCC = mean(angVelFCC,2);
stdFC = std(angVelFC,[],2);
stdFCC = std(angVelFCC,[],2);

figure('Position',[200 200 1000 600]),
subplot(1,2,1)
p1 = boundedline(remapPosToDegFC{1,1}(4:end),meanAngVelFC(4:end),stdFC(4:end)/sqrt(i),'k','alpha')
hold on
p2 = boundedline(remapPosToDegFCC{1,1}(1:end-3),meanAngVelFCC(1:end-3),stdFCC(1:end-3)/sqrt(i),'r','alpha')
plot([0 0], [-60 60],'k','handleVisibility','off')
plot([-180 180], [0 0], 'k','handleVisibility','off')
ylim([-60 60]); xlim([-180 180]);
title('Open-loop tracking responses');
ylabel('Angular velocity (deg/s)'); xlabel('Bar position (deg)');
legend([p1,p2],'Clockwise turns', 'Counterclockwise turns');


% Plot response to Fourier bars

for i = 1:length(FastClockOpenLoopIndex)
    angVelSC(:,i) = smoothedSC{1,i}.angularVel;
    angVelSCC(:,i) = smoothedSCC{1,i}.angularVel;
end
meanAngVelSC = mean(angVelSC,2);
meanAngVelSCC = mean(angVelSCC,2);
stdSC = std(angVelSC,[],2);
stdSCC = std(angVelSCC,[],2);

subplot(1,2,2)
p1 = boundedline(remapPosToDegSC{1,1}(4:end),meanAngVelSC(4:end),stdSC(4:end)/sqrt(i),'k','alpha')
hold on
p2 = boundedline(remapPosToDegSCC{1,1}(1:end-3),meanAngVelSCC(1:end-3),stdSCC(1:end-3)/sqrt(i),'r','alpha')
plot([0 0], [-60 60],'k','handleVisibility','off')
plot([-180 180], [0 0], 'k','handleVisibility','off')
ylim([-60 60]); xlim([-180 180]);
title('Open-loop tracking responses with Fourier bar');
ylabel('Angular velocity (deg/s)'); xlabel('Bar position (deg)');
legend([p1,p2],'Clockwise turns', 'Counterclockwise turns');

saveas(gcf,strcat(path,'AngVelOpenLoop',file(1:end-4),'.png'));

%% Adding Bahl's analyses

%(1) P = (Rcw+Rccw)/2

figure
subplot(1,2,1)
for i = 1:length(FastClockOpenLoopIndex)
    Pfast(:,i) = (smoothedFC{1,i}.angularVel(4:end) + flip(smoothedFCC{1,i}.angularVel(1:end-3)))/2;
    plot(remapPosToDegFC{1,i}(4:end),Pfast(:,i))
    hold on
end
title('P for open-loop bar stimuli');
ylim([-50, 50]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-50 50],'Color','k','LineStyle','--');
meanPfast = mean(Pfast,2);
plot(remapPosToDegFC{1,1}(4:end),meanPfast,'k','LineWidth',2)


%(2) M = (Rcw-Rccw)/2
subplot(1,2,2)
for i = 1:length(FastClockOpenLoopIndex)
    Mfast(:,i) = (flip(smoothedFCC{1,i}.angularVel(1:end-3)) - smoothedFC{1,i}.angularVel(4:end))/2;;
    plot(remapPosToDegFC{1,i}(4:end),Mfast(:,i))
    hold on
end
title('M for open-loop stimuli');
ylim([-50, 50]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-50 50],'Color','k','LineStyle','--');
meanMfast = mean(Mfast,2);
plot(remapPosToDegFC{1,1}(4:end),meanMfast,'k','LineWidth',2)

%%  For the fourier bar

figure
subplot(1,2,1)
for i = 1:length(FastClockOpenLoopIndex)
    Pfourier(:,i) = (smoothedSC{1,i}.angularVel(4:end) + flip(smoothedSCC{1,i}.angularVel(1:end-3)))/2;
    plot(remapPosToDegSC{1,i}(4:end),Pfourier(:,i))
    hold on
end
title('P for fourier bar');
ylim([-50, 50]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-50 50],'Color','k','LineStyle','--');
meanPfourier = mean(Pfourier,2);
plot(remapPosToDegSC{1,1}(4:end),meanPfourier,'k','LineWidth',2)


%(2) M = (Rcw-Rccw)/2
subplot(1,2,2)
for i = 1:length(FastClockOpenLoopIndex)
    Mfourier(:,i) = (smoothedSC{1,i}.angularVel(4:end) - flip(smoothedSCC{1,i}.angularVel(1:end-3)))/2;
    plot(remapPosToDegSC{1,i}(4:end),Mfourier(:,i))
    hold on
end
title('M for fourier stimuli');
ylim([-50, 50]);
hold on
line([-180 180],[0 0],'Color','k','LineStyle','--');
line([0 0],[-50 50],'Color','k','LineStyle','--');
meanMfourier = mean(Mfourier,2);
plot(remapPosToDegSC{1,1}(4:end),meanMfourier,'k','LineWidth',2)