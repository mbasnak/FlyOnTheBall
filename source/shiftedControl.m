function shiftedControl(distance,data,trials,rewardDimension)

%Tot's suggestion of shifting the data randomly and then looking at the
%median effect then
clear shiftedData shiftedMean

%I need to define the number of iterations for this and repeat it that
%amount of times
nIterations = 50;

for i = 1:nIterations
    %sample a vector of random shifts, of length equal to the number of trials, from 0 to total distance (48
    %dimensions)
    possibleShifts(i,:) = [-length(distance):1:length(distance)];
    randomDeltaT(i,:) = randsample(possibleShifts(i,:),length(data),true);
    %shift the actual data by that amount
    for j = 1:length(trials)
        shiftedData(i,:,j) = circshift(data(trials(j),:),randomDeltaT(i,:));
    end
    %get the median effect
    shiftedMean(i,:) = nanmean(shiftedData(i,:,:),3);
end

ShiftedMean = mean(shiftedMean);
stdShiftedMean = std(shiftedMean);

%get the confidence interval for the data
SEM = stdShiftedMean/sqrt(nIterations);               % Standard Error
ts = tinv([0.025  0.975],nIterations-1);      % T-Score
CIl = ShiftedMean + ts(1)*SEM;                      % Confidence Intervals
CIh = ShiftedMean + ts(2)*SEM;                      % Confidence Intervals
CI = [CIl',CIh'];

figure,
[opT] = boundedline(distance,nanmean(data(trials,:)),std(data(trials,:)),'-ro','alpha');
hold on
[prT] = boundedline(distance,ShiftedMean,stdShiftedMean,'-bo','alpha');
%legend([opT,prT],'Actual data', 'Shifted control')
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('Forward velocity (mm/s)');
xlabel('Distance (mm)');


figure,
[opT] = boundedline(distance,nanmean(data(trials,:)),std(data(trials,:)),'-ro','alpha');
hold on
plot(distance,ShiftedMean,'-bo')
plot(distance,CIl,'-b')
plot(distance,CIh,'-b')
%legend([opT,prT],'Actual data', 'Shifted control')
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('Forward velocity (mm/s)');
xlabel('Distance (mm)');


end