% This plots the "Distance to the goal" for every fly pooled


edges = [-180:10:180];

[counts] = histcounts(fly12,edges);
probabilities12 = counts./sum(counts);

degs = linspace(-180,180,length(counts));

figure,
histogram(fly7,edges,'Normalization','probability');
h.FaceColor = [0.2,0.5,1];
xlim([-180 180]);
title('Histogram of the distance to the goal');
ylabel('Probability'); xlabel('Distance to the goal (deg)');

figure,
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')
%suptitle(typeOfStim);
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Distribution of the distance to the goal');
ylabel('Probability density'); xlabel('Distance to the goal (deg)');



%% Using mean values

meanAllFlies = mean(allFlies,2);
stdAllFlies = std(allFlies,[],2);


[counts] = histcounts(meanAllFlies,edges);
probabilities = counts./sum(counts);
degs = linspace(-180,180,length(counts));

figure,histogram(meanAllFlies,edges)
xlim([-180 180]);
title('Histogram of the distance to the goal');
ylabel('Probability'); xlabel('Distance to the goal (deg)');

figure,
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')
%suptitle(typeOfStim);
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Distribution of the distance to the goal');
ylabel('Probability density'); xlabel('Distance to the goal (deg)');


figure, plot(degs,allProbabilities)
hold on
boundedline(degs,meanAllProb,stdAllProb,'k','alpha')
xlim([-180 180]);
title('Distribution of the distance to the goal for all flies');
ylabel('Probability density'); xlabel('Distance to the goal (deg)');

