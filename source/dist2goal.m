%Code to calculate the distance to the goal for each fly...

AngPos = {};
%create one vector of pos for each fly for every trial
for i = 1:length(angPosAllFlies)
    AngPos{i} = [angPosAllFlies{1,i}{1,1}, angPosAllFlies{1,i}{2,1}, angPosAllFlies{1,i}{3,1}, angPosAllFlies{1,i}{4,1}, angPosAllFlies{1,i}{5,1}, angPosAllFlies{1,i}{6,1}, angPosAllFlies{1,i}{7,1}, angPosAllFlies{1,i}{8,1}, angPosAllFlies{1,i}{9,1}, angPosAllFlies{1,i}{10,1}];
end

%wrap to 180
AngPos180 = cellfun(@wrapTo180, AngPos,'UniformOutput',0);

%calculate the mean which will be the goal for each fly
Goals = cellfun(@mean,AngPos180);

%generate a cell with distance to the goal
for i = 1:length(Goals)
    dist{i} = AngPos180{i}-Goals(i);
end

%plot the distance to the goal for each fly
edges = [-180:10:180];

%concatenate the distance for all flies
Distance = [dist{1},dist{2},dist{3},dist{4},dist{5},dist{6},dist{7},dist{8},dist{9},dist{10},dist{11},dist{12}];
stdAllFlies = std(Distance,[],2);
[allcounts] = histcounts(Distance,edges);
Allprobabilities = allcounts./sum(allcounts);
alldegs = linspace(-180,180,length(allcounts));

figure,   
for i = 1:length(Goals)
    [counts{i}] = histcounts(dist{i},edges);
    probabilities{i} = counts{i}./sum(counts{i});
    degs{i} = linspace(-180,180,length(counts{i}));
    hold on
    plot(degs{i},probabilities{i})
end

allFliesprob = [probabilities{1,1};probabilities{1,2};probabilities{1,3};probabilities{1,4};probabilities{1,5};probabilities{1,6};probabilities{1,7};probabilities{1,8};probabilities{1,9};probabilities{1,10};probabilities{1,11};probabilities{1,12}];
error = std(allFliesprob);

hold on
%plot(alldegs,Allprobabilities,'k','LineWidth',3)
boundedline(alldegs,Allprobabilities,error,'k','alpha')
xlabel('Distance to goal (deg)');
ylabel('Probability density');
xlim([-180 180]);
ylim([0, 0.25]);

saveas(gcf,'realdist2goal.svg');




