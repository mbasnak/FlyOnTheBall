%% Pooling experiments

% At this point, every fly is receiving different experiments of multiple
% trials each. This file will pool every experiment together for a given
% fly, group the similar trials and plot them

close all; clear all

% go inside a certain fly's folder
% mypath = uigetdir;

% for every dataExpNum, load the rawData as a field in a struct, the trials,
%and the typeOfStim

% Read the files inside the fly's folder
% name the files that have 'dataExp' in the name
Dir = dir('dataExp*');

% Initialize empty struct
data = struct;

% load the rawData into the struct
for i = 1:length(Dir)
    data(i).rawData = load(Dir(i).name,'rawData');
end 

% If the rawData has only one trial (ie it's not an experiment), dicard it
classes = {};
for j = 1:length(data)
    classes{j} = class(data(j).rawData.rawData);    
end
% keep only the 'cell' data (ie, having more than one trial)
Data = data(cellfun(@(x) isequal(x,'cell'), classes));

% Read the names of the files we are keeping
fileNames = Dir(cellfun(@(x) isequal(x,'cell'), classes));
% load the trials pattern for each as a struct row
for i = 1:length(fileNames)
    Data(i).trials = load(fileNames(i).name,'trials');
end 

%% Subset the 3 types of stimulus and reorganize into 3 structs (or cells?)

for j = 1:length(fileNames)
    DarkBarTrials{j} = cellfun(@(x) isequal(x,'DarkClosedLoopBar'), Data(j).trials.trials);
    LightBarTrials{j} = cellfun(@(x) isequal(x,'ClosedLoopBar'), Data(j).trials.trials);
    OpenLoopGratingTrials{j} = cellfun(@(x) isequal(x,'OpenLoopGrating'), Data(j).trials.trials);
end

for i = 1:length(fileNames)  
    darkBar{i} = Data(i).rawData.rawData(DarkBarTrials{i});
    lightBar{i} = Data(i).rawData.rawData(LightBarTrials{i});
    gratings{i} = Data(i).rawData.rawData(OpenLoopGratingTrials{i});
end
    
%reshape the arrays to loose the division by experiments
DarkBar = horzcat(darkBar{:});
LightBar = horzcat(lightBar{:});
Gratings = horzcat(gratings{:});

%% For dark closed-loop bars

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;

for i = 1:size(DarkBar,2)
 darkB.xPanelVolts{i} =  DarkBar{i} (:,xPanels); 
 VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
 maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
 darkB.xPanelPos{i} = round ((darkB.xPanelVolts{i}  * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels

 darkB.yPanelVolts{i} =  DarkBar{i} (:, yPanels);
 VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
 maxValY = 1;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2
 darkB.yPanelPos{i} = round ((darkB.yPanelVolts{i}  * maxValY) /VOLTAGE_RANGE);

%FicTrac data
 darkB.ficTracAngularPosition{i} = DarkBar{i} ( : , headingFly); 
 darkB.ficTracIntx{i} = DarkBar{i} ( : , xFly); 
 darkB.ficTracInty{i} = DarkBar{i} ( : , yFly); 
end

% Output in degrees of the Panels position

% Pos x=92 is 0 deg (ie facing the fly), I measured this empirically

pxToDeg = 360/97; % There are 97 possible positions (the last one = first one) and this represents 360 deg

for i = 1:size(DarkBar,2)
posToDeg{i} = zeros(1,length(darkB.xPanelPos{i}));

% Convert from xpos to degrees, knowing that xpos 5 = 0 deg
for j=1:length(darkB.xPanelPos{i})
    if isequal(darkB.xPanelPos{i}(j),93) | isequal(darkB.xPanelPos{i}(j),94) | isequal(darkB.xPanelPos{i}(j),95) | isequal(darkB.xPanelPos{i}(j),96) | isequal(darkB.xPanelPos{i}(j),97)
        posToDeg{i}(j) = (darkB.xPanelPos{i}(j)-92)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        posToDeg{i}(j) = (darkB.xPanelPos{i}(j)+4)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end

% Create a time vector in sec
time{i} = linspace(0,(size(DarkBar{i},1)/1000),size(DarkBar{i},1)); %1000 is our sampling rate

end

% How much is the fly moving?

% 1) Assess noise to determine the best criterion for whether the fly is
% moving or not

voltThresh = assessNoise;

for i = 1:size(DarkBar,2)
  voltThreshold{i} = voltThresh;
 [percentMoving(i), moving{i}] = IsFlyWalking(DarkBar{i},voltThresh);

%add a zero before moving start, for the frame 1 to have a "is not moving"
%assigned
 moving{i} = [0,moving{i}];

%change to boolean
 Moving{i} = true(size(moving{i},1),size(moving{i},2));

 for j = 1:length(moving{i})
    if moving{i}(1,j) == 0
        Moving{i}(j) = false;
    else
        Moving{i}(j) = true;
    end
 end

end

% Probability density keeping only those frames when the fly was moving

for i = 1:size(DarkBar,2)
 flyPosToDegMoving{i} = darkB.ficTracAngularPosition{i}(Moving{i}).*36; %if 10 V = 360 deg, then xV = 36x deg

 remapFlyPosToDegMoving{i} = flyPosToDegMoving{i};
 for j = 1:length(remapFlyPosToDegMoving{i})   
    if remapFlyPosToDegMoving{i}(j) > 180
        remapFlyPosToDegMoving{i}(j) = remapFlyPosToDegMoving{i}(j) - 360;  
    end   
 end

% Plot the histogram and probability density
 [countsFlyMoving{i}] = histcounts(remapFlyPosToDegMoving{i},20);
 probabilitiesFlyMoving{i} = countsFlyMoving{i}./sum(countsFlyMoving{i});
 degsFlyMoving{i} = linspace(-180,180,length(countsFlyMoving{i}));

end

% Compute angular velocity, specially for the open-loop grating trials

for i = 1:size(DarkBar,2)
    
 LOWPASS_FILTER_CUTOFF = 100; % Hz
 THRESHOLD_ANGULAR_VELOCITY = 800; % degrees / s  this is the max velocity that can be allowed into analysis
% decode angular velocity and accumulated position
 [angularVelocity{i},filteredFicTracPos{i}] = ficTracSignalDecoding(darkB.ficTracAngularPosition{i}, 1000, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);

end

%Plot
figure,

for i = 1:size(DarkBar,2)

        plot(degsFlyMoving{i},probabilitiesFlyMoving{i})
        xlim([-180 180])
        title('Probability density of the fly heading');
        ylabel('Probability density'); xlabel('Fly heading (deg)');
        hold on
        %add the mean and standard deviation
     
end

DarkBarProbabilities = cell2mat(probabilitiesFlyMoving);
DarkBarProbabilitiesReshaped = reshape(DarkBarProbabilities,[20,length(DarkBarProbabilities)/20]);
meanHeading = mean(DarkBarProbabilitiesReshaped,2);
plot(degsFlyMoving{1},meanHeading,'k','LineWidth',2)

%% For light, closed-loop bars


for i = 1:size(LightBar,2)
 lightB.xPanelVolts{i} =  LightBar{i} (:,xPanels); 
 VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
 maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
 lightB.xPanelPos{i} = round ((lightB.xPanelVolts{i}  * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels

 lightB.yPanelVolts{i} =  LightBar{i} (:, yPanels);
 VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
 maxValY = 1;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2
 lightB.yPanelPos{i} = round ((lightB.yPanelVolts{i}  * maxValY) /VOLTAGE_RANGE);

%FicTrac data
 lightB.ficTracAngularPosition{i} = LightBar{i} ( : , headingFly); 
 lightB.ficTracIntx{i} = LightBar{i} ( : , xFly); 
 lightB.ficTracInty{i} = LightBar{i} ( : , yFly); 
end

% Output in degrees of the Panels position

% Pos x=92 is 0 deg (ie facing the fly), I measured this empirically

pxToDeg = 360/97; % There are 97 possible positions (the last one = first one) and this represents 360 deg

for i = 1:size(LightBar,2)
posToDeg{i} = zeros(1,length(lightB.xPanelPos{i}));

% Convert from xpos to degrees, knowing that xpos 5 = 0 deg
for j=1:length(lightB.xPanelPos{i})
    if isequal(lightB.xPanelPos{i}(j),93) | isequal(lightB.xPanelPos{i}(j),94) | isequal(lightB.xPanelPos{i}(j),95) | isequal(lightB.xPanelPos{i}(j),96) | isequal(lightB.xPanelPos{i}(j),97)
        posToDeg{i}(j) = (lightB.xPanelPos{i}(j)-92)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        posToDeg{i}(j) = (lightB.xPanelPos{i}(j)+4)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end

% Create a time vector in sec
time{i} = linspace(0,(size(LightBar{i},1)/1000),size(LightBar{i},1)); %1000 is our sampling rate

end

% How much is the fly moving?

% 1) Assess noise to determine the best criterion for whether the fly is
% moving or not

voltThresh = assessNoise;

for i = 1:size(LightBar,2)
  voltThreshold{i} = voltThresh;
 [percentMoving(i), moving{i}] = IsFlyWalking(LightBar{i},voltThresh);

%add a zero before moving start, for the frame 1 to have a "is not moving"
%assigned
 moving{i} = [0,moving{i}];

%change to boolean
 Moving{i} = true(size(moving{i},1),size(moving{i},2));

 for j = 1:length(moving{i})
    if moving{i}(1,j) == 0
        Moving{i}(j) = false;
    else
        Moving{i}(j) = true;
    end
 end

end

% Probability density keeping only those frames when the fly was moving

for i = 1:size(LightBar,2)
 flyPosToDegMoving{i} = lightB.ficTracAngularPosition{i}(Moving{i}).*36; %if 10 V = 360 deg, then xV = 36x deg

 remapFlyPosToDegMoving{i} = flyPosToDegMoving{i};
 for j = 1:length(remapFlyPosToDegMoving{i})   
    if remapFlyPosToDegMoving{i}(j) > 180
        remapFlyPosToDegMoving{i}(j) = remapFlyPosToDegMoving{i}(j) - 360;  
    end   
 end

% Plot the histogram and probability density
 [countsFlyMoving{i}] = histcounts(remapFlyPosToDegMoving{i},20);
 probabilitiesFlyMoving{i} = countsFlyMoving{i}./sum(countsFlyMoving{i});
 degsFlyMoving{i} = linspace(-180,180,length(countsFlyMoving{i}));

end

% Compute angular velocity, specially for the open-loop grating trials

for i = 1:size(LightBar,2)
    
 LOWPASS_FILTER_CUTOFF = 100; % Hz
 THRESHOLD_ANGULAR_VELOCITY = 800; % degrees / s  this is the max velocity that can be allowed into analysis
% decode angular velocity and accumulated position
 [angularVelocity{i},filteredFicTracPos{i}] = ficTracSignalDecoding(lightB.ficTracAngularPosition{i}, 1000, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);

end

%Plot

figure,

for i = 1:size(LightBar,2)

        plot(degsFlyMoving{i},probabilitiesFlyMoving{i})
        xlim([-180 180])
        title('Probability density of the fly heading');
        ylabel('Probability density'); xlabel('Fly heading (deg)');
        hold on
        %add the mean and standard deviation
     
end

LightBarProbabilities = cell2mat(probabilitiesFlyMoving);
LightBarProbabilitiesReshaped = reshape(LightBarProbabilities,[20,length(DarkBarProbabilities)/20]);
meanHeading = mean(LightBarProbabilitiesReshaped,2);
plot(degsFlyMoving{1},meanHeading,'k','LineWidth',2)

