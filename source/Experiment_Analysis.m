%% Experiment analysis
% This code analyses multi-trial experiments

close all; clear all;

% prompt the user to select the file to open and load it.
[file,path] = uigetfile();
load([path,file]);

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;

%% Subset acquisition of x and y pos, as well as FicTrac data

for i = 1:size(rawData,2)
 %Panel's data   
 data.xPanelVolts{i} =  rawData{i} (:,xPanels); 
 VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
 maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
 data.xPanelPos{i} = round ((data.xPanelVolts{i} * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels

 data.yPanelVolts{i} =  rawData{i} (:, yPanels);
 VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
 maxValY = 1;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2
 data.yPanelPos{i} = round ((data.yPanelVolts{i} * maxValY) /VOLTAGE_RANGE);

%FicTrac data
 data.ficTracAngularPosition{i} = rawData{i} ( : ,headingFly); 
 data.ficTracIntx{i} = rawData{i} ( : , xFly); 
 data.ficTracInty{i} = rawData{i} ( : , yFly); 
end


%% Output in degrees of the Panels position

% Pos x=92 is 0 deg (ie facing the fly), I measured this empirically

pxToDeg = 360/97; % There are 97 possible positions (the last one = first one) and this represents 360 deg

for i = 1:size(rawData,2)
    
posToDeg{i} = zeros(1,length(data.xPanelPos{i}));

% Convert from xpos to degrees, knowing that xpos 92 = 0 deg
for j=1:length(data.xPanelPos{i})
    if isequal(data.xPanelPos{i}(j),93) | isequal(data.xPanelPos{i}(j),94) | isequal(data.xPanelPos{i}(j),95) | isequal(data.xPanelPos{i}(j),96) | isequal(data.xPanelPos{i}(j),97)
        posToDeg{i}(j) = (data.xPanelPos{i}(j)-92)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        posToDeg{i}(j) = (data.xPanelPos{i}(j)+4)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end

% Create a time vector in sec
time{i} = linspace(0,(size(rawData{i},1)/10000),size(rawData{i},1)); %10000 is our sampling rate

end


%% Probability density function of the stimulus position

% Remapping the positions to span -180 to 180 deg

remapPosToDeg = posToDeg;

for p = 1:size(remapPosToDeg,2)
    
 for i = 1:length(remapPosToDeg{p})   
    if remapPosToDeg{p}(i) > 180
       remapPosToDeg{p}(i) = remapPosToDeg{p}(i) - 360;  
    end   
 end
 
 % Calculate the probability densities
[counts{p}] = histcounts(remapPosToDeg{p},20);
probabilities{p} = counts{p}./sum(counts{p});
degs{p} = linspace(-180,180,length(counts{p}));

end

% remap starting position for the bar trials
remapStartPos = cell(1,size(rawData,2));
for p = 1:size(rawData,2)
    %transform and remap stim starting pos if closed-loop bar
 if isequal(trials{p}, 'ClosedLoopBar') | isequal(trials{p}, 'DarkClosedLoopBar')
    remapStartPos{p} = 0;
     if isequal(startPosition{p}(1),93) | isequal(startPosition{p}(1),94) | isequal(startPosition{p}(1),95) | isequal(startPosition{p}(1),96) | isequal(startPosition{p}(1),97)
        startingPos{p} = (startPosition{p}(1)-92)*pxToDeg; % Correct the offset and multiply by factor to get deg
     else
        startingPos{p} = (startPosition{p}(1)+4)*pxToDeg; % Correct the offset and multiply by factor to get deg
     end
    
     if startingPos{p} > 180
       remapStartPos{p} = startingPos{p} - 360;
     else
       remapStartPos{p} = startingPos{p};
     end 
 end
 
end


%%  How much is the fly moving?

% 1) Assess noise to determine the best criterion for whether the fly is
% moving or not
voltThresh = assessNoise;

% 2) Use the voltage threshold obtained by the empty trial to determne the
% frames in which the fly is moving
for i = 1:size(rawData,2)
  voltThreshold{i} = voltThresh;
 [percentMoving(i), moving{i}] = IsFlyWalking(rawData{i},voltThresh);

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


%% Probability density keeping only those frames when the fly was moving

for i = 1:size(rawData,2)
 flyPosToDegMoving{i} = data.ficTracAngularPosition{i}(Moving{i}).*36; %if 10 V = 360 deg, then xV = 36x deg

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


%% Compute angular velocity, specially for the open-loop grating trials

for i = 1:size(rawData,2)
    
 LOWPASS_FILTER_CUTOFF = 50; % Hz
 THRESHOLD_ANGULAR_VELOCITY = 200; % degrees / s  this is the max velocity that can be allowed into analysis
% decode angular velocity and accumulated position
 [angularVelocity{i},filteredFicTracPos{i}] = ficTracSignalDecoding(data.ficTracAngularPosition{i}, 10000, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);
 [forwardVelocity{i},filteredForwardPos{i}] = ficTracSignalDecoding(data.ficTracInty{i}, 10000, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);
end


%% Plot everything

for i = 1:size(rawData,2)
    figure,
    subplot(1,2,1)
    plot(time{i},posToDeg{i})
    ylim([0 365]);
    ylabel('Position of the stimulus (deg)');
    xlabel('Time (s)');
    title('Position of the stimulus in deg as a function of time', 'Interpreter', 'none');
 
% For the closed-loop bars, plot stim position and fly heading
    if isequal(trials{i}, 'ClosedLoopBar')
    subplot(1,2,2)
    plot(degsFlyMoving{i},probabilitiesFlyMoving{i},'k')
    xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving{i})+0.05]);
    title('Probability density of the fly heading');
    ylabel('Probability density'); xlabel('Fly heading (deg)');
    hold on
    line([remapStartPos{i} remapStartPos{i}],[0 max(probabilitiesFlyMoving{i})+0.05],'Color',[0 0 1])

    elseif isequal(trials{i}, 'DarkClosedLoopBar')
    subplot(1,2,2)
    plot(degsFlyMoving{i},probabilitiesFlyMoving{i},'b')
    xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving{i})+0.05]);
    title('Probability density of the fly heading');
    ylabel('Probability density'); xlabel('Fly heading (deg)');
    hold on
    line([remapStartPos{i} remapStartPos{i}],[0 max(probabilitiesFlyMoving{i})+0.05],'Color',[0 0 0])
        
        
% For the open-loop gratings, plot stim position and fly's angular velocity
    elseif isequal(trials{i}, 'OpenLoopGratingRight')
    subplot(1,2,2)
    plot(time{i},angularVelocity{i})
    ylabel('Angular velocity of the fly');
    xlabel('Time (s)');
    title('Angular velocity of the fly for Right Garting trials');
    
    else
    subplot(1,2,2)
    plot(time{i},angularVelocity{i})
    ylabel('Angular velocity of the fly');
    xlabel('Time (s)');
    title('Angular velocity of the fly for Left Grating trials');
        
    
    end
    
end

%% Plot by trial type

% For dark-closed loop bars

figure,

for i = 1:size(rawData,2)
    
    if isequal(trials{i},'DarkClosedLoopBar')
        plot(degsFlyMoving{i},probabilitiesFlyMoving{i})
        xlim([-180 180])
        title({'Probability density of the fly heading';'Dark closed-loop bar'});
        ylabel('Probability density'); xlabel('Fly heading (deg)');
        hold on
        %add the mean and standard deviation
    end
     
end

DarkBarTrials = cellfun(@(x) isequal(x,'DarkClosedLoopBar'), trials);
DarkBarProbabilities = cell2mat(probabilitiesFlyMoving(DarkBarTrials));
DarkBarProbabilitiesReshaped = reshape(DarkBarProbabilities,[20,length(DarkBarProbabilities)/20]);
meanHeading = mean(DarkBarProbabilitiesReshaped,2);
stdHeading = std(DarkBarProbabilitiesReshaped,[],2);
steHeading = stdHeading/(sqrt(sum(DarkBarTrials)));
[hl,hp]= boundedline(degsFlyMoving{1},meanHeading, steHeading,'k','alpha')
ho = outlinebounds(hl,hp);
set(ho, 'linestyle', ':', 'color', 'k');
saveas(gcf,strcat(path,'AllDarkBarTrials_PDHeading_ExpNum', file(end-4), '.png'));


% For light closed loop bars

figure,

for i = 1:size(rawData,2)
    
    if isequal(trials{i},'ClosedLoopBar')
        plot(degsFlyMoving{i},probabilitiesFlyMoving{i})
        xlim([-180 180])
        title({'Probability density of the fly heading';'Light closed-loop bar'});
        ylabel('Probability density'); xlabel('Fly heading (deg)');
        hold on
        %add the mean and standard deviation
    end
     
end

LightBarTrials = cellfun(@(x) isequal(x,'ClosedLoopBar'), trials);
LightBarProbabilities = cell2mat(probabilitiesFlyMoving(LightBarTrials));
LightBarProbabilitiesReshaped = reshape(LightBarProbabilities,[20,length(LightBarProbabilities)/20]);
meanHeading = mean(LightBarProbabilitiesReshaped,2);
stdHeading = std(LightBarProbabilitiesReshaped,[],2);
steHeading = stdHeading/(sqrt(sum(LightBarTrials)));
[hl,hp]= boundedline(degsFlyMoving{1},meanHeading, steHeading,'k','alpha')
ho = outlinebounds(hl,hp);
set(ho, 'linestyle', ':', 'color', 'k');
saveas(gcf,strcat(path,'AllLightBarTrials_PDHeading_ExpNum', file(end-4), '.png'));


% For Right side open-loop gratings

figure,

for i = 1:size(rawData,2)
    
    if isequal(trials{i},'OpenLoopGratingRight')
        plot(time{i},angularVelocity{i})
        ylabel('Angular velocity of the fly');
        xlabel('Time (s)');
        title({'Angular velocity of the fly';'Clockwise open-loop gratings'});
        hold on
    end
     
end

OpenLoopGratingTrials = cellfun(@(x) isequal(x,'OpenLoopGratingRight'), trials);
OpenLoopGratingVelocities = cell2mat(angularVelocity(OpenLoopGratingTrials));
meanVelocity = mean(OpenLoopGratingVelocities,2);
stdVelocity = std(OpenLoopGratingVelocities,[],2);
steVelocity = stdVelocity/(sqrt(sum(OpenLoopGratingTrials)));
[hl,hp]= boundedline(time{find(OpenLoopGratingTrials==1,1)},meanVelocity, steVelocity,'k','alpha')
ho = outlinebounds(hl,hp);
set(ho, 'linestyle', ':', 'color', 'k');
saveas(gcf,strcat(path,'AllRightGratingTrials_PDHeading_ExpNum', file(end-4), '.png'));

% For Left side open-loop gratings

figure,

for i = 1:size(rawData,2)
    
    if isequal(trials{i},'OpenLoopGratingLeft')
        plot(time{i},angularVelocity{i})
        ylabel('Angular velocity of the fly');
        xlabel('Time (s)');
        title({'Angular velocity of the fly';'Counter clock wise open-loop gratings'});
        hold on
    end
     
end

OpenLoopGratingTrialsLeft = cellfun(@(x) isequal(x,'OpenLoopGratingLeft'), trials);
OpenLoopGratingVelocitiesLeft = cell2mat(angularVelocity(OpenLoopGratingTrialsLeft));
meanVelocityLeft = mean(OpenLoopGratingVelocitiesLeft,2);
stdVelocityLeft = std(OpenLoopGratingVelocitiesLeft,[],2);
steVelocityLeft = stdVelocityLeft/(sqrt(sum(OpenLoopGratingTrialsLeft)));
[hl,hp]= boundedline(time{find(OpenLoopGratingTrialsLeft==1,1)},meanVelocityLeft, steVelocityLeft,'k','alpha')
ho = outlinebounds(hl,hp);
set(ho, 'linestyle', ':', 'color', 'k');
saveas(gcf,strcat(path,'AllLeftBarTrials_PDHeading_ExpNum', file(end-4), '.png'));