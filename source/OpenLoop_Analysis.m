%Code to analyze open-loop tracking behavior

%%% This code analyses the output of the panels and the FicTrac data
close all; clear all;

% prompt the user to select the file to open and load it.
% [file,path] = uigetfile();
% load([path,file]);

openLoop = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment11\*\*\dataOpenLoopBar*.mat');

%load the open loop data and determine velocity
vel = zeros(1,length(openLoop));
for i = 1:length(openLoop)
    openLoopData{i} = load(strcat(openLoop(i).folder,'\',openLoop(i).name));
    if openLoopData{1,i}.velocity == 60
        vel(i) = 60;
    else
        vel(i) = 20;
    end
end

%separate acording to velocity
openLoop20 = cell(size(openLoopData));
openLoop20 = openLoopData(vel==20);
openLoop60 = cell(size(openLoopData));
openLoop60 = openLoopData(vel==60);


%get only the mvt data
for i = 1:length(openLoop20)
    data20{i} = openLoop20{1,i}.rawData;
end

for i = 1:length(openLoop60)
    data60{i} = openLoop60{1,i}.rawData;
end

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;


%% Subset acquisition of x and y pos, as well as FicTrac data

for i = 1:length(data20)
data(i).xPanelVolts =  data20{1,i} (:,xPanels); 
VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
data(i).xPanelPos = round ((data(i).xPanelVolts  * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels

data(i).yPanelVolts =   data20{1,i} (:, yPanels);
VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
maxValY = 1;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2
data(i).yPanelPos = round ((data(i).yPanelVolts  * maxValY) /VOLTAGE_RANGE);

%FicTrac data
data(i).ficTracAngularPosition =  data20{1,i} ( : , headingFly); 
data(i).ficTracIntx = - data20{1,i} ( : , xFly); %the negative sign is necessary under my current conditions (z axis facing up)
data(i).ficTracInty =  data20{1,i}( : , yFly); %I think if I wanted to look into this one, I probably need to invert it too
end

for i = 1:length(data60)
Data60(i).xPanelVolts =  data60{1,i} (:,xPanels); 
Data60(i).xPanelPos = round ((Data60(i).xPanelVolts  * maxValX ) /VOLTAGE_RANGE); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels

Data60(i).yPanelVolts =   data60{1,i} (:, yPanels);
Data60(i).yPanelPos = round ((Data60(i).yPanelVolts  * maxValY) /VOLTAGE_RANGE);

%FicTrac data
Data60(i).ficTracAngularPosition =  data60{1,i} ( : , headingFly); 
Data60(i).ficTracIntx = - data60{1,i} ( : , xFly); %the negative sign is necessary under my current conditions (z axis facing up)
Data60(i).ficTracInty =  data60{1,i}( : , yFly); %I think if I wanted to look into this one, I probably need to invert it too
end

%% Downsample, unwrap and smooth position data, then get velocity and smooth

for i = 1:length(data20)
[smoothed20{i}] = posDataDecoding(data(i),1000);
end

for i = 1:length(data60)
[smoothed60{i}] = posDataDecoding(Data60(i),1000);
end

%Get the side of the turn to divide the plot
for i = 1:length(openLoop20)
   if openLoop20{1,i}.side == 1;
       side(i)=1;
   else
       side(i)=0;
   end    
end

for i = 1:length(openLoop60)
   if openLoop60{1,i}.side == 1;
       side60(i)=1;
   else
       side60(i)=0;
   end    
end



figure,
for i = 1:length(smoothed)
    if side(i) == 1
plot(smoothed{1,i}.angularVel,'k')
    else
plot(smoothed{1,i}.angularVel,'r')   
    end
hold on
end

figure,
for i = 1:length(smoothed60)
    if side60(i) == 1
plot(smoothed60{1,i}.angularVel,'k')
    else
plot(smoothed60{1,i}.angularVel,'r')   
    end
hold on
end

%angularVel
for i= 1:length(side60)
    
end