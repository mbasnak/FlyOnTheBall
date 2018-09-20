%% Code to acquire the pannel x and y information to later coordinate with fictrac output
clear all;
daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
niOI.DurationInSeconds = 15; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:5 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:5
    aI(i).InputType = 'SingleEnded';
end

%acquire data
rawData = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined


%define channels ID
headingFly = 1;
xFly = 2;
yFly = 3;
xPanels = 4;
yPanels = 5;

%% set acquistion for x and y data from the panels

data.xPanelVolts =  rawData (:,xPanels); 
   % Decode Xpos from voltage reading
% data.xPanelPos = processPanelDataX ( data.xPanelVolts , stimulus.panelParams ); %I can't find where stimulus.panelParams is defined.
VOLTAGE_RANGE = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
data.xPanelPos = round ((data.xPanelVolts  * maxValX ) /VOLTAGE_RANGE);

data.yPanelVolts =  rawData (:, yPanels);
   % Decode Ypos from voltage reading
% data.yPanelPos = processPanelDataY ( data.yPanelVolts , stimulus.panelParams );
VOLTAGE_RANGE = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
maxValY = 1;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2
data.yPanelPos = round ((data.yPanelVolts  * maxValY) /VOLTAGE_RANGE);

%set acquisition for fictrac data
data.ficTracAngularPosition = rawData ( : , headingFly); 
data.ficTracIntx = rawData ( : , xFly); 
data.ficTracInty = rawData ( : , yFly); 
