function RunExperiment(trialPattern)
% This function allows you to run an experiment of concatenated trials,
% whose settings you will load from a predetermined configuration

% INPUT
% trialPattern is a vector with the name of the types of trials to run on
% this experiment.

trials = trialPattern;

% set up NiDaq acquisition session
daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:5 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:5
    aI(i).InputType = 'SingleEnded';
end

% Initialize empty cell array for the data
rawData = {};

% for loop in which for every element of the vector of trials,
% the name of the type of stim will be read, the info will be loaded from
% a settings file, the NiDaq will be initialized and the data will be
% acquired as a separated field in the rawData struct?

for i = 1:size(trials,2)
    
% read settings file
run(['settings_',trials{1,i}])
    
Panel_com('start');
rawData{i} = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined
Panel_com('stop');

end

end