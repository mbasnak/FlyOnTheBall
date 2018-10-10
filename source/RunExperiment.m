function RunExperiment(trialPattern,flyNum,expNum)
% This function allows you to run an experiment of concatenated trials,
% whose settings you will load from a predetermined configuration

% INPUT
% trialPattern is a vector with the name of the types of trials to run on
% this experiment.

% shuffle the trials position
trials = trialPattern;%(randperm(numel(trialPattern)));

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

% Initialize empty cell arrays for the data
rawData = {};
startPosition = {};

% for loop in which for every element of the vector of trials,
% the name of the type of stim will be read, the info will be loaded from
% a settings file, the NiDaq will be initialized and the data will be
% acquired as a separated field in the rawData struct?

for i = 1:size(trials,2)
    
% read settings file and run the panels commands
run(['settings_',trials{1,i}])

% save the starting position for the stimuli
startPosition{i} = startPos;
    
Panel_com('start');
rawData{i} = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined
Panel_com('stop');

end

Panel_com('all_off'); %turn off the panels

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data';
    
if flyNum ==1 %if it's the first fly
   mkdir ([date]) %make a folder with today's date
end

if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end

save(strcat('dataExpNum',num2str(expNum),'.mat'),'rawData','startPosition','trials'); %save as dataExpNum

end