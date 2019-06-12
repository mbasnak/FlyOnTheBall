function RunMultitrialExperiment(flyNum,expNum,TrialNum,folder)
% This function allows you to run an experiment of different concatenated trials,
% whose settings you will load from a predetermined configuration, and have
% them randomly presented (but predefine in what proportion you're
% presenting each trial type)

% INPUT
% flyNum
% expNum
% TrialNum: how many times I am repeating my basic presentation of 14
% trials to the fly.
% folder: in which folder are you saving the data


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

trialTypes = {'lightClosedLoopBar','darkClosedLoopBar','fastClockwiseOpenLoop','fastCounterclockwiseOpenLoop','fastClockwiseOpenLoop','fastCounterclockwiseOpenLoop','fastClockwiseOpenLoop','fastCounterclockwiseOpenLoop','clockwiseFourier','clockwiseFourier','clockwiseFourier','counterclockwiseFourier','counterclockwiseFourier','counterclockwiseFourier'};
trials = repmat(trialTypes,1,TrialNum);
trials = trials(randperm(length(trials)));


% for loop in which for every element of the vector of trials,
% the name of the type of stim will be read, the info will be loaded from
% a settings file, the NiDaq will be initialized and the data will be
% acquired as a separated field in the data struct
for i = 1:length(trials)   
% read settings file and run the panels commands
    run(['settings_',trials{1,i}]) 
    startPosition{i} = startPos;
    Panel_com('start');
    rawData{i} = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined
    Panel_com('stop');
    Panel_com('all_off')
    pause(3);

end

Panel_com('all_off'); %turn off the panels

cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder)]);
    
if flyNum ==1 %if it's the first fly
   mkdir ([date]) %make a folder with today's date
end

if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end

% %save each trial type data separately
% for i = 1:length(trials)
%     save(strcat('data',trials{i},num2str(i),'.mat'),rawData{i}); %save each trial separately
% end

save(strcat('dataExpNum',num2str(expNum),'.mat'),'rawData','startPosition','trials'); %save as dataExpNum

% 
% %save the general info
% save('generalInfo.mat','trials','startpos');


end