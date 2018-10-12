function RunSetGoalExp(flyNum,expNum,trialNum,trialLength)
%%This function will run trials of a 2 px light bar on a dark number. They
%%will all be the same length, specified by "trialLength", and they will be
%%as many trials as is specified by "trialNum"

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

% for loop in which for every stripe fixation trial
% the NiDaq will be initialized and the data will be
% acquired as a separated field in the rawData struct?

for i = 1:trialNum
    
% save the starting position for the stimuli
startPosition{i} = [(round(rand*96)+1)];
% set the NiDaq acquisition to whatever length we determined in trialLength
niOI.DurationInSeconds = trialLength;

Panel_com('set_pattern_id', 11); %load the light stripe pattern
Panel_com('set_position',[startPosition{i} 1]);
pause(1); %I think I need this pause in order for the bar not to go right back to the other position with the closed-loop.
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('start');
rawData{i} = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined
Panel_com('stop');
Panel_com('all_off'); %I don't know if I want to keep this

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

save(strcat('dataExpNum',num2str(expNum),'.mat'),'rawData','startPosition','trialLength'); %save as dataExpNum


end