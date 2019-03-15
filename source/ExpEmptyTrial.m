% run the NiDaq and no visual stimulation
% this is to record the fly's behavior in the absense of visual stim.


function ExpEmptyTrial(flyNum,expNum)

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment7\';

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
niOI.DurationInSeconds = 200; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:6 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:6
    aI(i).InputType = 'SingleEnded';
end


%acquire data
rawData = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined


if flyNum ==1 %if it's the first fly
   mkdir ([date]) %make a folder with today's date
end


if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment7\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment7\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment7\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end


save(strcat('EmptyTrialExpNum',num2str(expNum),'.mat'),'rawData'); %save


end