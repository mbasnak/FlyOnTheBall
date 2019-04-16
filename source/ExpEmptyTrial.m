% run the NiDaq and no visual stimulation
% this is to record the fly's behavior in the absense of visual stim.


function ExpEmptyTrial(flyNum,expNum,folder)

%flyNum is the number of the fly that you are running in the day
%expNum refers to the experiment number that this represents of the ones
%this fly received thus far
%folder refers to the Experiment folder number in which this file should
%go. If it belongs in 'Experiment10', then folder will be 10.

cd(strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\'));

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


if flyNum == 1 && expNum == 1 %if it's the first fly and the first experiment
   mkdir ([date]) %make a folder with today's date
end


if expNum == 1 %if it's the first experiment for this fly
   cd(strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\',date)); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd(strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\',date,'\flyNum',num2str(flyNum)));
   getFlyInfo() %get fly's details
else
   cd(strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\',date,'\flyNum',num2str(flyNum))); %otherwise move to this fly's folder
end


save(strcat('EmptyTrialExpNum',num2str(expNum),'.mat'),'rawData'); %save


end