% run the NiDaq and a closed-loop bar

function ExpVirtualHallway(flyNum,expNum,time,folder)

cd(strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\'));

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 4000;% set sample rate
niOI.DurationInSeconds = time; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:7 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:7
    aI(i).InputType = 'SingleEnded';
end

%set a rendom starting point for the stim
pattern = 34;

Panel_com('set_pattern_id', pattern); %load the light stripe pattern
Panel_com('set_mode', [0 3]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('start');


%acquire data
rawData = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined

pause(2);
Panel_com('stop');
Panel_com('all_off');
Panel_com('set_AO',[3 0]);


if flyNum == 1 && expNum == 1%if it's the first fly and the first experiment
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


save(strcat('dataVirtualHallway',num2str(expNum),'.mat'),'rawData','pattern'); %save as 


end