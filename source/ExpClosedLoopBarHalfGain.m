% run the NiDaq and a closed-loop bar

function ExpClosedLoopBarHalfGain(flyNum,expNum,time)

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3';

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
niOI.DurationInSeconds = time; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:6 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:6
    aI(i).InputType = 'SingleEnded';
end

%set a rendom starting point for the stim
Panel_com('set_pattern_id', 16); %load the light stripe pattern
pause(0.01)
Panel_com('send_gain_bias', [10 0 0 0]);
pause(0.01)
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('start');


%acquire data
rawData = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined

pause(2);
Panel_com('stop');
Panel_com('all_off');
Panel_com('set_AO',[3 0]);


if flyNum ==1 %if it's the first fly
   mkdir ([date]) %make a folder with today's date
end


if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment3\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end


save(strcat('dataCloseLoopBarHalfGain',num2str(expNum),'.mat'),'rawData'); %save as 


end