% run the NiDaq and an openLoop grating

function ExpOpenLoopGrating(flyNum,expNum,time)

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data';

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
niOI.DurationInSeconds = time; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:5 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:5
    aI(i).InputType = 'SingleEnded';
end

Panel_com('set_pattern_id', 4); %the pattern 2 are the gratings
Panel_com('set_mode', [4 4]); %set the mode to be controlled by an outside function
Panel_com('set_position',[5 1]); %we should run this if we want to stimuli to start centered (doesn't make a lot of sense for the gratings though)
Panel_com('set_posfunc_id',[1 2]);
Panel_com('set_posfunc_id',[2 2]);
Panel_com('start');

%acquire data
rawData = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined

pause(2);
Panel_com('stop');
Panel_com('all_off');

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


typeOfStim = 'open_loop_grating';


save(strcat('dataExpNum',num2str(expNum),'.mat'),'rawData','typeOfStim'); %save as dataExpNum

end