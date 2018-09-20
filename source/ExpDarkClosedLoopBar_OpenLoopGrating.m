% run the NiDaq and a dark closed-loop bar followed by a short open-loop
% grating

function ExpDarkClosedLoopBar_OpenLoopGrating(flyNum,expNum)

cd 'Z:\Wilson Lab\Mel\codes\behavior';

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
niOI.DurationInSeconds = 20; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:5 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:5
    aI(i).InputType = 'SingleEnded';
end

%set a rendom starting point for the stim
startPos = [round(rand*97) 1];

Panel_com('set_pattern_id', 3); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_position',startPos); %we can also comment this out, or start at [5 1]
Panel_com('start');
%acquire data
rawData.trial1 = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined
Panel_com('stop');

niOI.DurationInSeconds = 2; %set duration in seconds
Panel_com('set_pattern_id', 2); %the pattern 2 are the gratings
Panel_com('set_mode', [4 4]); %set the mode to be controlled by an outside function
Panel_com('set_position',[5 1]); %we should run this if we want to stimuli to start centered (doesn't make a lot of sense for the gratings though)
Panel_com('set_posfunc_id',[1 3]);
Panel_com('set_posfunc_id',[2 3]);
Panel_com('start');
rawData.trial2 = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined
pause(1);
Panel_com('stop');
Panel_com('all_off');


if flyNum ==1 %if it's the first fly
   mkdir ([date]) %make a folder with today's date
end


if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\codes\behavior\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\codes\behavior\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\codes\behavior\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end


typeOfStim = 'dark_closed_loop_bar/open_loop_grating';

save(strcat('dataExpNum',num2str(expNum),'.mat'),'rawData','typeOfStim','startPos'); %save as dataExpNum


end