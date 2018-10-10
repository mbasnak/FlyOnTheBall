% run the NiDaq and a dark closed-loop bar

function ExpDarkClosedLoopBar(flyNum,expNum,time)

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

%set a rendom starting point for the stim
startPos = [round(rand*97) 1];

Panel_com('set_pattern_id', 9); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_position',startPos); %we can also comment this out, or start at [5 1]
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


typeOfStim = 'dark_closed_loop_bar';

save(strcat('dataExpNum',num2str(expNum),'.mat'),'rawData','typeOfStim','startPos'); %save as dataExpNum


end