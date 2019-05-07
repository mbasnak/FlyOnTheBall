% run the NiDaq and an open-loop bar at different speeds

function ExpOpenLoopBar(flyNum,expNum,folder)

cd(strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\'));

%define the velocity of the experiment
velocity = round(rand)*40 + 20; %I set it to either 20 or 60 deg/s

%set the acquisition time acording to the velocity.
%if bar moves at 20 deg/s, then 18 sec. Else,6 sec
if velocity == 20
    time = 18;
else
    time = 6;
end

%define the side of the turn
side = round(rand);
%side ==1 will be clockwise

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
niOI.DurationInSeconds = time + 1; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:6 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:6
    aI(i).InputType = 'SingleEnded';
end

%define the mvtFunction acording to the velocity and turn side:
if velocity == 20 && side == 1
    mvtFunct = 50;
elseif velocity == 20 && side == 0
    mvtFunct = 51;
elseif velocity == 60 && side == 1
    mvtFunct = 54;
else
    mvtFunct = 55;
end

%set a rendom starting point for the stim
Panel_com('set_pattern_id', 18); %load a 4 px dark stripe
pause(0.01)
Panel_com('set_position',[20,2]);
Panel_com('set_mode', [4 1]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_funcx_freq', 50);
Panel_com('set_posfunc_id',[1 mvtFunct]);
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


save(strcat('dataOpenLoopBar',num2str(expNum),'.mat'),'rawData','velocity','side'); %save as 


end