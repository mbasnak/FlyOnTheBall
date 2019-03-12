function [daq_data] = ExpYoked_LowError(flyNum,expNum)

% This function runs an experiment in which the experience of a low error
% block from a master fly is presented to a yolked fly.

%INPUT
%flyNum : which number of fly are you running
%expNum : which experiment number are you running

%OUTPUT
%daq_data : matrix that returns the voltage data acquired in the NiDaq

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
aI = niOI.addAnalogInputChannel( devID , 1:6 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:6
    aI(i).InputType = 'SingleEnded';
end
niOI.DurationInSeconds = 2420;

fid = fopen('log.dat','w+'); %this opens a csv file named "log",creating it for writing (and overwriting existing filed) with the flag "w+"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground(); %start acquiring

startPos = 2; %to match the starting position of the Y pattern.
jumpFunction = randperm(4,1)+19 %get a random number from 1 to 4 to determine the pos function

if jumpFunction == 20
    jumps = [-90,-90,90,-90,90,90,90,90,-90,-90,90,-90];
elseif jumpFunction == 21
    jumps = [90,-90,-90,90,-90,-90,90,90,-90,90,-90,90];
elseif jumpFunction == 22
    jumps = [-90,90,-90,-90,90,90,-90,-90,90,90,-90,90];
elseif jumpFunction == 23
    jumps = [90,-90,-90,90,-90,-90,90,90,-90,90,-90,90];
end
%not sure about these...

%%%%%% Run the panels %%%%%%
Panel_com('set_pattern_id', 1); %set the bar
Panel_com('set_mode', [4 0]); 
pause(0.03)
Panel_com('send_gain_bias', [0 0 0 0]);
pause(0.03)
Panel_com('set_position',[startPos 1]);
pause(0.03)
Panel_com('set_funcx_freq', 25);
pause(0.03)
Panel_com('set_posfunc_id',[1 jumpFunction]); %set it to be yoked, to one of the 4 functions made from master flies.
pause(0.03)
Panel_com('set_AO',[3 32767]);
Panel_com('start');
pause(2404) %record for the time it takes to span the number of trials requested
Panel_com('stop');
Panel_com('set_AO',[3 0]);
Panel_com('all_off');

% add a signal to make sure that the panels have stopped running

niOI.IsDone % will report 0
niOI.wait() % wait and keep acquiring until the acquisition time is over
niOI.IsDone % will report 1

delete(lh) % delete the listener handle

[daq_data] = loadFromLogFile('log.dat',6); %load the data just saved to the dat file, using Stephen's function

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment6\' %move to our data directory
    
% if flyNum ==1 %if it's the first fly
%    mkdir ([date]) %make a folder with today's date
% end

if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment6\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment6\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment6\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end

save(strcat('Yoked_LowErrorBlockExp',num2str(expNum),'.mat'),'daq_data','startPos','jumps'); %save daq data and starting positions as a .mat file


end