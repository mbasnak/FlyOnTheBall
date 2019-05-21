function [daq_data] = ExpHighErrorBlock(flyNum,expNum)

% This function runs an experiment in which a bar is presented in closed-loop and
% every 200 s, it jumps to a position at either 90 or -90 deg from the
% current one.

%INPUT
%flyNum : which number of  fly are you running
%expNum : which experiment number are you running
%TrialNum : how many closed-loop trials are you running

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
niOI.DurationInSeconds = 1000; 

fid = fopen('log.dat','w+'); %this opens a csv file named "log",creating it for writing (and overwriting existing filed) with the flag "w+"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground(); %start acquiring

startPos = 2; %to match the starting position of the Y pattern.

jumpFunction = randperm(5,1)+28 %get a random number from 1 to 5 to determine the pos function

if jumpFunction == 29
    jumps = [90,90,90,90,-90,-90,-90,-90,-90,-90,90,90,-90,90,-90,-90,90,90,-90,-90,90,-90,90,90,90,90,90,90,-90,-90,-90,-90,-90,-90,90,90,-90,90,-90,-90,90,90,-90,-90,90,-90,90,90];
elseif jumpFunction == 30
    jumps = [-90,-90,-90,90,90,90,90,-90,-90,90,90,-90,-90,90,-90,90,-90,90,-90,-90,90,-90,90,90,-90,-90,-90,90,90,90,90,-90,-90,90,90,-90,-90,90,-90,90,-90,90,-90,-90,90,-90,90,90];
elseif jumpFunction == 31
    jumps = [-90,-90,90,-90,-90,90,90,-90,90,-90,-90,90,-90,90,90,90,-90,90,90,-90,90,-90,90,-90,-90,-90,90,-90,-90,90,90,-90,90,-90,-90,90,-90,90,90,90,-90,90,90,-90,90,-90,90,-90];
elseif jumpFunction == 32
    jumps = [90,90,-90,90,90,-90,-90,-90,90,-90,90,-90,-90,90,-90,90,90,90,-90,-90,-90,-90,90,90,90,90,-90,90,90,-90,-90,-90,90,-90,90,-90,-90,90,-90,90,90,90,-90,-90,-90,-90,90,90];
elseif jumpFunction == 33
    jumps = [-90,90,-90,-90,90,90,-90,90,90,-90,-90,90,90,-90,-90,-90,-90,90,-90,90,90,90,90,-90,-90,90,-90,-90,90,90,-90,90,90,-90,-90,90,90,-90,-90,-90,-90,90,-90,90,90,90,90,-90];
end

%%%%%% Run the panels %%%%%%
Panel_com('set_pattern_id', 1); %set the bar
Panel_com('set_mode', [3 4]); %set it to closed-loop mode in the x dimension and to be controlled by a function in the y dimension 
pause(0.03)
Panel_com('set_position',[startPos 1]);
pause(0.03)
Panel_com('set_funcy_freq', 5);
pause(0.03)
Panel_com('set_posfunc_id',[2 jumpFunction]); %set it to jump every 20 sec, to one of the 5 jumping functions created.
pause(0.03)
Panel_com('set_AO',[3 32767]);
Panel_com('start');
pause(984) %record for the time it takes to span the number of trials requested
Panel_com('stop');
Panel_com('set_AO',[3 0]);
Panel_com('all_off');

% add a signal to make sure that the panels have stopped running

niOI.IsDone % will report 0
niOI.wait() % wait and keep acquiring until the acquisition time is over
niOI.IsDone; % will report 1

delete(lh); % delete the listener handle

[daq_data] = loadFromLogFile('log.dat',6); %load the data just saved to the dat file, using Stephen's function

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\' %move to our data directory
    
% if flyNum ==1 %if it's the first fly
%    mkdir ([date]) %make a folder with today's date
% end

if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment13\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end

save(strcat('HighErrorBlockExp',num2str(expNum),'.mat'),'daq_data','startPos','jumps', 'jumpFunction'); %save daq data and starting positions as a .mat file


end