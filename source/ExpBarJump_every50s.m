function [daq_data] = ExpBarJump_every50s(flyNum,expNum,TrialNum)

% This function runs an experiment in which a bar is presented in closed-loop and
% every 50 s, it jumps to a random position

%INPUT
%flyNum : which number of fly are you running
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
niOI.DurationInSeconds = TrialNum*50 + 5; %set the duration to the total panel display time + 5 extra seconds

fid = fopen('log.dat','w+'); %this opens a csv file named "log",creating it for writing (and overwriting existing filed) with the flag "w+"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground(); %start acquiring

startPos = round(rand*96)+1; % generate a random starting position

%%%%%% Run the panels %%%%%%
Panel_com('set_pattern_id', 12); %set the bar
Panel_com('set_mode', [3 4]); %set it to closed-loop mode in the x dimension and to be controlled by a function in the y dimension 
Panel_com('set_position',[startPos 1]);
Panel_com('set_posfunc_id',[2 4]); %set it to jump every 10 sec.
Panel_com('set_AO',[3 32767]);
Panel_com('start');
pause(TrialNum*50) %record for the time it takes to span the number of trials requested
Panel_com('stop');
Panel_com('set_AO',[3 0]);
Panel_com('all_off');

% add a signal to make sure that the panels have stopped running

niOI.IsDone % will report 0
niOI.wait() % wait and keep acquiring until the acquisition time is over
niOI.IsDone % will report 1

delete(lh) % delete the listener handle

[daq_data] = loadFromLogFile('log.dat',6); %load the data just saved to the dat file, using Stephen's function

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\' %move to our data directory
    
if flyNum ==1 %if it's the first fly
   mkdir ([date]) %make a folder with today's date
end
%For some reason this isn't working

if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end

save(strcat('dataExpNum',num2str(expNum),'.mat'),'daq_data','startPos','TrialNum'); %save daq data and starting positions as a .mat file


end