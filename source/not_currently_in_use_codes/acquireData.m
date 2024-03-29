function [daq_data] = acquireData(flyNum,expNum,TrialNum,TrialLength)

% This function acquires data from the NiDaq in the background, allowing
% other matlab code to execute in the meantime.
% It is for when you want to change the stimulation and record the signals
% during those transitions

%INPUT
%flyNum : which number of fly are you running
%expNum : which experiment number are you running
%TrialNum : how many closed-loop trials are you running
%TrialLength : how long is each trial

%OUTPUT
%daq_data : matrix that returns the voltage data acquired in the NiDaq


daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
aI = niOI.addAnalogInputChannel( devID , 1:5 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:5
    aI(i).InputType = 'SingleEnded';
end
niOI.DurationInSeconds = (TrialNum*TrialLength)+10; %set the duration to the total panel display time + 10 extra seconds

fid = fopen('log.dat','w+'); %this opens a csv file named "log",creating it for writing (and overwriting existing filed) with the flag "w+"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground(); %start acquiring


%%%%%% Run the panels %%%%%%
StartPosx = (round(rand(TrialNum,1)*96)+1); %generate a vector of as many random starting positions as trials we are running

Panel_com('set_pattern_id', 11); %set the bar
Panel_com('set_mode', [3 0]); %set it to closed-loop mode

for i = 1:TrialNum %for each trial     
Panel_com('set_position',[StartPosx(i) 1]); %set the position to one of the random ones from the vector
pause(1)
Panel_com('start');
pause(TrialLength) %go on for as long as your trials last
Panel_com('stop');
end

Panel_com('all_off');

% add a signal to make sure that the panels have stopped running

niOI.IsDone % will report 0
niOI.wait() % wait and keep acquiring until the acquisition time is over
niOI.IsDone % will report 1

delete(lh) % delete the listener handle

[daq_data] = loadFromLogFile('log.dat',5); %load the data just saved to the dat file, using Stephen's function

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

save(strcat('dataExpNum',num2str(expNum),'.mat'),'daq_data','StartPosx'); %save daq data and starting positions as a .mat file

end
