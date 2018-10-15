function [daq_data] = acquireData(flyNum,expNum,TrialNum,TrialLength)

% This function acquires data from the NiDaq in the background, allowing
% other matlab code to execute in the meantime.
% It is for when you want to change the stimulation and record the signals
% during those transitions

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
aI = niOI.addAnalogInputChannel( devID , 1:5 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:5
    aI(i).InputType = 'SingleEnded';
end
niOI.DurationInSeconds = (TrialNum*TrialLength)+10;

fid = fopen('log.dat','w+'); %this opens a csv file named "log", 
%creating it for writing with the flag "w"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground();

%%%%%% Run the panels %%%%%%
StartPosx = (round(rand(TrialNum,1)*96)+1);

Panel_com('set_pattern_id', 11);
Panel_com('set_mode', [3 0]);

for i = 1:TrialNum
Panel_com('set_position',[StartPosx(i) 1]);
Panel_com('start');
pause(TrialLength)
Panel_com('stop');
end

Panel_com('all_off');
% add a signal to make sure that the panels have stopped running

niOI.IsDone % will report 0
niOI.wait() % rather than while
niOI.IsDone % will report 

delete(lh)

[daq_data] = loadFromLogFile('log.dat',5);

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

save(strcat('dataExpNum',num2str(expNum),'.mat'),'daq_data','StartPosx'); %save as dataExpNum

end
