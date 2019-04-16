function [daq_data] = ExpOptomotor(flyNum,expNum,TrialNum)

% This function runs an experiment that alternates short stripe fixation
% trials with short optomotor response trials

%INPUT
%flyNum : which number of  fly are you running
%expNum : which experiment number are you running
%TrialNum : how many stripe-fixation+optomotor response trials are you
%running

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
niOI.DurationInSeconds = 14*TrialNum; 

fid = fopen('log.dat','w+'); %this opens a csv file named "log",creating it for writing (and overwriting existing filed) with the flag "w+"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground(); %start acquiring

startPos = 2; %to match the starting position of the Y pattern.

%%%%%% Run the panels in a loop %%%%%%
count = 1;
jumpFunction = zeros(1,TrialNum); %I'll store the jumpFunction for the optomotor response to know the direction

for count = 1:TrialNum
    %stripe-fixation trial
    Panel_com('set_pattern_id', 1); %set the bar
    Panel_com('set_mode', [3 1]); %set it to closed-loop mode in the x dimension 
    pause(0.03)
    Panel_com('set_position',[startPos 1]);
    pause(0.03)
    Panel_com('set_AO',[3 32767]);
    Panel_com('start');
    pause(8) %record for the time it takes to span the number of trials requested
    Panel_com('stop');
    Panel_com('set_AO',[3 0]);
    Panel_com('all_off');
    pause(1)
    
    %optomotor response trial
    jumpFunction(count) = round(rand+50); %toggle between clockwise (function 50) and counterclockwise (function 51)
    Panel_com('set_pattern_id', 15); %set the grating
    Panel_com('set_mode', [4 1]); %set it to function-driven in the x dimension 
    pause(0.03)
    Panel_com('set_position',[1 1]);
    pause(0.03)
    Panel_com('set_funcx_freq', 50);
    pause(0.03)
    Panel_com('set_posfunc_id',[1 jumpFunction(count)]);
    Panel_com('set_AO',[3 32767]);
    Panel_com('start');
    pause(3) %record for the time it takes to span the number of trials requested
    Panel_com('stop');
    Panel_com('set_AO',[3 0]);
    Panel_com('all_off');
    pause(1) 

count = count+1;
end

% add a signal to make sure that the panels have stopped running
niOI.IsDone % will report 0
niOI.wait() % wait and keep acquiring until the acquisition time is over
niOI.IsDone; % will report 1

delete(lh); % delete the listener handle

[daq_data] = loadFromLogFile('log.dat',6); %load the data just saved to the dat file, using Stephen's function

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment9\' %move to our data directory
    
% if flyNum ==1 %if it's the first fly
%    mkdir ([date]) %make a folder with today's date
% end

if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment9\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment9\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment9\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end

save(strcat('ExpOptomotor',num2str(expNum),'.mat'),'daq_data','TrialNum', 'jumpFunction'); %save daq data and starting positions as a .mat file


end