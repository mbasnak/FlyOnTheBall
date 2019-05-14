function [daq_data] = ExpHighErrorRandomJumps(flyNum,expNum,TrialNum,folder)

% This function runs an experiment in which a grating is presented in closed-loop and
% every 30 s, the loop opens for 3 sec and the grating rotates randomly
% clockwise or counterclockwise

%INPUT
%flyNum : which number of fly are you running
%expNum : which experiment number are you running
%TrialNum : how many closed-loop trials + optomotor trials are you running
%folder: Experiment folder number in which this file should
%go. If it belongs in 'Experiment10', then folder will be 10.

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
niOI.DurationInSeconds = TrialNum*23+10; %set the duration to the total panel display time + 10 extra seconds

fid = fopen('log.dat','w+'); %this opens a csv file named "log",creating it for writing (and overwriting existing filed) with the flag "w+"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground(); %start acquiring


%%%%%% Run the panels in a loop %%%%%%
count = 1;
startPos = round(rand*96); %start stim at a random position


%stripe-fixation trial
    Panel_com('set_pattern_id', 13); %set the light bar as a pattern
    Panel_com('set_mode', [3 1]); %set it to closed-loop mode in the x dimension 
    pause(0.03)
    Panel_com('set_position',[startPos 2]);
    pause(0.03)
    Panel_com('set_AO',[3 32767]);
    Panel_com('start');
    pause(20) %time we want the bar in closed-loop (20 sec).

%code the jumps in a loop
for count = 1:TrialNum
    startPosition(count) = round(rand*96);
    Panel_com('stop');
    tic; %save the time when the panels stop
    Panel_com('set_AO',[3 16384]); %Im trying to set the V in the ON/OFF channel to an intermediate level to have a signal for when
    %there is a change in trial, but to not turn the pannels off.
    Panel_com('set_pattern_id', 13); %set the light bar as a pattern
    Panel_com('set_position',[startPosition(count) 2]); %the bar will jump to a new position
    %maybe I should add a function to keep the panels still, and that would
    %'reset' the closed loop
    pauseEndCLtoOpto(count) = toc;%save the time when the panels start to calculate how long they are stopped
    pause(3) %the panels stay frozen for 3 sec
    
    Panel_com('set_AO',[3 32767]);
    Panel_com('start');
    pause(20) %time we want the bar in closed-loop (20 sec).  

count = count+1;
end

Panel_com('stop');
Panel_com('all_off');
Panel_com('set_AO',[3 0]);
    
niOI.IsDone % will report 0
niOI.wait() % wait and keep acquiring until the acquisition time is over
niOI.IsDone % will report 1

delete(lh) % delete the listener handle

[daq_data] = loadFromLogFile('log.dat',6); %load the data just saved to the dat file, using Stephen's function

cd(strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\')); %move to our data directory
    
if flyNum == 1 && expNum == 1 %if it's the first fly and the first experiment
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

save(strcat('dataExpNum',num2str(expNum),'.mat'),'daq_data','TrialNum','pauseEndCLtoOpto','startPos'); %save daq data and starting positions as a .mat file


end