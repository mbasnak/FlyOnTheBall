function [daq_data] = ExpYawStabilization_grating_yoked(flyNum,expNum,TrialNum,folder)

% This function runs an experiment in which the closed-loop experience that a master fly had is
% presented in open-loop for 30 sec, and then optomotor response trialss
% 3 sec long are interspersed, rotating randomly clockwise or
% counterclockwise

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
niOI.DurationInSeconds = TrialNum*33+10; %set the duration to the total panel display time + 10 extra seconds

fid = fopen('log.dat','w+'); %this opens a csv file named "log",creating it for writing (and overwriting existing filed) with the flag "w+"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground(); %start acquiring


%%%%%% Run the panels in a loop %%%%%%
count = 1;
jumpFunction = zeros(1,TrialNum); %I'll store the jumpFunction for the optomotor response to know the direction

masterFly = round(rand+1); %randomly choose masterFly1 or masterFly2
%determine the yokedFunction matrix based on the master fly
if masterFly == 1
    yokedFunction = [56:105];
else
    yokedFunction = [106:155];
end

for count = 1:TrialNum
    %yoked open-loop trial
    Panel_com('set_pattern_id', 22); %set the grating as a pattern
    Panel_com('set_mode', [4 0]); %set it to closed-loop mode in the x dimension 
    pause(0.03)
    Panel_com('set_posfunc_id',[1 yokedFunction(count)])
    pause(0.03)
    Panel_com('set_position',[1,2]);
    Panel_com('set_funcx_freq', 25);
    pause(0.03)
    Panel_com('set_AO',[3 32767]);
    Panel_com('start');
    pause(30) %record for the time we want the grating in closed-loop (30 sec).
    Panel_com('stop');
    tic; %save the time when the panels stop
    Panel_com('set_AO',[3 16384]); %Im trying to set the V in the ON/OFF channel to an intermediate level to have a signal for when
    %there is a change in trial, but to not turn the pannels off.
    
    %optomotor response trial
    jumpFunction(count) = round(rand+52); %toggle between clockwise (function 50) and counterclockwise (function 51)
    Panel_com('set_pattern_id', 19); %set the grating
    Panel_com('set_mode', [4 1]); %set it to function-driven in the x dimension 
    Panel_com('set_position',[1 1]);
    Panel_com('set_funcx_freq', 50);
    Panel_com('set_posfunc_id',[1 jumpFunction(count)]);
    Panel_com('set_AO',[3 32767]);
    pauseEndCLtoOpto(count) = toc; %save the time when the panels start to calculate how long they are stopped
    Panel_com('start');
    pause(3) %record for the time we want the grating in open-loop
    Panel_com('stop');
    Panel_com('set_AO',[3 16384]);

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

save(strcat('yokedDataExpNum',num2str(expNum),'.mat'),'daq_data','TrialNum','jumpFunction','pauseEndCLtoOpto','masterFly','yokedFunction'); %save daq data and starting positions as a .mat file


end