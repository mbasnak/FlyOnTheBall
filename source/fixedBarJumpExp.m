function fixedBarJumpExp(flyNum,expNum)


%%%%%%%%  Initialize the DAC session %%%%%%%%
daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
aI = niOI.addAnalogInputChannel( devID , 1:6 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:6
    aI(i).InputType = 'SingleEnded';
end
% niOI.DurationInSeconds = 10*200 + 10; %set the duration to the total panel display time + 5 extra seconds
niOI.DurationInSeconds = 270; %the duration of 12 trials of 200 sec plus 10 extra sec

fid = fopen('log.dat','w+'); %this opens a csv file named "log",creating it for writing (and overwriting existing filed) with the flag "w+"
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));
%This listener tells matlab to save the data in the file "fid" once the data is available

niOI.startBackground(); %start acquiring



%%%%%% Run the panels %%%%%%

jumps = repelem([-90,90,180],4);
jumps = jumps(randperm(length(jumps))); %randomize the order of the 3 possible jumps
bias = jumps/1.8; %assign the bias applying the proper transformation to the jumps vector
bias = [0,bias]; %add a 0 at the beginning so that the first time, there is no offset

% Initialize the pattern
Panel_com('set_pattern_id', 6); %set the bar
Panel_com('set_mode', [3 0]); %set it to closed-loop mode in the x dimension and to be controlled by a function in the y dimension 
startPos = round(rand*96);
Panel_com('set_position',[startPos 1]);
Panel_com('set_AO',[3 32767]);

% run in a for loop the 'trials'
for i = 1:13
    Panel_com('send_gain_bias', [10 bias(i) 0 0]);
    pause(0.03)
    Panel_com('start');
    pause(20);
    %pause(200)
    Panel_com('stop');
end

% Stop the panels
Panel_com('set_AO',[3 0]);
Panel_com('all_off');


niOI.IsDone % will report 0
niOI.wait() % wait and keep acquiring until the acquisition time is over
niOI.IsDone % will report 1

delete(lh) % delete the listener handle



%%%%%%%%% Recover and save the data in the right folders %%%%%%%%%%%
[daq_data] = loadFromLogFile('log.dat',6); %load the data just saved to the dat file, using Stephen's function

cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment2\' %move to our data directory
    
if flyNum ==1 %if it's the first fly
   mkdir ([date]) %make a folder with today's date
end


if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment2\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment2\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment2\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end

save(strcat('dataExpNum',num2str(expNum),'.mat'),'daq_data','startPos','jumps'); %save daq data and starting positions as a .mat file

end