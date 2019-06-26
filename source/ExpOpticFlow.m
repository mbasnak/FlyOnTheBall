function ExpOpticFlow(flyNum,expNum,time,folder,pattern)

%code to make a grating move from front to back when the animal moves
%forward and sideways when the animal rotates.

%INPUTS
%flyNum : number of fly you're running this day
%expNum : number that this experiment represents for this fly
%time : how long, in seconds, is the experiment running for
%folder : what folder number is the data going to be saving in
%pattern : what pattern number are you running, knowing that:
    %pattern 29 = grating with translational optic flow in 26 deg windows 180 deg apart
    %pattern 30 = starfield with translational optic flow in 26 deg windows 180 deg apart
    %pattern 31 = grating with translational optic flow in 52 deg windows 180 deg apart
    %pattern 32 = starfield with translational optic flow in 52 deg windows 180 deg apart
    %pattern 33 = starfield with translational optic flow in 52 deg windows
    %180 apart, and 2 'dark ends'
    


cd(strcat('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment',num2str(folder),'\'));

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
niOI.DurationInSeconds = time; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:7 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:7
    aI(i).InputType = 'SingleEnded';
end

Panel_com('set_pattern_id', pattern); %load the grating optic flow
pause(0.01)
Panel_com('set_mode', [3 3]); %set both dimensions to be in closed-loop with the animal's mvts.
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


save(strcat('dataOpticFlowHallwayWindow',num2str(expNum),'.mat'),'rawData','pattern'); %save as 


end