function RunFizedBarJumps(flyNum,expNum,trialNum,trialLength)
%%This function will run trials of a 2 px light bar on a dark number. They
%%will all be the same length, specified by "trialLength", and they will be
%%as many trials as is specified by "trialNum". The bar will jump by 180
%%degress from the animal's current position at the moment of the jump

% set up NiDaq acquisition session
daqreset 
devID = 'Dev1';  
niOI = daq.createSession('ni');
niOI.Rate = 1000;
aI = niOI.addAnalogInputChannel( devID , 1:6 , 'Voltage' );
for i = 1:5
    aI(i).InputType = 'SingleEnded';
end

niOI.DurationInSeconds = trialNum*trialLength + 10;

fid = fopen('log.dat','w+');
lh = niOI.addlistener('DataAvailable',@(src,event)logDaqData(fid,event));

niOI.startBackground();


startPosition = {};
rawData = {};


% Run the panels
Panel_com('set_pattern_id', 11); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop

% Set a for loop for making the bar jump 180 deg from current position
for i = 1:trialNum 

% I need to read the current position live and make it jump 180
startPosition{i} = [event.EventData(4,end) + 48]; % This will likely be a problem when the current position is bigger than 48?
% I think in background acquisition the data is stored in the struct called
% event, but I don't think I'll be able to read it that way.

Panel_com('set_position',[startPosition{i} 1]);
Panel_com('start');

end

Panel_com('stop');
Panel_com('set_A0',[3 0]);
Panel_com('all_off'); %turn off the panels


% Terminate the niDaq acquisition
niOI.IsDone
niOI.wait() 
niOI.IsDone 
delete(lh) 

[daq_data] = loadFromLogFile('log.dat',6); 


% Save the data
cd 'Z:\Wilson Lab\Mel\FlyOnTheBall\data';
    
if flyNum ==1 %if it's the first fly
   mkdir ([date]) %make a folder with today's date
end

if expNum == 1 %if it's the first experiment for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date]); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date,'\flyNum',num2str(flyNum)])
   getFlyInfo() %get fly's details
else
   cd (['Z:\Wilson Lab\Mel\FlyOnTheBall\data\',date,'\flyNum',num2str(flyNum)]) %otherwise move to this fly's folder
end

save(strcat('dataExpNum',num2str(expNum),'.mat'),'rawData','startPosition','trialLength'); %save as dataExpNum


end