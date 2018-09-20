%% analysisOfWalkingSpeed_Vm
clear all;
close all;
ephysSettings;

% pull out file names for the trials where the wanted stimulus was shown:
includeIfNameContainsString = true;

%wantedStimulusName = 'closedLoop_2bars_180deg' ;
%wantedStimulusName = 'closedLoop_2bars_90deg' ;
wantedStimulusName = 'closedLoop_vertStripeON' ;
%wantedStimulusName = 'constantCurrent' ;
%wantedStimulusName = 'closedLoop' ;

trialFilesList = extractTrialsWithCertainStimulusName( wantedStimulusName , includeIfNameContainsString);
%
%% Loop over all the files and extract the data for each epoch of closedloop stimulus
dataWholeExp = struct();

for fileNum = 1 : length ( trialFilesList )    
    % load current file for current trial
    load( trialFilesList(fileNum).name ); 
    
    
    if( fileNum == 1) % first trial set up dataWholeExp variable with subfields that match data
        dataWholeExp = data;
    else
    % concatenate all the fields in data.___ into a new large struct will
    % whole data set in it still organized by trial order
        dataWholeExp = [ dataWholeExp , data ];

    end
end

% Concatinate variables from whole single experiment 
voltage = cat(1, dataWholeExp(:).voltage);
xPanelPos = cat(1, dataWholeExp(:).xPanelPos);
ficTracAngularPosition = cat(1, dataWholeExp(:).ficTracAngularPosition);

% extract and concatinate Integrated X if exists
if(isfield( dataWholeExp ,'ficTracIntx') )
    ficTracIntx = cat(1, dataWholeExp(:).ficTracIntx);
end

% Median Filter Vm and plot data
MEDIAN_FILTER_WIDTH_SEC = 0.04; % sec
ORDER_MEDFILTER = MEDIAN_FILTER_WIDTH_SEC * settings.sampRate;   

medFilteredVoltage = medfilt1( voltage, ORDER_MEDFILTER , 'truncate' ); % Median filtering of the trace
medFilteredPanelPos = medfilt1( xPanelPos, ORDER_MEDFILTER, 'truncate' ); % Median filtering of the trace

% %% Plot data from whole experiment:
% FigHand = figure('Position',[50, 50, 1800, 800]);
% set(gcf, 'Color', 'w');
% 
 timeArray = (1  :  length(voltage) ) / settings.sampRate; % seconds
% 
% ax(1) = subplot(3,1,1);
% plot( timeArray,  voltage ,'-k', 'DisplayName' , 'membrane Voltage' ); hold on;
% %plot( timeArray, medFilteredVoltage, '-k');
% 
% ax(2) = subplot(3,1,2);
% plot( timeArray,  xPanelPos , 'DisplayName' , 'panel position' ); hold on;
% plot( timeArray, medFilteredPanelPos, '-k');
% 
% % plot position data from fictrac (0 - 10 Volts)
% ax(3) = subplot(3,1,3);
% plot( timeArray,  ficTracAngularPosition , 'DisplayName' , 'ficTrac heading (volts)' ); hold on;
% 
% linkaxes(ax,'x');
% legend('show')
%% :: ANGULAR VELOCITY :: Find and Plot Vm closed loop angular velocity
LOWPASS_FILTER_CUTOFF= 100; % Hz
THRESHOLD_ANGULAR_VELOCITY = 1500; % degrees / s  this is the max velocity that can be allowed into analysis
% decode angular velocity and accumulated position
[ angularVelocity , filteredFicTracPos ] = ficTracSignalDecoding( ficTracAngularPosition, settings.sampRate, LOWPASS_FILTER_CUTOFF, THRESHOLD_ANGULAR_VELOCITY);

% 
FigHand = figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');
 cx(1) = subplot(3, 1, 1);
 plot( timeArray, voltage, '-k', 'DisplayName' , 'membrane Voltage'  ); hold on;
 box off
% 
cx(2) = subplot(3, 1, 2);
 plot( timeArray, angularVelocity,  'DisplayName' , 'angular velocity'   ); hold on;
 box off
 
 cx(2) = subplot(3, 1, 3);
 plot( timeArray, medFilteredPanelPos,  'DisplayName' , 'panel position'   ); hold on;
 box off
 xlabel('sec')
 
 
linkaxes(cx,'x');
legend('show')

% Vm as a function of angular velocity
BIN_SIZE_FOR_ANGULAR_VELOCITY = 25;% degrees/s

possibleAngularVelocity =  min(angularVelocity) : BIN_SIZE_FOR_ANGULAR_VELOCITY : max(angularVelocity);

meanVoltageByVel = [];
stdVoltageByVel = [];
currLogical = [];
numberOfSamplesInduced = [];
secondsSpendAtVelocity = [];

for i = 1: (length ( possibleAngularVelocity ) -1)
    currLogical = ( possibleAngularVelocity(i) < angularVelocity) & (angularVelocity < possibleAngularVelocity (i + 1)) ;
    
    % store number of samples that exist for current bin
    numberOfSamplesInduced(i) = sum( currLogical );
    secondsSpendAtVelocity(i) = numberOfSamplesInduced(i) / settings.sampRate ; % convert samples to seconds
    
    % store mean voltage values
    meanVoltageByVel(i) = mean( medFilteredVoltage( currLogical ) );  
end

middleOfAngularVelocityValues = mean([possibleAngularVelocity(1:end-1);possibleAngularVelocity(2:end)]);
% account for sign flip from ficTrac directionality 
middleOfAngularVelocityValues = -1 * middleOfAngularVelocityValues;

%%
 meanVoltageByVelAboveThreshold = meanVoltageByVel;
figure;
set(gcf, 'Color', 'w');
%plot( middleOfAngularVelocityValues, meanVoltageByVel,'LineWidth', 1.5, 'Color',[0.5,0,0.1]); hold on;

DATA_THRESHOLD_High = 0.25 * settings.sampRate; % data from more than 0.25 sec of recording
DATA_THRESHOLD_Med = 0.05 * settings.sampRate; % data from more than 0.05 sec of recording
DATA_THRESHOLD_Low = 0.005 * settings.sampRate; % data from more than 0.005 sec of recording
%DATA_THRESHOLD =  5.0000e-05 * settings.sampRate; % data from more than 0.005 sec of recording

% % only plot values more often than threshold

% meanVoltageByVelAboveThreshold( numberOfSamplesInduced < DATA_THRESHOLD) = nan;
scatter(middleOfAngularVelocityValues(numberOfSamplesInduced> DATA_THRESHOLD_Low), meanVoltageByVelAboveThreshold(numberOfSamplesInduced> DATA_THRESHOLD_Low), 'filled', 'c',  'DisplayName' , 'data threshold 0.005s' ); hold on;
scatter(middleOfAngularVelocityValues(numberOfSamplesInduced> DATA_THRESHOLD_Med), meanVoltageByVelAboveThreshold(numberOfSamplesInduced> DATA_THRESHOLD_Med),'filled', 'b', 'DisplayName' , 'data threshold 0.05s' ); hold on;
scatter(middleOfAngularVelocityValues(numberOfSamplesInduced>DATA_THRESHOLD_High), meanVoltageByVelAboveThreshold(numberOfSamplesInduced> DATA_THRESHOLD_High),'filled', 'k',  'DisplayName' , 'data threshold 0.25s'); hold on;


%scatter ( middleOfAngularVelocityValues, meanVoltageByVelAboveThreshold, 'filled');
%scatter ( middleOfAngularVelocityValues, meanVoltageByVelAboveThreshold, 'filled');
 legend('show')
ylabel( 'Vm')
xlabel( 'angular velocity (deg / s)');
box off
niceaxes;
title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum)  ' stim: ' wantedStimulusName] );

%% :: ANGULAR SPEED :: Vm as a function of angular SPEED!
binSizeForAngularSpeed = 25;% degrees/s

angularSpeed = abs( angularVelocity );

possibleAngularSpeed =  min(angularSpeed): binSizeForAngularSpeed : max(angularSpeed);

meanVoltageBySpeed = [];
currLogical = [];
numberOfSamplesInduced = [];
secondsSpendAtSpeed = [];

for i = 1: (length ( possibleAngularSpeed ) -1)
    currLogical = ( possibleAngularSpeed(i) < angularSpeed) & (angularSpeed < possibleAngularSpeed (i + 1)) ;
    
    % store number of samples that exist for current bin
    numberOfSamplesInduced(i) = sum( currLogical );
    secondsSpendAtSpeed(i) = numberOfSamplesInduced(i) / settings.sampRate ; % convert samples to seconds
    
    % store mean voltage values
    meanVoltageBySpeed(i) = mean( medFilteredVoltage( currLogical ) );
end

middleOfAngularSpeedValues = mean([possibleAngularSpeed(1:end-1); possibleAngularSpeed(2:end)]);
%
figure; set(gcf, 'Color', 'w');

DATA_THRESHOLD = 0.25 * settings.sampRate; % data from more than 0.25 sec of recording
% only plot values more often than threshold
meanVoltageBySpeedAboveThreshold = meanVoltageBySpeed;
meanVoltageBySpeedAboveThreshold( numberOfSamplesInduced < DATA_THRESHOLD) = nan;

scatter ( middleOfAngularSpeedValues, meanVoltageBySpeedAboveThreshold, 'filled');  hold on;

box off
niceaxes;

% remove Nans:
y = meanVoltageBySpeedAboveThreshold( ~isnan(meanVoltageBySpeedAboveThreshold) );
x = middleOfAngularSpeedValues( ~isnan(meanVoltageBySpeedAboveThreshold) );

%Linear regression analysis:
[fitobj, gof ] = fit( x' ,y' , 'poly1' ); %lineaer Y = p1*x + p2

fitobj
gof

% plot linear fit on scatter
plot(fitobj); 
legend off;

ylabel( 'Vm')
xlabel( 'angular speed (deg / s)');
title( [ 'R-sq:' num2str( round(gof.rsquare,2)) '  ' num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum)] );


%% PLOT angular speed on the trace:
FigHand = figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');
 cx(1) = subplot(2, 1, 1);
 plot( timeArray, voltage, '-k', 'DisplayName' , 'membrane Voltage'  ); hold on;
 box off
% 
cx(2) = subplot(2, 1, 2);
 plot( timeArray, angularSpeed,  'DisplayName' , 'angular Speed (deg/s)'   ); hold on;
 box off
 xlabel('sec')
linkaxes(cx,'x');
legend('show')


%% :: FORWARD VELOCITY :: Find and Plot Vm closed loop forward velocity
LOWPASS_FILTER_CUTOFF= 100; % Hz
THRESHOLD_FORWARD_VELOCITY = 1500; % degrees / s  this is the max velocity that can be allowed into analysis

% decode forward velocity and accumulated X position
[ forwardVelocity , filteredIntX ] = ficTracSignalDecoding( ficTracIntx , settings.sampRate, LOWPASS_FILTER_CUTOFF, THRESHOLD_FORWARD_VELOCITY);

% Vm as a function of forward velocity
binSizeForForwardVelocity = 25;% degrees/s

possibleForwardVelocity =  min(forwardVelocity): binSizeForForwardVelocity : max(forwardVelocity);

meanVoltageByVel = [];
stdVoltageByVel = [];
currLogical = [];
numberOfSamplesInduced = [];
secondsSpendAtVelocity = [];

for i = 1: (length ( possibleForwardVelocity ) -1)
    currLogical = ( possibleForwardVelocity(i) < forwardVelocity) & (forwardVelocity < possibleForwardVelocity (i + 1)) ;
    
    % store number of samples that exist for current bin
    numberOfSamplesInduced(i) = sum( currLogical );
    secondsSpendAtVelocity(i) = numberOfSamplesInduced(i) / settings.sampRate ; % convert samples to seconds
    
    % store mean voltage values
    meanVoltageByVel(i) = mean( medFilteredVoltage( currLogical ) );
end

middleOfForwardVelocityValues = mean([possibleForwardVelocity(1:end-1);possibleForwardVelocity(2:end)]);
% account for sign flip from ficTrac directionality 
middleOfForwardVelocityValues = -1 * middleOfForwardVelocityValues;

DATA_THRESHOLD = 0.25 * settings.sampRate; % data from more than 0.25 sec of recording
% only plot values more often than threshold
meanVoltageByVelAboveThreshold = meanVoltageByVel;
meanVoltageByVelAboveThreshold( numberOfSamplesInduced < DATA_THRESHOLD) = nan;

%
figure; set(gcf, 'Color', 'w');
scatter ( middleOfForwardVelocityValues , meanVoltageByVelAboveThreshold, 'filled');

ylabel( 'Vm')
xlabel( 'forward velocity (deg / s)');
box off
niceaxes;
title( [ num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum)  ' stim: ' wantedStimulusName] );
%% :::  TRANSLATIONAL SPEED ::: speed: both forward and backward count:
binSizeForForwardSpeed = 25;% degrees/s

forwardSpeed = abs( forwardVelocity );

possibleForwardSpeed =  min(forwardSpeed): binSizeForForwardSpeed : max(forwardSpeed);

meanVoltageBySpeed = [];
currLogical = [];
numberOfSamplesInduced = [];
secondsSpendAtSpeed = [];

for i = 1: (length ( possibleForwardSpeed ) -1)
    currLogical = ( possibleForwardSpeed(i) < forwardSpeed) & (forwardSpeed < possibleForwardSpeed (i + 1)) ;
    
    % store number of samples that exist for current bin
    numberOfSamplesInduced(i) = sum( currLogical );
    secondsSpendAtSpeed(i) = numberOfSamplesInduced(i) / settings.sampRate ; % convert samples to seconds
    
    % store mean voltage values
    meanVoltageBySpeed(i) = mean( medFilteredVoltage( currLogical ) );
end

middleOfForwardSpeedValues = mean([possibleForwardSpeed(1:end-1); possibleForwardSpeed(2:end)]);


figure; set(gcf, 'Color', 'w');

DATA_THRESHOLD = 0.25 * settings.sampRate; % data from more than 0.25 sec of recording
% only plot values more often than threshold
meanVoltageBySpeedAboveThreshold = meanVoltageBySpeed;
meanVoltageBySpeedAboveThreshold( numberOfSamplesInduced < DATA_THRESHOLD) = nan;

scatter ( middleOfAngularSpeedValues, meanVoltageBySpeedAboveThreshold, 'filled');  hold on;


box off
niceaxes;

% remove Nans:
y = meanVoltageBySpeedAboveThreshold( ~isnan(meanVoltageBySpeedAboveThreshold) );
x = middleOfAngularSpeedValues( ~isnan(meanVoltageBySpeedAboveThreshold) );

%Linear regression analysis:
[fitobj, gof ] = fit( x' ,y' , 'poly1' ); %lineaer Y = p1*x + p2

fitobj
gof

% plot linear fit on scatter
plot(fitobj); 
legend off;

ylabel( 'Vm')
xlabel( 'forward speed (deg / s)');
title( [ 'R-sq:' num2str( round(gof.rsquare,2)) '  ' num2str(exptInfo.dNum) ' fly#: ' num2str(exptInfo.flyNum) ' cell#: '  num2str(exptInfo.cellNum) ' expt#: ' num2str(exptInfo.cellExpNum)] );

%%  Plotting slope of relationship, color coded for R2 value DATA FIRST
slopeRegressionWithAngularSpeed =100*[0.00161300000000000;-0.000504600000000000;0.00111000000000000;-0.00587500000000000;0.00484100000000000;-0.00150400000000000;0.00404800000000000;-0.00504500000000000;0.000537300000000000;0.000594700000000000;0.000765100000000000;0.000878600000000000;0.00217500000000000;0.00428900000000000;0.00559400000000000;0.000117900000000000;-0.000600100000000000;0.00119600000000000;0.00370800000000000];
Rsq_AngularSpeed = [0.980000000000000;0.626700000000000;0.790000000000000;0.821300000000000;0.993100000000000;0.552600000000000;0.955800000000000;0.899500000000000;0.566100000000000;0.598200000000000;0.958100000000000;0.555700000000000;0.994700000000000;0.959000000000000;0.991000000000000;0.232800000000000;0.831600000000000;0.988900000000000;0.963600000000000];

slopeRegressionWithTransSpeed = 100*[-0.000842300000000000;0.00160200000000000;-0.00561100000000000;0.00249600000000000;-0.00351300000000000;0.00355000000000000;-0.00246500000000000;0.0213200000000000;0.00175300000000000;0.00107400000000000;0.000258400000000000;0.00314200000000000];
Rsq_TransSpeed = [0.784400000000000;0.513600000000000;0.784100000000000;0.977900000000000;0.856400000000000;0.994200000000000;0.907300000000000;0.969900000000000;0.990400000000000;0.950400000000000;0.0776000000000000;0.995100000000000];
%%
figure;set(gcf, 'Color', 'w');
% positive slope, above Rsquared cut off
RSQUARED_CUTOFF = 0.5;

positiveSlope = slopeRegressionWithAngularSpeed( slopeRegressionWithAngularSpeed > 0);
RsquaredPos = Rsq_AngularSpeed( slopeRegressionWithAngularSpeed > 0);

postiveSlopeHighRsq = positiveSlope(RsquaredPos > RSQUARED_CUTOFF);
RsquaredPos = RsquaredPos(RsquaredPos > RSQUARED_CUTOFF);


s(2)= scatter( ones( 1, length( postiveSlopeHighRsq )) , postiveSlopeHighRsq,  50, RsquaredPos , 'filled'); hold on
s(2).MarkerEdgeColor = 'k'; hold on;
scatter( 1, mean(postiveSlopeHighRsq), 'k', 'filled');

colormap('parula')
box off;
niceaxes


c = colorbar;
c.Label.String = 'R-squared';
ylabel ('Slope  (mV / 100 deg/s)')

%% negative slope
figure;set(gcf, 'Color', 'w');
% positive slope, above Rsquared cut off
RSQUARED_CUTOFF = 0.5;

negativeSlope = slopeRegressionWithAngularSpeed( slopeRegressionWithAngularSpeed < 0);
RsquaredNeg = Rsq_AngularSpeed( slopeRegressionWithAngularSpeed < 0);

negativeSlopeHighRsq = negativeSlope (RsquaredNeg > RSQUARED_CUTOFF);
RsquaredNeg = RsquaredNeg( RsquaredNeg > RSQUARED_CUTOFF);

s2= scatter( ones( 1, length( negativeSlopeHighRsq )) , negativeSlopeHighRsq,  50, RsquaredNeg , 'filled'); hold on;
s2.MarkerEdgeColor = 'k';

scatter( 1, mean(negativeSlopeHighRsq), 'k', 'filled');

colormap('parula')
box off; niceaxes

c = colorbar;
c.Label.String = 'R-squared';
ylabel ('Slope  (mV / 100 deg/s)')

%%
figure;set(gcf, 'Color', 'w');
plot([ones(1, length( slopeRegressionWithAngularSpeed )); 2*ones(1, length( slopeRegressionWithTransSpeed ))], [ slopeRegressionWithAngularSpeed'  ; slopeRegressionWithTransSpeed' ], 'k');

hold on
s2= scatter( ones( 1, length( slopeRegressionWithAngularSpeed )) , slopeRegressionWithAngularSpeed,  140, Rsq_AngularSpeed , 'filled');
s2.MarkerEdgeColor = 'k';



hold on;
s3=scatter( 2*ones( 1, length( slopeRegressionWithTransSpeed )) , slopeRegressionWithTransSpeed,  140, Rsq_TransSpeed , 'filled');
s3.MarkerEdgeColor = 'k';
colormap('parula')
box off;
niceaxes

xlim([0 3])

c = colorbar;
c.Label.String = 'R-squared';
ylabel ('Slope  (mV / 100 deg/s)')
