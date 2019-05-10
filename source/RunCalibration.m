function [data] = RunCalibration(time)

%This function runs a calibration to acquire data from a rotary encoder
%that carries a signal from a ball being rotated by a servo-motor, as well
%as the fictrac output data, to see how well they match.


%INPUTS
%time = how long, in seconds, you are going to run the calibration for

%OUTPUTS
%data = array with the voltage data from both fictrac and the rotary
%encoder

% set up NiDaq acquisition session
daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:5 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:5
    aI(i).InputType = 'SingleEnded';
end

niOI.DurationInSeconds = time;

data = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined


end