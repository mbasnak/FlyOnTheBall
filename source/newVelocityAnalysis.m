function [smoothed] = newVelocityAnalysis(data, sampleRate)
% This function processes the raw position data obtained from FicTrac to
% give back the velocity in degrees.

% It takes the rawData (data) and the sample rate (sampleRate) as inputs
% and gives back a struct names "smoothed" with the 3 final velocities for
% x, y and angular velocity.

% Unlike my previous function to obtain the velocity, this one doesn't
% smooth or removes extreme data


%% Downsample the position data to match FicTrac's output


    % Downsample to match FicTrac's output
    downsampled.Intx = downsample(data.ficTracIntx,sampleRate/25); %For a 1000 rate acquisition frame rate from the NiDaq, downsampling to 25 Hz equals taking 1 every 40 frames
    %downsampled.Inty = downsample(data.ficTracInty,sampleRate/25); %For a 1000 rate acquisition frame rate from the NiDaq, downsampling to 25 Hz equals taking 1 every 40 frames  
    downsampled.angularPosition = downsample(data.ficTracAngularPosition,sampleRate/25);
    

%% Tranform signal from voltage to radians for unwrapping

    downsRad.Intx = downsampled.Intx .* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq
    %downsRad.Inty = downsampled.Inty .* 2 .* pi ./ 10;
    downsRad.angularPosition = downsampled.angularPosition .* 2 .* pi ./ 10;

%% Unwrapping; I'm not smoothing, but I'm calling the data 'smoothed' to maintain the old name convention

    smoothed.Intx = unwrap(downsRad.Intx);
    %smoothed.Inty = unwrap(downsRad.Inty);
    smoothed.angularPosition = unwrap(downsRad.angularPosition);
    
     
%% Transform to useful systems 
    
    deg.Intx = smoothed.Intx * 4.5; % wer tranform the pos to mm by scaling the value by the sphere's radius
    %deg.Inty = smoothed.Inty * 4.5;
    deg.angularPosition = (smoothed.angularPosition / (2*pi)) * 360; % we transform the angular position to degrees
    
%% Take the derivative. Again, I'm changing the data name to 'smoothed' even though I'm not smoothing

    smoothed.xVel = gradient(deg.Intx).* 25; %we multiply by 25 because we have downsampled to 25 Hz
    %smoothed.yVel = gradient(deg.Inty).* 25; 
    smoothed.angularVel = gradient(deg.angularPosition).* 25;



end