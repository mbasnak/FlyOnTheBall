function [smoothed] = singleTrialVelocityAnalysis(data, sampleRate)
% This function processes the raw position data obtained from FicTrac to
% give back the velocity in degrees.

% It takes the rawData (data) and the sample rate (sampleRate) as inputs
% and gives back a struct names "smoothed" with the 3 final velocities for
% x, y and angular velocity.


%% Downsample the position data to match FicTrac's output


    % Downsample to match FicTrac's output
    downsampled.Intx = downsample(data.ficTracIntx,sampleRate/25); %For a 1000 rate acquisition frame rate from the NiDaq, downsampling to 25 Hz equals taking 1 every 40 frames
    downsampled.Inty = downsample(data.ficTracInty,sampleRate/25);
    downsampled.angularPosition = downsample(data.ficTracAngularPosition,sampleRate/25);
       
 
% The output is downsampled. It isn't noticeable when plotting solid lines, and
% it is barely noticeable when plotting dotted lines.

%% Tranform signal from voltage to radians for unwrapping

    downsRad.Intx = downsampled.Intx .* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq
    downsRad.Inty = downsampled.Inty .* 2 .* pi ./ 10;
    downsRad.angularPosition = downsampled.angularPosition .* 2 .* pi ./ 10;

% Now the position is going between 0 and 2 pi.

%% Unwrapping 

    unwrapped.Intx = unwrap(downsRad.Intx);
    unwrapped.Inty = unwrap(downsRad.Inty);
    unwrapped.angularPosition = unwrap(downsRad.angularPosition);

% Now the position is unwrapped, so it doesn't jump when moving from 0 to
% 2pi and vice versa

%% Smooth the data

    smoothed.Intx = smoothdata(unwrapped.Intx,10); %smooth using a 10 number window
    smoothed.Inty = smoothdata(unwrapped.Inty,10);
    smoothed.angularPosition = smoothdata(unwrapped.angularPosition,10);
    
%     smoothed.Intx = smooth(unwrapped.Intx,10); %smooth using a 10 number window
%     smoothed.Inty = smooth(unwrapped.Inty,10);
%     smoothed.angularPosition = smooth(unwrapped.angularPosition,10);
%     
    
    
    
%% Transform to useful systems 
    
    deg.Intx = smoothed.Intx * 3; % wer tranform the pos to mm by scaling the value by the sphere's radius
    deg.Inty = smoothed.Inty * 3;
    deg.angularPosition = (smoothed.angularPosition / (2*pi)) * 360; % we transform the angular position to degrees

    
%% Take the derivative

    diff.Intx = gradient(deg.Intx).* 25; %we multiply by 25 because we have downsampled to 25 Hz
    diff.Inty = gradient(deg.Inty).* 25;
    diff.angularPosition = gradient(deg.angularPosition).* 25;


    %%  Smooth again
 
    smoothed.xVel = smoothdata(diff.Intx,10);
    smoothed.yVel = smoothdata(diff.Inty,10);
    smoothed.angularVel = smoothdata(diff.angularPosition,5);
    
%     smoothed.xVel = smooth(diff.Intx,10);
%     smoothed.yVel = smooth(diff.Inty,10);
%     smoothed.angularVel = smooth(diff.angularPosition,5);

end