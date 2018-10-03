
function [smoothed] = velocityAnalysis(data, sampleRate)
% This function processes the raw position data obtained from FicTrac to
% give back the velocity in degrees.

% It takes the rawData (data) and the sample rate (sampleRate) as inputs
% and gives back a struct names "smoothed" with the 3 final velocities for
% x, y and angular velocity.


%% Downsample the position data to match FicTrac's output

%figure,

for i = 1:size(data.ficTracIntx,2)
    % Downsample to match FicTrac's output
    downsampled.Intx{i} = downsample(data.ficTracIntx{i},sampleRate/25); %For a 1000 rate acquisition frame rate from the NiDaq, downsampling to 25 Hz equals taking 1 every 40 frames
    downsampled.Inty{i} = downsample(data.ficTracInty{i},sampleRate/25);
    downsampled.angularPosition{i} = downsample(data.ficTracAngularPosition{i},sampleRate/25);
       
% subplot(1,3,1)
% plot(downsampled.Intx{i})
% hold on
% title('Downsampled x position')
% subplot(1,3,2)
% plot(downsampled.Inty{i})
% hold on
% title('Downsampled y position')
% subplot(1,3,3)
% plot(downsampled.angularPosition{i})
% hold on
% title('Downsampled angular position')
% hold on
    
end   

% The output is downsampled. It isn't noticeable when plotting solid lines, and
% it is barely noticeable when plotting dotted lines.

%% Tranform signal from voltage to radians for unwrapping

%figure,

for i = 1:size(data.ficTracIntx,2)
    
    downsRad.Intx{i} = downsampled.Intx{i} .* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq
    downsRad.Inty{i} = downsampled.Inty{i} .* 2 .* pi ./ 10;
    downsRad.angularPosition{i} = downsampled.angularPosition{i} .* 2 .* pi ./ 10;
    
% subplot(1,3,1)
% plot(downsRad.Intx{i})
% hold on
% title('Downsampled radians x position')
% subplot(1,3,2)
% plot(downsRad.Inty{i})
% hold on
% title('Downsampled radians y position')
% subplot(1,3,3)
% plot(downsRad.angularPosition{i})
% hold on
% title('Downsampled radians angular position')
% hold on   
    
end

% Now the position is going between 0 and 2 pi.

%% Unwrapping 

%figure,

for i = 1:size(data.ficTracIntx,2)

    unwrapped.Intx{i} = unwrap(downsRad.Intx{i});
    unwrapped.Inty{i} = unwrap(downsRad.Inty{i});
    unwrapped.angularPosition{i} = unwrap(downsRad.angularPosition{i});
    
% subplot(1,3,1)
% plot(unwrapped.Intx{i})
% hold on
% title('Unwrapped x position')
% subplot(1,3,2)
% plot(unwrapped.Inty{i})
% hold on
% title('Unwrapped y position')
% subplot(1,3,3)
% plot(unwrapped.angularPosition{i})
% hold on
% title('Unwrapped angular position')
% hold on  
    
end

% Now the position is unwrapped, so it doesn't jump when moving from 0 to
% 2pi and vice versa

%% Smooth the data

%figure,

for i = 1:size(data.ficTracIntx,2)
    smoothed.Intx{i} = smooth(unwrapped.Intx{i},10); %smooth using a 10 number window
    smoothed.Inty{i} = smooth(unwrapped.Inty{i},10);
    smoothed.angularPosition{i} = smooth(unwrapped.angularPosition{i},10);
    
% subplot(1,3,1)
% plot(smoothed.Intx{i})
% hold on
% title('Smoothed x position')
% subplot(1,3,2)
% plot(smoothed.Inty{i})
% hold on
% title('Smoothed y position')
% subplot(1,3,3)
% plot(smoothed.angularPosition{i})
% hold on
% title('Smoothed angular position')
% hold on 
end
    
    
%% Transform back to degrees 

%figure, 

for i = 1:size(data.ficTracIntx,2)
    
    deg.Intx{i} = (smoothed.Intx{i} / (2*pi)) * 360;
    deg.Inty{i} = (smoothed.Inty{i} / (2*pi)) * 360;
    deg.angularPosition{i} = (smoothed.angularPosition{i} / (2*pi)) * 360;

% subplot(1,3,1)
% plot(deg.Intx{i})
% hold on
% title('Unwrapped x position in degrees')
% subplot(1,3,2)
% plot(deg.Inty{i})
% hold on
% title('Unwrapped y position in degrees')
% subplot(1,3,3)
% plot(deg.angularPosition{i})
% hold on
% title('Unwrapped angular position in degrees')
% hold on 

end

    
%% Take the derivative
    
%figure,

for i = 1:size(data.ficTracIntx,2)
    
    diff.Intx{i} = gradient(deg.Intx{i}).* 25; %we multiply by 25 because we have downsampled to 25 Hz
    diff.Inty{i} = gradient(deg.Inty{i}).* 25;
    diff.angularPosition{i} = gradient(deg.angularPosition{i}).* 25;
    
% subplot(1,3,1)
% plot(diff.Intx{i})
% hold on
% title('x Velocity in degrees/s')
% subplot(1,3,2)
% plot(diff.Inty{i})
% hold on
% title('y velocity in degrees/s')
% subplot(1,3,3)
% plot(diff.angularPosition{i})
% hold on
% title('angular velcoty in degrees/s')
% hold on    
%     
end

    %%  Smooth again
 
%figure,

for i=1:size(data.ficTracIntx,2)
    
    smoothed.xVel{i} = smooth(diff.Intx{i},10);
    smoothed.yVel{i} = smooth(diff.Inty{i},10);
    smoothed.angularVel{i} = smooth(diff.angularPosition{i},10);
  
% subplot(1,3,1)
% plot(smoothed.xVel{i})
% hold on
% title('Smoothed x Velocity in degrees/s')
% subplot(1,3,2)
% plot(smoothed.xVel{i})
% hold on
% title('Smoothed y velocity in degrees/s')
% subplot(1,3,3)
% plot(smoothed.xVel{i})
% hold on
% title('Smoothed angular velcoty in degrees/s')
% hold on  

end

end