function [smoothed] = posDataDecoding(data, sampleRate)

%% (1) downsample

 downs.angularPosition = downsample(data.ficTracAngularPosition,sampleRate/25);
 downs.Intx = downsample(data.ficTracIntx,sampleRate/25);
 downs.Inty = downsample(data.ficTracInty,sampleRate/25);
% downs.yaw = downsample(data.yaw, sampleRate/25);


%% (2) Transform units

 downsRad.angularPosition = downs.angularPosition .* 2 .* pi ./ 10;
 %downsRad.yaw = downs.yaw .* 2 .* pi ./ 10;
 downsRad.Intx =  downs.Intx.* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq
 downsRad.Inty = downs.Inty .* 2 .* pi ./ 10;
 %% 2) unwrap
 
 unwrapped.angularPosition = unwrap(downsRad.angularPosition);
 %unwrapped.yaw = unwrap(downsRad.yaw);
 unwrapped.Intx = unwrap(downsRad.Intx);
 unwrapped.Inty = unwrap(downsRad.Inty);
 

    
    %% (3) low pass filter
 filteredPosition = lowPassFilter(unwrapped.angularPosition, 25, 1000);
 %filteredYaw = lowPassFilter(unwrapped.yaw, 25, 1000);
 filteredPositionx = lowPassFilter(unwrapped.Intx, 25, 1000);
 filteredPositiony = lowPassFilter(unwrapped.Inty, 25, 1000);

%(4) low pass filter again
filteredPosition2 = lowPassFilter(filteredPosition, 25, 1000);
%filteredYaw2 = lowPassFilter(filteredYaw, 25, 1000);
filteredPosition2x = lowPassFilter(filteredPositionx, 25, 1000);
filteredPosition2y = lowPassFilter(filteredPositiony, 25, 1000);


 smoothed.angularPosition = filteredPosition2;
 %smoothed.yaw = filteredYaw2;
 smoothed.Intx = filteredPosition2x;
 smoothed.Inty = filteredPosition2y;
%% (5) transform to degrees

deg.angularPosition = (filteredPosition2 / (2*pi)) * 360;
%deg.yaw = (filteredYaw2 / (2*pi)) * 360;
deg.xPosition = filteredPosition2x*4.5;
deg.yPosition = filteredPosition2y*4.5;

%% (6) take the derivative

smoothed.angularVel = gradient(deg.angularPosition).* 25;
%smoothed.yawVel = gradient(deg.yaw).* 25;
smoothed.xVel = gradient(deg.xPosition).* 25;
smoothed.yVel = gradient(deg.yPosition).* 25;


     %% Repeat the previous processes with the uncorrected heading for the trajectories
     
    if (isfield(data,'fictracAngularPosition') == 1)
        downs.AngularPosition = downsample(data.fictracAngularPosition,sampleRate/25);
        downsRad.AngularPosition =  downs.AngularPosition.* 2 .* pi ./ 10;
        unwrapped.AngularPosition = unwrap(downsRad.AngularPosition);
        filteredAngularPosition = lowPassFilter(unwrapped.AngularPosition, 25, 1000);
        filteredAngularPosition2 = lowPassFilter(filteredAngularPosition, 25, 1000);
        smoothed.AngularPosition = filteredAngularPosition2; 
        deg.AngularPosition = (filteredAngularPosition2 / (2*pi)) * 360;
        smoothed.AngularVel = gradient(deg.AngularPosition).* 25;
    end


end
