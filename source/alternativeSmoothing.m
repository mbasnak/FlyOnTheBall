function [alternativeSmoothed] = alternativeSmoothing(data, sampleRate)

downsampled.Intx = downsample(data.ficTracIntx,sampleRate/25); %For a 1000 rate acquisition frame rate from the NiDaq, downsampling to 25 Hz equals taking 1 every 40 frames
downsampled.Inty = downsample(data.ficTracInty,sampleRate/25);
downsampled.angularPosition = downsample(data.ficTracAngularPosition,sampleRate/25);

downsRad.Intx = downsampled.Intx .* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq
downsRad.Inty = downsampled.Inty .* 2 .* pi ./ 10;
downsRad.angularPosition = downsampled.angularPosition .* 2 .* pi ./ 10;

unwrapped.Intx = unwrap(downsRad.Intx);
unwrapped.Inty = unwrap(downsRad.Inty);
alternativeSmoothed.angularPosition = unwrap(downsRad.angularPosition);

deg.Intx = unwrapped.Intx * 4.75; % wer tranform the pos to mm by scaling the value by the sphere's radius
deg.Inty = unwrapped.Inty * 4.75;
deg.angularPosition = (alternativeSmoothed.angularPosition / (2*pi)) * 360;

alternativeSmoothed.xVel = filter(-smooth_diff(10),1,deg.Intx)/0.04;
alternativeSmoothed.yVel = filter(-smooth_diff(10),1,deg.Inty)/0.04;
alternativeSmoothed.angularVel = filter(-smooth_diff(10),1,deg.angularPosition)/0.04;

end