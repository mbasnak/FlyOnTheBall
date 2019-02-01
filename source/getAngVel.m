function [angularVelocity] = getAngVel(data)

downsampled.aP = downsample(data,1000/25); %downsample
downsRad.aP = downsampled.aP .* 2 .* pi ./ 10; %change to radians
unwrapped.aP = unwrap(downsRad.aP); %unwrap
smoothed.aP = smoothdata(unwrapped.aP,10); %smooth
deg.aP = rad2deg(smoothed.aP); %change to degrees
diff.aP = gradient(deg.aP).* 25; %take the derivative
angularVelocity = smoothdata(diff.aP,10); %smooth again

end