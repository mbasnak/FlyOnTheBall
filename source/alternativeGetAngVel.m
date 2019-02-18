function [angularVelocity] = alternativeGetAngVel(data)

downsampled.aP = downsample(data,1000/25); %downsample
downsRad.aP = downsampled.aP .* 2 .* pi ./ 10; %change to radians
unwrapped.aP = unwrap(downsRad.aP); %unwrap
deg.aP = rad2deg(unwrapped.aP); %change to degrees
diff.aP = filter(-smooth_diff(10),1,deg.aP)/0.04;
angularVelocity = smoothdata(diff.aP,10); %smooth again

end