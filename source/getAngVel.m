function [angularVelocity] = getAngVel(data)

downsampled.aP = downsample(data,1000/25); %downsample
downsRad.aP = downsampled.aP .* 2 .* pi ./ 10; %change to radians
unwrapped.aP = unwrap(downsRad.aP); %unwrap
smoothed.aP = smoothdata(unwrapped.aP,'rlowess',25); %smooth
deg.aP = rad2deg(smoothed.aP); %change to degrees
diff.aP = gradient(deg.aP).* 25; %take the derivative

percentile25 = prctile(diff.aP,2.5);
percentile975 = prctile(diff.aP,97.5);
boundedDiffAngularPos = diff.aP;
boundedDiffAngularPos(boundedDiffAngularPos<percentile25 | boundedDiffAngularPos>percentile975) = NaN;
    
    [pointsVectorAV] = find(~isnan(boundedDiffAngularPos));
    valuesVectorAV = boundedDiffAngularPos(pointsVectorAV);
    xiAV = 1:length(boundedDiffAngularPos);
    interpAngVel = interp1(pointsVectorAV,valuesVectorAV,xiAV);
    
angularVelocity = smoothdata(interpAngVel,'rlowess',15); %smooth again

end