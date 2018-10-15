function [voltThresh] = assessNoise()

%This function reads a file that was ran without a fly on the ball to
%measure the noise in the voltage reading under that day's experimental
%conditions and gives back a reading of the fluctuations to later compare
%with the actual signals of the fly walking of the ball


% prompt the user to open the corresponding noise file
[file,path] = uigetfile();
data = load([path,file],'rawData');
data = data.rawData;

% calculate the differences between consecutive values in the voltage
% vectors to see the change per frame
changes = diff(data);

% calculate the mean difference per channel
meanChanges = mean(changes);

% calculate the std of the difference per channel
stdChanges = std(changes);

% determine the voltage change threshold as the mean + 2 std?
% and take the maximum value from the 3 data channels (1,2 and 3)

voltThresh = meanChanges+(2*stdChanges);
voltThresh = max(voltThresh(1:3));

end