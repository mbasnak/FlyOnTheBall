function headingFig = analisisFicTracOutput(data)
%% this code reads the output.dat file generated by FicTrac and plots the heading of the fly over time

%INPUT
%data = .dat file MBoutput that came from FicTrac

flyHeading = table2array(data(:,17)); %convert to array the column where the heading data is
time = table2array(data(:,22)); %gives you back the time in absolute seconds acording to the processor
timeNorm = time-time(1); %gives you back the time in 0, taking the first element as 0 s.

headingFig = figure,
plot(timeNorm, flyHeading)
title('Heading of the fly in lab coordinates');
xlabel('Time (s)'); ylabel('Heading (lab coordinates)');

%I need to figure out what the lab coordinates actually mean
% I presume that heading = 0 will be given by wherever we've set our 
% front reference point.
% It's seems as though the heading is computed in rads with a max of 2pi


%% From here on I am trying to "unwrap" the position like Yvette does it, but I need to check this

unwrappedPos = unwrap(flyHeading);
% Q = unwrap(P) corrects the radian phase angles in a vector P by adding multiples
% of �2? when absolute jumps between consecutive elements of P are greater than or equal to the default jump tolerance of ? radians.

% find indexes where the unwrapping happened (tolerance = pi)
upwrappedIndexes = find ( abs( diff( flyHeading )) > pi); 
%I'm not following this

NUM_SAMPLES_FROM_WRAP_TO_REPLACE = 2;
% handle edge case so we don't fill off the edge of the trace
upwrappedIndexes = upwrappedIndexes( upwrappedIndexes > NUM_SAMPLES_FROM_WRAP_TO_REPLACE & upwrappedIndexes < (length ( unwrappedPos ) - NUM_SAMPLES_FROM_WRAP_TO_REPLACE) ); 

cleanedPos = unwrappedPos;
% replace potentially problematic indexes with Nan
for i = 1: length ( upwrappedIndexes )
    index_start = upwrappedIndexes(i) -  NUM_SAMPLES_FROM_WRAP_TO_REPLACE ; 
    index_end = upwrappedIndexes(i) +  NUM_SAMPLES_FROM_WRAP_TO_REPLACE ; 
    
    cleanedPos( index_start : index_end ) = NaN;
end

% replace NaN values with the last preceding value that was a real number
nanIDX = find( isnan( cleanedPos ) ); % find NaN indexes
% replace with preceeding value
while( ~isempty( nanIDX ) )
    cleanedPos(nanIDX) = cleanedPos(nanIDX - 1);
    
    % find any remaining NaN
    nanIDX  = find( isnan(cleanedPos) );
end

% transform from radians into degrees, send to user
accumulatedPositionOut = ( cleanedPos / (2*pi) ) * 360;

% take derivative and ajust for sample rate to solve for deg/s
%velocityOut = diff( accumulatedPositionOut ) .* sampleRate ; % degees / sec
velocityOut = gradient( accumulatedPositionOut ) .* 1000 ; % degees / sec

%low pass filter the velocity signal
%velocityOut = lowPassFilter( velocityOut, lowPassFilterCutOff, sampleRate );

% remove velocity values that are too large to be possible for the fly
velocityOut = replaceValuesOutsideThresholdBound( velocityOut, 900);

figure, plot(accumulatedPositionOut)

end