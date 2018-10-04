%% Setings 3 px, darker Open loop grating to the left

startPos = [round(rand*97) 1];
niOI.DurationInSeconds = 1;

Panel_com('set_pattern_id', 3); %the pattern 4 are the gratings
Panel_com('set_mode', [4 4]); %set the mode to be controlled by an outside function
Panel_com('set_position',startPos); %we should run this if we want to stimuli to start centered (doesn't make a lot of sense for the gratings though)
Panel_com('set_posfunc_id',[1 2]); %this function moves the grating to the right at 20 Hz.
Panel_com('set_posfunc_id',[2 2]);