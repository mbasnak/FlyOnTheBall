%% Setings starfield to the left

startPos = [(round(rand*96)+1) 1];
niOI.DurationInSeconds = 1;

Panel_com('set_pattern_id', 17); %the pattern 4 are the gratings
Panel_com('set_mode', [4 4]); %set the mode to be controlled by an outside function
Panel_com('set_position',startPos); %we should run this if we want to stimuli to start centered (doesn't make a lot of sense for the gratings though)
Panel_com('set_posfunc_id',[1 2]); %This function moves the pattern to the left at 20 Hz
Panel_com('set_posfunc_id',[2 2]);