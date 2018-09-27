%% Setings Open loop grating to the right

startPos = [round(rand*97) 1];
niOI.DurationInSeconds = 0.5;

Panel_com('set_pattern_id', 4); %the pattern 4 are the gratings
Panel_com('set_mode', [4 4]); %set the mode to be controlled by an outside function
Panel_com('set_position',startPos); %we should run this if we want to stimuli to start centered (doesn't make a lot of sense for the gratings though)
Panel_com('set_posfunc_id',[1 4]);
Panel_com('set_posfunc_id',[2 4]);