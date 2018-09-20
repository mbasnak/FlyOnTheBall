function runOpenLoopGrating(time)

Panel_com('set_pattern_id', 2); %the pattern 2 are the gratings
Panel_com('set_mode', [4 4]); %set the mode to be controlled by an outside function
Panel_com('set_position',[5 1]); %we should run this if we want to stimuli to start centered (doesn't make a lot of sense for the gratings though)
Panel_com('set_posfunc_id',[1 2]);
Panel_com('set_posfunc_id',[2 2]);
Panel_com('start');
pause(time);
Panel_com('stop');
Panel_com('all_off');

end