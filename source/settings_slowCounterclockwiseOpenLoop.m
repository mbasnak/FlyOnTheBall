startPos = 50;
niOI.DurationInSeconds = 18;
Panel_com('set_pattern_id', 13);
pause(0.01)
Panel_com('set_position',[startPos,2]);
Panel_com('set_mode', [4 1]); 
Panel_com('set_funcx_freq', 50);
Panel_com('set_posfunc_id',[1 51]);
Panel_com('start');