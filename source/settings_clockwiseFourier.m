startPos = 21;
niOI.DurationInSeconds = 9.3;
Panel_com('set_pattern_id', 23);
pause(0.01)
Panel_com('set_position',[startPos,2]);
Panel_com('set_mode', [4 1]); 
Panel_com('set_funcx_freq', 50);
Panel_com('set_posfunc_id',[1 157]);
Panel_com('start');