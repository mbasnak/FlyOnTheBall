function runClosedLoopToOpenLoop(time)

Panel_com('set_pattern_id', 4); %load the light stripe pattern
Panel_com('set_mode', [4 4]); 
Panel_com('set_posfunc_id',[1 2]);
Panel_com('set_posfunc_id',[2 2]);
Panel_com('set_position',[round(rand*97) 1]); %we can also comment this out, or start at [5 1]
Panel_com('start');
pause(time/2);
Panel_com('stop');
Panel_com('set_mode', [3 0]);
Panel_com('start');
pause(time/2);
Panel_com('stop');
Panel_com('all_off');

end