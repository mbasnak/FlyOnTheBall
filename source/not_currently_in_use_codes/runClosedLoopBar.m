function runClosedLoopBar(time)

Panel_com('set_pattern_id', 4); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_position',[round(rand*97) 1]); %we can also comment this out, or start at [5 1]
Panel_com('start');
pause(time);
Panel_com('stop');
Panel_com('all_off');

end