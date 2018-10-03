%% Settings 3 px dark closed-loop bar, darker background

startPos = [round(rand*97) 1];
niOI.DurationInSeconds = 20;

Panel_com('set_pattern_id', 12); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_position',startPos);