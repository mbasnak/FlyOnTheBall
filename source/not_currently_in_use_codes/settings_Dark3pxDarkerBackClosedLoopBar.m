%% Settings 3 px dark closed-loop bar, darker background

startPos = [(round(rand*96)+1) 1];
niOI.DurationInSeconds = 40;

Panel_com('set_pattern_id', 13); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_position',startPos);