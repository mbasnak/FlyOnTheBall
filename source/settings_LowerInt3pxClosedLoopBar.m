%% Settings light, lower intensity, 3 px, closed-loop bar

startPos = [round(rand*97) 1];
niOI.DurationInSeconds = 40;

Panel_com('set_pattern_id', 14); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_position',startPos);