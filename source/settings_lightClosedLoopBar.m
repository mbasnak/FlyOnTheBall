%% Settings dark 4 px closed-loop bar

startPos = 1;
niOI.DurationInSeconds = 100;
Panel_com('set_pattern_id', 13); %load the light stripe pattern
Panel_com('set_position',[startPos,2]);  %set the starting position to a random location
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
