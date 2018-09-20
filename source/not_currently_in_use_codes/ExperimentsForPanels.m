function runPanels(stimulus,time)
%experiments I will want to run with the panels
% I would like to make this a function where one can run the panels and
% choose the stimulus and underneath each stimulus is detailes.
% Something along the lines of runPanels(stimulus,time)
% if stimulus = 'open_loop_grating"...
%...Panel_com('set_pattern_id',2) etc


%1-open loop gratings for optomotor response
if stimulus == 'open_loop_grating'
    
Panel_com('set_pattern_id', 2); %the pattern 2 are the gratings
Panel_com('set_mode', [4 4]); %set the mode to be controlled by an outside function
Panel_com('set_position',[5 1]); %we should run this if we want to stimuli to start centered (doesn't make a lot of sense for the gratings though)
%maybe I should start at a random position? 
%we could do this by doing:
%Panel_com('set_position',[round(rand*97) 1]); this might make sense for
%the bar
Panel_com('start');
Pause(time); %this can later be changed by "time"
Panel_com('stop');


%2-closed loop light bar for stripe fixation
elseif stimulus == 'closed_loop_lightbar'

Panel_com('set_pattern_id', 4); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_position',[round(rand*97) 1]); %we can also comment this out, or start at [5 1]
Panel_com('start');
Pause(time);
Panel_com('stop')

end

end