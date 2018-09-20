function runPanels(stimulus,time)
%experiments I will want to run with the panels
%for some reason it will only run the first stimulus and not the others
stim = stimulus;

%1-open loop gratings for optomotor response
if isequal(stim, 'open_loop_gratings')
    
Panel_com('set_pattern_id', 2); %the pattern 2 are the gratings
Panel_com('set_mode', [4 4]); %set the mode to be controlled by an outside function
Panel_com('set_position',[5 1]); %we should run this if we want to stimuli to start centered (doesn't make a lot of sense for the gratings though)
Panel_com('set_posfunc_id',[1 2]);
Panel_com('set_posfunc_id',[2 2]);
Panel_com('start');
pause(time);
Panel_com('stop');
Panel_com('all_off');


%2-closed loop light bar for stripe fixation
elseif isequal(stim, 'closed_loop_lightbar')

Panel_com('set_pattern_id', 4); %load the light stripe pattern
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('set_position',[round(rand*97) 1]); %we can also comment this out, or start at [5 1]
Panel_com('start');
pause(time);
Panel_com('stop');
Panel_com('all_off');



%3 open-loop bar

elseif isequal(stim, 'open_loop_lightbar')

Panel_com('set_pattern_id', 4); %load the light stripe pattern
Panel_com('set_mode', [4 4]); 
Panel_com('set_posfunc_id',[1 2]);
Panel_com('set_posfunc_id',[2 2]);
Panel_com('set_position',[round(rand*97) 1]); %we can also comment this out, or start at [5 1]
Panel_com('start');
pause(time);
Panel_com('stop');
Panel_com('all_off');



%4 from closed-loop to open loop bar

elseif isequal(stim, 'openToclosed_loop_lightbar')

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
    
%maybe I should also make the NI-DAQ acquisition into a function and add
%the time as I call the function. That way, I can run the runPanels first,
%and ask for a little more time, and then run the NI-DAQ and ask for a
%little less. Before all that, I need to start the FicTrac.

%I could think of a way of integrating both by running a single code? Will
%they both run at once otherwise?
end

end