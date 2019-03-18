%%This is a function to be able to plot the 2D trajectory of the fly
% I basically used the Cpp code to figure out how to obtain the variables
% posx and posy (integrated x/y position taking into account the changes in
% heading) using intx and inty (integrates x/y position without taking into
% account the heading) and the heading position.


function [posx,posy]=FlyTrajectory(intx,inty,heading)

heading_step = [0;diff(heading)];
heading_step = heading_step/4;
velx = [0;diff(intx)];
vely = [0;diff(inty)];
dir = [velx, vely];

phi = heading+heading_step/2.0;

speed = sqrt(velx.*velx+vely.*vely);
step = speed/4;

dir2x = cos(phi).*velx -sin(phi).*vely;
dir2y = sin(phi).*velx + cos(phi).*vely;

posx = zeros(1,length(intx));
posy = zeros(1,length(inty));

for i=1:length(intx)
posx(i+1) = posx(i)+step(i)*dir2x(i);
posy(i+1) = posy(i)+step(i)*dir2y(i);
end


end