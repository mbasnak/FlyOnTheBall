function FlyData = getFlyInfo(varargin)
%{
GETFLYDETAILS Used in the case of a new fly to get information from the user about this
fly/experimental parameters/particulars of the dissection

INPUT



OUTPUT

FlyData (struct)
.line (genotype/Gal4/effector...)
.eclosionDate

All of these subfields are obtained from the user and saved
%}

%% Ask user for input
FlyData.line = input('Line: ','s');

% Get eclosion date
FlyData.eclosionDate = input('Ecclosion date: ','s');

% Get starvation time
FlyData.sex = input('Starved since: ','s');

% Get head status
FlyData.sex = input('Glued head?: ','s');

% Get state of the wings
FlyData.sex = input('Wings state?: ','s');

% Save
save('flyData','FlyData')
