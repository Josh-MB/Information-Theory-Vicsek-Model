% Startup stuff

global vicsek_version;
vicsek_version.major = 0;
vicsek_version.minor = 1;

fprintf('[vicsek startup] Initialising vicsek version %d.%d\n', vicsek_version.major, vicsek_version.minor);

% Add vicsek root dir + appropriate subdirs to path

global vicsek_root;
vicsek_root = fileparts(mfilename('fullpath')); % directory containing this file
addpath(vicsek_root);
fprintf('[vicsek startup] Added path %s\n',vicsek_root);

%{
% Initialize mvgc library

global mvgc_root;
mvgc_root = '../../ncomp/mvgc'; % amend to location of mvgc library
cd(mvgc_root);
startup;
cd(vicsek_root);
%}

%{
% Initialize utils library

global utils_root;
utils_root = '../../utils'; % amend to location of utils
addpath(utils_root);
%}

%set(0,'DefaultFigureToolBar','none','DefaultFigureMenuBar','none'); % this stuff just takes up space
set(0,'DefaultFigureToolBar','none');

fprintf('[vicsek startup] Initialised (you may re-run ''startup'' at any time)\n');
