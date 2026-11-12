% setup_paths.m  -  add the model source and test folders to the MATLAB path.
% Run once per MATLAB session from this folder:   >> setup_paths
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'src')));
addpath(fullfile(here, 'tests'));
fprintf('pbmr1d (MATLAB) paths added: %s\n', here);
