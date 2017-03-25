%adds all of the quadruped code to the path
display('adding code to path');
root = fileparts(mfilename('fullpath'));
addpath(genpath(root))
%addpath_gurobi
addpath_mosek
%addpath_iris