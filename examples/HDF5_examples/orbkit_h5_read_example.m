%  This file is part of orbkit. See the main program or documentation 
%  for information on the license.
%  
%  Example file that shows how to read an HDF5 orbkit output with MATLAB.


fid = 'h2o.h5';

% Load grid
x = h5read(fid,'/x');
y = h5read(fid,'/y');
z = h5read(fid,'/z');

%  Load the density
rho = h5read(fid,'/rho');
% Attention! MATLAB interchanges the order of the indices
rho = permute(rho,[3,2,1])

%  Do something with the data
d3r = (x(2,1)-x(1,1)) * (y(2,1)-y(1,1)) * (z(2,1)-z(1,1));
disp(sprintf('We have %.3f electrons.',(sum(sum(sum(rho)))*d3r)))
