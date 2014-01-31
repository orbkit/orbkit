%  MDC
%  Axel Schild (axel.schild [at] fu-berlin.de)
%  Gunter Hermann
%  Vincent Pohl
%  
%  
%  This file is part of MDC.
%  
%  MDC is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as 
%  published by the Free Software Foundation, either version 3 of 
%  the License, or any later version.
%  
%  MDC is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%  
%  You should have received a copy of the GNU Lesser General Public 
%  License along with MDC.  If not, see <http://www.gnu.org/licenses/>.


fid = 'h2o.h5';

% Load Grid
x = h5read(fid,'/x');
y = h5read(fid,'/y');
z = h5read(fid,'/z');

%  Load the density calculated
rho = h5read(fid,'/rho');

%  Do something with the data
d3r = (x(2,1)-x(1,1)) * (y(2,1)-y(1,1)) * (z(2,1)-z(1,1));
disp(sprintf('We have %.3f electrons.',(sum(sum(sum(rho)))*d3r)))