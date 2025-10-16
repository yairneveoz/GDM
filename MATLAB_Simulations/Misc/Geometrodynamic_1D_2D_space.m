% Geometrodynamic_1D_space
clear
clc

inputs.grid_size = 41; % number of points in each dimention.
inputs.a0 = 1; % lattice constant at rest
inputs.rho_0 = 1;
inputs.rho_max = 4;

inputs.cT = 3e10; % cm/sec
inputs.alpha = 1/137.04;
inputs.cL = (pi/2)*(1 + pi*inputs.alpha)*inputs.cT; % cm/sec
 
% Gaussians:
inputs.sigma_r = inputs.grid_size/20;
inputs.Nc = 64;
% (1 - rho_0/rho_p) = -(1 - rho_0/rho_n);
% 2/rho_0 = 1/rho_n + 1/rho_p
% 
oneDoutputs = get1D(inputs);
plot1D(inputs,oneDoutputs)

twoDoutputs = get2D(inputs);
plot2D(inputs,twoDoutputs)
% plot2D_dynamics_XY(inputs)
% plot2D_dynamics_YZ(inputs)
%
%%