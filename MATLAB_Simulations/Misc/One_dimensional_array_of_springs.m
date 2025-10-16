% One_dimensional_array_of_springs

%{
Tracks the dynamics (locations over time) of an array of
points connected by 'springs' in one dimension.
%}

clear
clc

dx = 1; 
%{
Distance between two points when the spring is at 
its resting length
%}
dt = 0.01; % sec.
Nx = 100; % Number of points in x direction.

c_T = 3; % Transverse light velocity.
c_L = c_T*1.608; % Longitudinal light velocity.

size_x = 30;
x0 = -size_x*dx:dx:size_x*dx;
y0 = zeros(size(x0));

ux = 0.7;
uy = 0.7;
rho_0 = 1;

figure(3)
plot(x0, y0, '.-')







