% Geometrodynamic_space
clear
clc

grid_size = 10; % number of points in each dimention.
a0 = 1; % lattice constant at rest
[grid_x, grid_y, grid_z] = meshgrid(1:grid_size);

% Gaussians:
sigma_x = 2;
sigma_y = 2;
sigma_z = 2;

mu_x = grid_size/2;
mu_y = grid_size/2;
mu_z = grid_size/2;

gauss_x = exp(-0.5*((grid_x - mu_x)/sigma_x).^2);
gauss_y = exp(-0.5*((grid_y - mu_y)/sigma_y).^2);
gauss_z = exp(-0.5*((grid_z - mu_z)/sigma_z).^2);

points_x = grid_x + 0.5*(gauss_y.^2 + gauss_z.^2);
points_y = grid_y + 0.5*(gauss_x.^2 + gauss_z.^2);
points_z = grid_z + 0.5*(gauss_x.^2 + gauss_y.^2);
% Unit cell


% Shift (x,y,z)
% grid_dx = diff(grid_x,1,2);
% grid_dy = diff(grid_y,1,1);
% grid_dz = diff(grid_z,1,3);
%
%%%

figure(17)
% scatter3(points_x(:), points_y(:), points_z(:),'.')
plot3(points_x(:), points_y(:), points_z(:),'.')
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal