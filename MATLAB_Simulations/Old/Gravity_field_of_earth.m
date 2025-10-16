% Load the Earth Gravitational Model coefficients (EGM2008)
clear
clc

degree = (2:2159)'; % Degrees from 2 to 2159
order = ones(length(degree), 1); % Placeholder for order (set to 1 for simplicity)
Cnm = randn(length(degree), 1); % Random coefficients for Cnm
Snm = randn(length(degree), 1); % Random coefficients for Snm

% Save the simulated coefficients as a MATLAB .mat file
save('EGM2008coeff_estimated.mat', 'degree', 'order', 'Cnm', 'Snm');
% load('EGM2008coeff_estimated.mat'); % You'll need to have the coefficients file

% Set up a grid for latitude and longitude
n = 360; % Number of points for latitude and longitude
[lat, lon] = meshgrid(linspace(-90, 90, n), linspace(-180, 180, n));

% Convert latitude and longitude to spherical coordinates
r = 6378137; % Earth's mean radius in meters
phi = deg2rad(lat); % Convert latitude to radians
lambda = deg2rad(lon); % Convert longitude to radians

% Calculate the gravitational potential using spherical harmonics
V = 0; % Initialize potential
for n = 2:2159 % Iterate through the coefficients
%     Cnm = EGM2008coeff(n, 3); % Extract coefficients
%     Snm = EGM2008coeff(n, 4);
    Pnm0 = legendre(n-1, sin(phi(:)), 'sch')'; % Associated Legendre functions
    Pnm = reshape(Pnm0, 360,[]);
    V = V + r^n*(Cnm.*cos(lambda(:)*(n-1))+Snm.*sin(lambda(:)*(n-1))).*Pnm;
end

V = reshape(V, size(lat)); % Reshape potential to match the grid

% Plot the curvature of the gravity field
figure(4);
surf(lon, lat, -V); % Negative potential for visualization
title('Curvature of Earth''s Gravity Field (EGM2008)');
xlabel('Longitude');
ylabel('Latitude');
zlabel('Negative Gravitational Potential');
