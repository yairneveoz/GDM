% Space_contraction_using_quaternions2

clear
clc

% Parameters for the grid
grid_size = 40;  % Number of points in each direction
x = linspace(-5, 5, grid_size);
y = linspace(-5, 5, grid_size);

% Create a meshgrid
[X, Y] = meshgrid(x, y);

% Define the scaling function (radial scaling example)
% Scaling factor: smaller near the origin, larger away
alpha = -0.5;   % Strength of scaling (negative for contraction)
sigma = 3;      % Range of influence for scaling
scaling_function = @(x, y) 1 + alpha * exp(-(x.^2 + y.^2) / sigma^2);

% Compute the scaling factors for each point in the grid
S = scaling_function(X, Y);

% Apply the scaling to the grid
X_scaled = S .* X;
Y_scaled = S .* Y;

% Plot the original and scaled grids
figure(13);
clf

% Original grid
subplot(1, 2, 1);
plot(X, Y, 'k'); hold on;
plot(X', Y', 'k'); % Transpose to plot vertical lines
title('Original Grid');
axis equal; grid on;

% Scaled grid
subplot(1, 2, 2);
plot(X_scaled, Y_scaled, 'r'); hold on;
plot(X_scaled', Y_scaled', 'r'); % Transpose to plot vertical lines
title('Scaled Grid (Localized Contraction)');
axis equal; grid on;

% Labels for clarity
sgtitle('Localized Contraction of a 2D Grid');
