% Space_contraction_using_quaternions1

clear
clc

% Parameters for the grid
grid_size = 21;  % Number of points in each direction
x = linspace(-5, 5, grid_size);
y = linspace(-5, 5, grid_size);

% Create a meshgrid
[X, Y] = meshgrid(x, y);

% Define the scaling function (radial scaling example)
% Scaling factor: larger near the origin, decaying away
alpha = 1.5;   % Strength of scaling
sigma = 3;     % Range of influence for scaling
scaling_function = @(x, y) 1 + alpha * exp(-(x.^2 + y.^2) / sigma^2);

% Compute the scaling factors for each point in the grid
S = scaling_function(X, Y);

% Apply the scaling to the grid
X_scaled = S .* X;
Y_scaled = S .* Y;

% Plot the original and scaled grids
figure(12);
clf
% Original grid
subplot(1, 2, 1);
plot(X, Y, 'k'); hold on;
plot(X', Y', 'k'); % Transpose to plot vertical lines
title('Original Grid');
axis equal; grid on;

% Scaled grid
subplot(1, 2, 2);
plot(X_scaled, Y_scaled, 'b'); hold on;
plot(X_scaled', Y_scaled', 'b'); % Transpose to plot vertical lines
title('Scaled Grid (Localized Scaling)');
axis equal; grid on;

% Labels for clarity
sgtitle('Localized Contraction/Inflation of a 2D Grid');
