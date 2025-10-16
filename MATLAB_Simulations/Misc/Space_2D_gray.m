% Space_2D_gray
clc; clear; close all;

% Grid size
grid_size = 40;

% Define contraction center
center = [grid_size / 2, grid_size / 2];

% Create original grid points
[x, y] = meshgrid(0:grid_size, 0:grid_size);

% Compute distance from center
distances = sqrt((x - center(1)).^2 + (y - center(2)).^2);
max_distance = max(distances(:)); % Normalize distances

% Define grayscale intensity reference
rho_0 = 0.5;  % Medium gray

% Define radial contraction function (for rho_+ → White)
contraction_factor = exp(-3 * (distances / max_distance).^2); % Strongest at center

% Define radial inflation function (for rho_- → Black)
expansion_factor = 1 - exp(-3 * (distances / max_distance).^2); % Strongest at center

% Compute grayscale values
grayscale_contracted = rho_0 + (1 - rho_0) * contraction_factor; % White at center
grayscale_inflated = rho_0 - rho_0 * contraction_factor; % Black at center

% Create figure
figure(2);
tiledlayout(1,2); % Side-by-side plots

% Plot contracted space (goes to white at the center)
nexttile;
imagesc(grayscale_contracted);
colormap(gray);
axis equal;
xlim([1 grid_size+1]);
ylim([1 grid_size+1]);
title('Contracted Space (\rho_+ → White)');
set(gca, 'XTick', [], 'YTick', []);

% Plot inflated space (goes to black at the center)
nexttile;
imagesc(grayscale_inflated);
colormap(gray);
axis equal;
xlim([1 grid_size+1]);
ylim([1 grid_size+1]);
title('Inflated Space (\rho_- → Black)');
set(gca, 'XTick', [], 'YTick', []);
