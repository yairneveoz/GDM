% Space_2D
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

% Define density values (adjustable)
rho_0 = 1;   % Reference density
rho_minus = 3; % Contracted space density (higher means stronger contraction)
rho_plus = 1 / (1/rho_0 - 1/rho_minus); % Compute inflated space density

% Define radial contraction function
contraction_factor = 1 - 0.8 * exp(-3 * (distances / max_distance).^2);
expansion_factor = 1 + 0.8 * exp(-3 * (distances / max_distance).^2);

% Apply radial contraction
x_contracted = center(1) + (x - center(1)) .* contraction_factor;
y_contracted = center(2) + (y - center(2)) .* contraction_factor;

% Apply radial inflation
x_inflated = center(1) + (x - center(1)) .* expansion_factor;
y_inflated = center(2) + (y - center(2)) .* expansion_factor;

% Create figure
figure(1);
tiledlayout(1,2); % Side-by-side plots

% Plot contracted space (blue grid)
nexttile;
hold on;
axis equal;
xlim([0 grid_size]);
ylim([0 grid_size]);
title('Radial Contracted Space');
set(gca, 'XTick', [], 'YTick', []);
for i = 1:size(x, 1)
    plot(x_contracted(i, :), y_contracted(i, :), 'b', 'LineWidth', 1.2);
end
for j = 1:size(y, 2)
    plot(x_contracted(:, j), y_contracted(:, j), 'b', 'LineWidth', 1.2);
end
hold off;

% Plot inflated space (red grid)
nexttile;
hold on;
axis equal;
xlim([0 grid_size]);
ylim([0 grid_size]);
title('Radial Inflated Space');
set(gca, 'XTick', [], 'YTick', []);
for i = 1:size(x, 1)
    plot(x_inflated(i, :), y_inflated(i, :), 'r', 'LineWidth', 1.2);
end
for j = 1:size(y, 2)
    plot(x_inflated(:, j), y_inflated(:, j), 'r', 'LineWidth', 1.2);
end
hold off;
