% space_2D_gray_and_grid
clc; clear; close all;

% Grid size
grid_size = 10;

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

% Apply radial contraction (shrinking space)
x_contracted = center(1) + (x - center(1)) .* contraction_factor;
y_contracted = center(2) + (y - center(2)) .* contraction_factor;

% Apply radial inflation (expanding space)
x_inflated = center(1) + (x - center(1)) .* expansion_factor;
y_inflated = center(2) + (y - center(2)) .* expansion_factor;

% Create figure
figure;
tiledlayout(1,2); % Side-by-side plots

% Plot contracted space (grid + grayscale)
nexttile;
hold on;
imagesc(grayscale_contracted);
colormap(gray);
axis equal;
xlim([1 grid_size+1]);
ylim([1 grid_size+1]);
title('Contracted Space (\rho_+ → White)');
set(gca, 'XTick', [], 'YTick', []);
for i = 1:size(x, 1)
    plot(x_contracted(i, :), y_contracted(i, :), 'k', 'LineWidth', 1.2);
end
for j = 1:size(y, 2)
    plot(x_contracted(:, j), y_contracted(:, j), 'k', 'LineWidth', 1.2);
end
hold off;

% Plot inflated space (grid + grayscale)
nexttile;
hold on;
imagesc(grayscale_inflated);
colormap(gray);
axis equal;
xlim([1 grid_size+1]);
ylim([1 grid_size+1]);
title('Inflated Space (\rho_- → Black)');
set(gca, 'XTick', [], 'YTick', []);
for i = 1:size(x, 1)
    plot(x_inflated(i, :), y_inflated(i, :), 'k', 'LineWidth', 1.2);
end
for j = 1:size(y, 2)
    plot(x_inflated(:, j), y_inflated(:, j), 'k', 'LineWidth', 1.2);
end
hold off;
