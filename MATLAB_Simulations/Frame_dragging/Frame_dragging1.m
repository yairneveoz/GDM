% Frame_dragging1

clear
clc

% Parameters
R0 = 1; % Radius of the sphere
nTheta = 20; % Number of points along theta
nPhi = 20;   % Number of points along phi

% Spherical coordinates
theta = linspace(0, 2*pi, nTheta); % Azimuthal angle (0 to 2*pi)
phi = linspace(0, pi, nPhi);       % Polar angle (0 to pi)
[Theta, Phi] = meshgrid(theta, phi);

% Convert to Cartesian coordinates (sphere surface)
X = R0 * sin(Phi) .* cos(Theta);
Y = R0 * sin(Phi) .* sin(Theta);
Z = R0 * cos(Phi);

% Arrow components (radial direction)
Ux = X; % X-component of the arrow
Uy = Y; % Y-component of the arrow
Uz = Z; % Z-component of the arrow

% Plotting
figure(14);
hold on;
% Plot the sphere
surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'k'); % Transparent sphere
colormap('winter');
% Plot the arrows (radially outward)
quiver3(X, Y, Z, Ux, Uy, Uz, 0.5, 'r'); % Scale factor of 0.5 for arrow length
axis equal
view(30,30)









