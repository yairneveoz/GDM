% Black_hole1

clear 
clc
% MATLAB Simulation: Black Hole in Cellular Space

% Simulation Parameters
nx = 200; % Grid size (x-axis)
ny = 200; % Grid size (y-axis)
dx = 1;   % Spatial step size
dt = 0.1; % Time step size
C = 1;    % Constant of proportionality for wave speed

% Time parameters
t_end = 200; % Simulation end time
num_steps = round(t_end / dt);

% Initialize the cellular space
l0 = 1; % Initial uniform cell size
l = l0 * ones(nx, ny); % Cellular grid

% Seed a localized contraction
x_center = nx / 2; y_center = ny / 2; % Center of the contraction
radius = 20; % Radius of contraction
l_contracted = 0.1; % Minimum cell size (black hole core)

for i = 1:nx
    for j = 1:ny
        r = sqrt((i - x_center)^2 + (j - y_center)^2);
        if r <= radius
            l(i, j) = l_contracted + (l0 - l_contracted) * (r / radius); % Smooth contraction
        end
    end
end

% Wave speed (v = C * l)
v = C * l;

% Initialize wave field
u = zeros(nx, ny); % Wave amplitude
u_prev = u; % Previous time step
u_next = u; % Next time step

% Initialize a wave disturbance
disturbance_radius = 10;
for i = 1:nx
    for j = 1:ny
        r = sqrt((i - x_center - 30)^2 + (j - y_center)^2); % Offset disturbance
        if r <= disturbance_radius
            u(i, j) = exp(-r^2 / (2 * disturbance_radius^2));
        end
    end
end

% Simulation loop
figure(9);
for t = 1:num_steps
    % Update wave field using finite-difference wave equation
    for i = 2:nx-1
        for j = 2:ny-1
            laplacian = (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4 * u(i, j)) / dx^2;
            u_next(i, j) = 2 * u(i, j) - u_prev(i, j) + (dt^2 * v(i, j)^2) * laplacian;
        end
    end
    
    % Boundary conditions (reflective)
    u_next(1, :) = u_next(2, :);
    u_next(nx, :) = u_next(nx-1, :);
    u_next(:, 1) = u_next(:, 2);
    u_next(:, ny) = u_next(:, ny-1);
    
    % Update wave fields
    u_prev = u;
    u = u_next;
    
    % Visualization
    if mod(t, 10) == 0 % Plot every 10 time steps
        imagesc(u, [-1 1]); % Display wave amplitude
        colormap(jet);
        colorbar;
        title(['Time Step: ', num2str(t)]);
        xlabel('x');
        ylabel('y');
        axis equal tight;
        drawnow;
    end
end
