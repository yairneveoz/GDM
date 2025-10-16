+9+% GDM_1D_Dynamics_Simulation
 +*9-+
clear
clc
 
% Spatial grid
x = linspace(0, 1e25, 500);  % cm

% Mass distribution: Gaussian at center
mass_center = 5e24;          % cm
mass_width = 5e23;
mass_distribution = exp(-((x - mass_center) / mass_width).^2);

% Mass gradient (as a proxy for attractive force)
dx = x(2) - x(1);
mass_gradient = -gradient(mass_distribution, dx);

% Time settings
dt = 1e14;        % seconds
num_steps = 200;

% Particles: Initial positions and velocities
num_particles = 100;
particle_positions = linspace(0, 1e25, num_particles)';
particle_velocities = zeros(num_particles, 1);

% Record positions over time
positions_over_time = zeros(num_particles, num_steps);
positions_over_time(:,1) = particle_positions;

% Repulsion field function
repulsion_field = @(pos) (pos - mass_center) / 1e25;  % outward force

% Simulation loop
for t = 2:num_steps
    % Interpolate attraction and repulsion at particle positions
    attr_force = interp1(x, mass_gradient, particle_positions, 'linear', 'extrap');
    repul_force = repulsion_field(particle_positions);
    
    % Net acceleration
    acceleration = attr_force + repul_force;
    
    % Update velocity and position
    particle_velocities = particle_velocities + acceleration * dt;
    particle_positions = particle_positions + particle_velocities * dt;
    
    % Store
    positions_over_time(:,t) = particle_positions;
end

% Plot animation
figure(1);
for t = 1:num_steps
    plot(linspace(0, 1e25, num_particles), positions_over_time(:,t), 'bo');
%     xlim([0, 1e25]);
%     ylim([0, 1e25]);
    xlabel('Initial Position (cm)');
    ylabel('Current Position (cm)');
    title(sprintf('1D GDM Simulation, Step %d', t));
    pause(0.1)
    drawnow;
end
