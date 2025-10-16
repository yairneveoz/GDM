% Material properties
poisson_ratio = 0.1839;  % Poisson's ratio
density = 1.0;  % Density of the material
elastic_modulus = 1.0;  % Elastic modulus (arbitrary value for the simulation)

% Calculate wave velocities
v_longitudinal = sqrt((elastic_modulus * (1 - poisson_ratio)) / density);
v_transverse = sqrt(elastic_modulus / density);

% Simulation parameters
time = linspace(0, 2, 50);  % Time array
distance = linspace(0, 2, 50);  % Distance array

% Simulation of wave propagation
longitudinal_wave = sin(2*pi*v_longitudinal*time' - 2*pi*distance);
transverse_wave = sin(2*pi*v_transverse*time' - 2*pi*distance);

% Plotting the waves
figure(1);
clf

plot(distance, longitudinal_wave, 'LineWidth', 1.5);
hold on;
plot(distance, transverse_wave, 'LineWidth', 1.5);
xlabel('Distance');
ylabel('Amplitude');
title('Simulation of Wave Propagation');
legend('Longitudinal Wave', 'Transverse Wave');
grid on;
hold off