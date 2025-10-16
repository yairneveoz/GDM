% Solitons1

clear
clc
% T = 200;

% Parameters
L = 100;           % Domain size (NxNxN grid)
Nx = 100;          % Grid resolution in x
Ny = 100;          % Grid resolution in y
Nz = 100;          % Grid resolution in z
dx = L / Nx;       % Grid spacing
dy = L / Ny;
dz = L / Nz;
dt = 0.1;          % Time step
T = 200;           % Total simulation time
c0 = 1;            % Scaling constant for wave speed

% Initialize the spatially varying cell size l(x, y, z)
[X, Y, Z] = ndgrid(linspace(-L/2, L/2, Nx), ...
                   linspace(-L/2, L/2, Ny), ...
                   linspace(-L/2, L/2, Nz));
l = 1 + 5 * exp(-0.01 * (X.^2 + Y.^2 + Z.^2)); % Gaussian cell size variation

% Compute wave speed as v(x, y, z) = c0 * l(x, y, z)
v = c0 * l;

% Initialize wave fields
u = zeros(Nx, Ny, Nz);      % Wave amplitude at current time step
u_prev = zeros(Nx, Ny, Nz); % Wave amplitude at previous time step
u_next = zeros(Nx, Ny, Nz); % Wave amplitude at next time step

% Add initial wave disturbance (e.g., Gaussian pulse)
x0 = 0; y0 = 0; z0 = 0; % Pulse center
sigma = 5; % Width of initial pulse
u = exp(-((X-x0).^2 + (Y-y0).^2 + (Z-z0).^2) / (2 * sigma^2));

gamma = 1; % Nonlinearity strength

figure(10)

% Time evolution
for t = 1:T
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 2:Nz-1
                % Laplacian term
                laplacian = (u(i+1,j,k) - 2*u(i,j,k) + u(i-1,j,k)) / dx^2 + ...
                            (u(i,j+1,k) - 2*u(i,j,k) + u(i,j-1,k)) / dy^2 + ...
                            (u(i,j,k+1) - 2*u(i,j,k) + u(i,j,k-1)) / dz^2;
                
                % Update wave equation with nonlinearity
                u_next(i,j,k) = 2*u(i,j,k) - u_prev(i,j,k) + ...
                                (dt^2 * v(i,j,k)^2) * laplacian + gamma * dt^2 * u(i,j,k)^3;
            end
        end
    end
    
    % Apply boundary conditions (e.g., zero at edges)
    u_next(1,:,:) = 0; u_next(end,:,:) = 0;
    u_next(:,1,:) = 0; u_next(:,end,:) = 0;
    u_next(:,:,1) = 0; u_next(:,:,end) = 0;

    % Update fields for the next iteration
    u_prev = u;
    u = u_next;

    % Visualization every few time steps
    if mod(t, 10) == 0
        sliceX = squeeze(u(:, round(Ny/2), :));
        imagesc(sliceX);
        colormap(jet); colorbar;
        title(['Wave Propagation at Time Step ', num2str(t)]);
        xlabel('Z'); ylabel('X');
        pause(0.1);
    end
end
