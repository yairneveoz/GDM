% gdm_warp_water_demo.m
% GDM warp-like field from a single H2O dipole (2D slice)
% H's -> contraction (+Δρ), O -> dilation (−Δρ)
% c'(x,y) = c*sqrt(rho0./rho),   g = -(1/2) c ∇c'

clear; 
close all; 
clc;

% -------- Constants / scales (normalized units) --------
rho0 = 1;      % undeformed space density
c    = 1;      % set c=1 for normalized units

% Geometry of water (Å) & orientation
r_OH      = 0.9572;        % O-H bond length
theta_HOH = 104.5;         % angle in degrees
phi_deg   = 0;             % rotate the molecule in-plane (0 = dipole along +x)

% Gaussian "charge" spreads (Å) and amplitudes
% + = contraction (H), − = dilation (O). Tune amplitudes to roughly neutralize.
sigma_H = 0.30;
sigma_O = 0.25;
A_H     = +0.12;           % each H contributes +A_H
A_O     = -0.24;           % O contributes A_O (≈ -2*A_H to keep near-neutral)

% -------- Grid (Å) --------
Lx=6; Ly=4; Nx=601; Ny=401;
x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(-Ly/2, Ly/2, Ny);
[XX,YY] = meshgrid(x,y);
dx = x(2)-x(1); dy = y(2)-y(1);

% -------- Place atoms (O at origin; H's on the bisector) --------
half = theta_HOH/2;
% Unrotated positions
pO  = [0; 0];
pH1 = [ r_OH*cosd(half);  r_OH*sind(half)];
pH2 = [ r_OH*cosd(half); -r_OH*sind(half)];
% Rotate molecule by phi_deg
R   = [cosd(phi_deg) -sind(phi_deg); sind(phi_deg) cosd(phi_deg)];
pO  = R*pO;  pH1 = R*pH1;  pH2 = R*pH2;

% -------- Space-density field  ρ = ρ0 * (1 + Σ A_i * Gaussian_i) --------
rO2  = (XX - pO(1)).^2  + (YY - pO(2)).^2;
rH12 = (XX - pH1(1)).^2 + (YY - pH1(2)).^2;
rH22 = (XX - pH2(1)).^2 + (YY - pH2(2)).^2;

rho = rho0 * ( 1 ...
    + A_H *exp(-rH12/(2*sigma_H^2)) ...
    + A_H *exp(-rH22/(2*sigma_H^2)) ...
    + A_O *exp(-rO2 /(2*sigma_O^2)) );

% Keep ρ positive (dilation is bounded in GDM)
rho = max(rho, 0.30*rho0);

% -------- Local light speed and g-field --------
cp = c*sqrt(rho0./rho);                  % c'(x,y)
[ dcp_dy, dcp_dx ] = gradient(cp, dy, dx);
Egx = -(0.5*c)*dcp_dx;                   % g = -(1/2) c ∇c'
Egy = -(0.5*c)*dcp_dy;

% g at the oxygen (molecular center)
[~,ix0] = min(abs(x-0)); [~,iy0] = min(abs(y-0));
g_center = [Egx(iy0,ix0), Egy(iy0,ix0)];
fprintf('g at O center = [%.3e, %.3e] (units: c^2/length)\n', g_center(1), g_center(2));

% -------- Plots --------
figure('Color','w','Position',[80 80 1200 760]);

subplot(2,2,1);
imagesc(x,y,(rho-rho0)/rho0); axis image; set(gca,'YDir','normal');
colorbar; title('Relative space-density  (\rho-\rho_0)/\rho_0'); hold on;
plot(pO(1), pO(2), 'wo', 'MarkerFaceColor','w', 'MarkerSize',6);
plot([pH1(1),pH2(1)],[pH1(2),pH2(2)],'w.','MarkerSize',18);
legend('','O','H','Location','southoutside'); legend boxoff;

subplot(2,2,2);
imagesc(x,y,cp); axis image; set(gca,'YDir','normal'); colorbar;
title('Local light speed  c''(x,y)'); hold on;
s = 12; xs = 1:s:Nx; ys = 1:s:Ny;
quiver(XX(ys,xs), YY(ys,xs), dcp_dx(ys,xs), dcp_dy(ys,xs), 'k');

subplot(2,2,3);
imagesc(x,y,hypot(Egx,Egy)); axis image; set(gca,'YDir','normal'); colorbar;
title('|g(x,y)| = (1/2) c | \nabla c'' |'); hold on;
s2 = 8; xs2 = 1:s2:Nx; ys2 = 1:s2:Ny;
quiver(XX(ys2,xs2), YY(ys2,xs2), Egx(ys2,xs2), Egy(ys2,xs2), 'k');

subplot(2,2,4);
% Plot along dipole axis (phi=0 uses y=0). For general phi, keep phi=0 for clarity.
iy_mid = find(abs(y-0) == min(abs(y-0)),1);
plot(x, cp(iy_mid,:), 'LineWidth',1.5); hold on; grid on;
yyaxis right; plot(x, Egx(iy_mid,:), 'LineWidth',1.5);
xlabel('x (dipole axis)'); legend('c''(x,0)','g_x(x,0)');
title('Centerline profiles');

sgtitle(sprintf('GDM water-dipole demo: A_H=%.2f, A_O=%.2f, \\sigma_H=%.2fÅ, \\sigma_O=%.2fÅ, \\phi=%g°', ...
    A_H, A_O, sigma_H, sigma_O, phi_deg));
