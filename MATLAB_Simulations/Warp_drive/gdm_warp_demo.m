% gdm_warp_demo.m
% GDM warp-bubble conditions (2-D toy model, no helper functions needed)
% Asymmetric space-density: contraction ahead (+Δρ), dilation behind (−Δρ)
% c'(x,y) = c*sqrt(rho0./rho),   g = -(1/2) c ∇c'

clear; close all; clc;

% ---------- Units & constants (normalized) ----------
rho0 = 1;      % undeformed space density
c    = 1;      % speed of light scale

% ---------- Grid ----------
Lx=6; Ly=4; Nx=601; Ny=401;
x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(-Ly/2, Ly/2, Ny);
[XX,YY] = meshgrid(x,y);
dx = x(2)-x(1); dy = y(2)-y(1);

% ---------- Bubble layout ----------
R     = 0.8;          % half-separation of lobes ("bubble radius")
sigma = 0.20;         % wall thickness
A_front = +0.06;      % + contraction ahead  (↑ρ → ↓c')
A_back  = -0.06;      % - dilation   behind (↓ρ → ↑c')
x_front = +R;  x_back = -R;  y0 = 0;

% ---------- Space density field ----------
rho = rho0*(1 ...
    + A_front*exp(-((XX-x_front).^2 + (YY-y0).^2)/(2*sigma^2)) ...
    + A_back *exp(-((XX-x_back ).^2 + (YY-y0).^2)/(2*sigma^2)));
rho = max(rho, 0.05*rho0);   % keep positive

% ---------- Local light speed and g-field ----------
cp = c*sqrt(rho0./rho);               % c'(x,y) per GDM
[ dcp_dy, dcp_dx ] = gradient(cp, dy, dx);  % built-in; sizes match cp

Egx = -(0.5*c)*dcp_dx;                % g = -(1/2) c ∇c'
Egy = -(0.5*c)*dcp_dy;

% Value at bubble center
[~,ix0] = min(abs(x-0)); [~,iy0] = min(abs(y-0));
g_center = [Egx(iy0,ix0), Egy(iy0,ix0)];
fprintf('Center acceleration g = [%.3e, %.3e] (units: c^2/length)\n', g_center(1), g_center(2));

% ---------- Plots ----------
figure('Color','w','Position',[80 80 1200 760]);

subplot(2,2,1);
imagesc(x,y,(rho-rho0)/rho0); axis image; set(gca,'YDir','normal');
colorbar; title('Relative space-density  (\rho-\rho_0)/\rho_0'); hold on;
th = linspace(0,2*pi,360); plot(R*cos(th), R*sin(th),'w--','LineWidth',0.8);

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
plot(x, cp(iy0,:), 'LineWidth',1.5); hold on; grid on;
yyaxis right; plot(x, Egx(iy0,:), 'LineWidth',1.5);
xlabel('x (forward)'); legend('c''(x,0)','g_x(x,0)');
title('Centerline profiles');

sgtitle(sprintf('GDM warp-bubble demo: A_f=%.3f (front contract), A_b=%.3f (rear dilate), R=%.2f, \\sigma=%.2f', ...
    A_front, A_back, R, sigma));

% --------- Try reversing thrust: swap signs ----------
% A_front = -A_front; A_back = -A_back;  % then re-run the block above
