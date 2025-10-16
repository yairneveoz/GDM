function gdm_warp_strong_demo
% Strong GDM warp configuration (2-D toy model)
% Bow (contraction) ahead + wide tail (dilation) behind + side cancelers.
% c'(x,y) = c*sqrt(rho0./rho),  g = -(1/2) c ∇c'

clear; close all; clc;

% -------- Base scales --------
rho0 = 1;      % undeformed space density
c    = 1;      % normalized
rho_min = 0.20;% clamp to keep rho positive

% -------- Grid --------
Lx=8; Ly=5; Nx=801; Ny=501;
x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(-Ly/2, Ly/2, Ny);
[XX,YY] = meshgrid(x,y);
dx = x(2)-x(1); dy = y(2)-y(1);

% -------- "Strong forward warp" lobe layout --------
% Design idea:
%  - Main front bow: 1 strong, 2 auxiliary contractive lobes -> steep ∂x c' ahead
%  - Rear tail: broad dilation (and two mild wings) -> keeps c' higher behind
%  - Side cancelers (contractive) near the flanks at x≈0 to kill lateral g
%
% Each row: [A, x0, y0, sx, sy, theta_deg]
%  A>0 contraction  (↑ρ → ↓c')
%  A<0 dilation     (↓ρ → ↑c')
cfg = [ ...
  % ---- Front bow (contraction) ----
   +0.14, +1.15,  0.00, 0.35, 0.22,   0;   % main bow
   +0.08, +0.95, +0.38, 0.32, 0.22, -12;   % upper bow
   +0.08, +0.95, -0.38, 0.32, 0.22, +12;   % lower bow
  % ---- Rear tail (dilation) ----
   -0.16, -1.10,  0.00, 0.55, 0.35,   0;   % broad core
   -0.06, -1.35, +0.45, 0.42, 0.26,  +8;   % upper wing
   -0.06, -1.35, -0.45, 0.42, 0.26,  -8;   % lower wing
  % ---- Side cancelers (contraction) ----
   +0.05,  0.10, +0.85, 0.40, 0.25,   0;   % top stabilizer
   +0.05,  0.10, -0.85, 0.40, 0.25,   0];  % bottom stabilizer

% Global amplitude scale (try 1.2–1.5, but keep rho positive)
amp_scale = 1.00;
cfg(:,1) = amp_scale * cfg(:,1);

% -------- Build rho from elliptical Gaussians --------
rho = rho0 * ones(Ny,Nx);
for k = 1:size(cfg,1)
    A  = cfg(k,1);
    x0 = cfg(k,2);  y0 = cfg(k,3);
    sx = cfg(k,4);  sy = cfg(k,5);
    th = deg2rad(cfg(k,6));
    ct = cos(th); st = sin(th);
    % rotate grid around (x0,y0): x' = c*(x-x0)+s*(y-y0); y' = -s*(x-x0)+c*(y-y0)
    dxk = XX - x0; dyk = YY - y0;
    xpr =  ct*dxk + st*dyk;
    ypr = -st*dxk + ct*dyk;
    rho = rho .* (1 + A * exp(-0.5*( (xpr./sx).^2 + (ypr./sy).^2 )));
end
rho = max(rho, rho_min*rho0);   % positivity clamp

% -------- Fields --------
cp = c*sqrt(rho0./rho);                 % local light speed
[dcp_dy, dcp_dx] = gradient(cp, dy, dx);
Egx = -(0.5*c)*dcp_dx;                  % g = -(1/2) c ∇c'
Egy = -(0.5*c)*dcp_dy;

% -------- Metrics (thrust proxies) --------
[~,ix0] = min(abs(x-0)); [~,iy0] = min(abs(y-0));
g_center = [Egx(iy0,ix0), Egy(iy0,ix0)];
% average g_x in a small disk around the craft
rW = 0.50;
W   = ((XX-x(ix0)).^2 + (YY-y(iy0)).^2) <= rW^2;
g_avg = sum(Egx(W),'all')/nnz(W);
fprintf('Center g = [%.3e, %.3e]  |  Area-avg g_x (r<=%.2f) = %.3e\n', ...
        g_center(1), g_center(2), rW, g_avg);

% -------- Plots --------
figure('Color','w','Position',[60 60 1250 780]);

subplot(2,2,1);
imagesc(x,y,(rho-rho0)/rho0); axis image; set(gca,'YDir','normal');
colorbar; title('Relative space-density  (\rho-\rho_0)/\rho_0'); hold on;
viscraft();

subplot(2,2,2);
imagesc(x,y,cp); axis image; set(gca,'YDir','normal'); colorbar;
title('Local light speed  c''(x,y)'); hold on;
s=12; xs=1:s:Nx; ys=1:s:Ny;
quiver(XX(ys,xs),YY(ys,xs),dcp_dx(ys,xs),dcp_dy(ys,xs),'k');

subplot(2,2,3);
imagesc(x,y,hypot(Egx,Egy)); axis image; set(gca,'YDir','normal'); colorbar;
title('|g(x,y)| = (1/2) c |\nabla c''|'); hold on;
s2=8; xs2=1:s2:Nx; ys2=1:s2:Ny;
quiver(XX(ys2,xs2),YY(ys2,xs2),Egx(ys2,xs2),Egy(ys2,xs2),'k');
viscraft();

subplot(2,2,4);
iy = iy0;
plot(x, cp(iy,:), 'LineWidth',1.5); hold on; grid on;
yyaxis right; plot(x, Egx(iy,:), 'LineWidth',1.5);
xlabel('x (forward)'); legend('c''(x,0)','g_x(x,0)');
title('Centerline profiles (y=0)');

sgtitle(sprintf('Strong 2-D GDM warp: amp_scale=%.2f, rho_{min}=%.2f', amp_scale, rho_min));

% ---- optional: flip thrust (reverse signs) ----
% cfg(:,1) = -cfg(:,1);  % then rebuild rho... (not done here)

% ---- nested utility to draw craft marker ----
    function viscraft
        plot(0,0,'wo','MarkerFaceColor','k','MarkerSize',6);
    end
end
