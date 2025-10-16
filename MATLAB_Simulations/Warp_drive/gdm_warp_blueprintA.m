function gdm_warp_blueprintA
% Strong asymmetric rho-pattern (Bow–Tail + Side-Cancelers) in 2D.
% Builds rho -> c'(x,y)=c*sqrt(rho0./rho) -> g = -(1/2) c grad c'.
% A tiny sweep picks a strong global amplitude that keeps rho >= rho_min*rho0.

clear; close all; clc;

% ---------- Base scales ----------
rho0    = 1;          % undeformed density
c_SI    = 3e8;        % for unit conversion note only
c       = 1;          % normalized c used in formulas
rho_min = 0.20;       % positivity clamp (fraction of rho0)

% ---------- Grid ----------
Lx=8; Ly=5; Nx=801; Ny=501;
x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(-Ly/2, Ly/2, Ny);
[XX,YY] = meshgrid(x,y);
dx = x(2)-x(1); dy = y(2)-y(1);

% ---------- Lobe configuration (Blueprint A) ----------
% Each row: [A, x0, y0, sx, sy, theta_deg]
cfg = [ ...
  % ---- Front bow (contraction) ----
   +0.22, +1.20,  0.00, 0.25, 0.18,   0;   % main bow (tight, strong)
   +0.12, +1.00, +0.45, 0.28, 0.20, -12;   % upper auxiliary
   +0.12, +1.00, -0.45, 0.28, 0.20, +12;   % lower auxiliary
  % ---- Rear tail (dilation, broad) ----
   -0.18, -1.25,  0.00, 0.65, 0.40,   0;   % broad core
   -0.07, -1.55, +0.55, 0.50, 0.30,  +8;   % upper wing
   -0.07, -1.55, -0.55, 0.50, 0.30,  -8;   % lower wing
  % ---- Side cancelers (small contraction) ----
   +0.06, +0.15, +0.95, 0.40, 0.25,   0;   % top stabilizer
   +0.06, +0.15, -0.95, 0.40, 0.25,   0];  % bottom stabilizer

% ---------- Tiny amplitude sweep (keeps rho positive) ----------
scales = linspace(0.7,1.6,10);
best = struct('score',-inf,'scale',NaN);   % <= keep 'score' field!

for s = scales
    [rho,~] = build_rho(cfg, s, XX,YY, rho0, rho_min);
    if min(rho,[],'all') < rho_min*rho0, continue; end

    cp = c*sqrt(rho0./rho);
    [dcp_dy,dcp_dx] = gradient(cp, dy, dx);
    Egx = -(0.5*c)*dcp_dx;  Egy = -(0.5*c)*dcp_dy;

    [~,ix0] = min(abs(x-0)); [~,iy0] = min(abs(y-0));
    rW = 0.50;  W = ((XX-x(ix0)).^2 + (YY-y(iy0)).^2) <= rW^2;
    gavg = mean(Egx(W),'all');                 % thrust proxy

    if gavg > best.score
        best.score = gavg; best.scale = s;
        best.rho = rho; best.cp = cp;
        best.Egx = Egx; best.Egy = Egy;
        best.ix0 = ix0; best.iy0 = iy0;
    end
end

if ~isfinite(best.score)
    error('All scales violated rho_min. Reduce amplitudes or raise rho_min.');
end

rho = best.rho; cp = best.cp; Egx = best.Egx; Egy = best.Egy;
ix0 = best.ix0; iy0 = best.iy0;

% ---------- Report ----------
g_center = [Egx(iy0,ix0), Egy(iy0,ix0)];
fprintf('Chosen amp_scale = %.3f (rho_min clamp = %.2f rho0)\n', best.scale, rho_min);
fprintf('Center g = [%.3e, %.3e] (normalized)\n', g_center(1), g_center(2));
fprintf('Area-avg g_x (r<=0.5) = %.3e (normalized)\n', best.score);
fprintf('To get SI-like units, multiply by (c_SI/c)^2 = %.3e.\n', (c_SI/c)^2);

% ---------- Plots ----------
figure('Color','w','Position',[60 60 1250 780]);

subplot(2,2,1);
imagesc(x,y,(rho-rho0)/rho0); axis image; set(gca,'YDir','normal');
colorbar; title('Relative space-density  (\rho-\rho_0)/\rho_0'); hold on; plot_craft();

subplot(2,2,2);
imagesc(x,y,cp); axis image; set(gca,'YDir','normal'); colorbar;
title('Local light speed  c''(x,y)'); hold on;
s = 12; xs=1:s:Nx; ys=1:s:Ny;
[ dcp_dy, dcp_dx ] = gradient(cp, dy, dx);
quiver(XX(ys,xs),YY(ys,xs),dcp_dx(ys,xs),dcp_dy(ys,xs),'k');

subplot(2,2,3);
imagesc(x,y,hypot(Egx,Egy)); axis image; set(gca,'YDir','normal'); colorbar;
title('|g(x,y)| = (1/2) c |\nabla c''|'); hold on;
s2 = 8; xs2=1:s2:Nx; ys2=1:s2:Ny;
quiver(XX(ys2,xs2),YY(ys2,xs2),Egx(ys2,xs2),Egy(ys2,xs2),'k'); plot_craft();

subplot(2,2,4);
plot(x, cp(iy0,:), 'LineWidth',1.5); hold on; grid on;
yyaxis right; plot(x, Egx(iy0,:), 'LineWidth',1.5);
xlabel('x (forward)'); legend('c''(x,0)','g_x(x,0)');
title('Centerline profiles (y=0)');

sgtitle(sprintf('Blueprint A (Bow–Tail+Cancelers)  |  amp\\_scale=%.2f  |  min(\\rho)/\\rho_0=%.2f', ...
    best.scale, min(rho,[],'all')/rho0));

% ---------- helpers ----------
function [rho, rho_factors] = build_rho(cfg, scale, XX,YY, rho0, rho_min)
    rho = rho0 * ones(size(XX));
    rho_factors = cell(size(cfg,1),1);
    for k = 1:size(cfg,1)
        A  = scale*cfg(k,1);
        x0 = cfg(k,2);  y0 = cfg(k,3);
        sx = cfg(k,4);  sy = cfg(k,5);
        th = deg2rad(cfg(k,6));
        ct = cos(th); st = sin(th);
        dxk = XX - x0; dyk = YY - y0;
        xpr =  ct*dxk + st*dyk;
        ypr = -st*dxk + ct*dyk;
        Gk  = 1 + A * exp(-0.5*((xpr./sx).^2 + (ypr./sy).^2));
        rho = rho .* Gk;
        rho_factors{k} = Gk; %#ok<NASGU>
    end
    rho = max(rho, rho_min*rho0); % bound from below
end

function plot_craft, plot(0,0,'wo','MarkerFaceColor','k','MarkerSize',6); end
end
