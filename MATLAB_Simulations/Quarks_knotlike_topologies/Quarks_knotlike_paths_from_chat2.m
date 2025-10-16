% GDM smooth 3-lobed closed path with constant speed (no sharp bends)
% Yair Neve-Oz — Aug 2025

clear; clc;

%% Geometry (toroidal 3-fold corrugation; all smooth sinusoids)
R = 1.0;      % major radius of the torus (overall size)
r = 0.32;     % minor corrugation amplitude (keep <= ~0.4*R for smoothness)
m = 1;        % base winding around the major circle (1 = single loop)
n = 3;        % 3 lobes -> "quark" regions
phi = 0.0;    % phase offset of the corrugation (rotate lobes)

% dense parameter for initial construction
N0 = 4000;
t0 = linspace(0, 2*pi, N0+1); t0(end) = [];

% smooth toroidal curve (C^∞)
x0 = (R + r*cos(n*t0 + phi)).*cos(m*t0);
y0 = (R + r*cos(n*t0 + phi)).*sin(m*t0);
z0 =  r*sin(n*t0 + phi);

% first and second derivatives (central differences on periodic curve)
wrap = @(v) [v(end), v, v(1)];        % helper for periodic padding
dt = t0(2)-t0(1);

dx = (wrap(x0(3:end)) - wrap(x0(1:end-2)))/(2*dt); dx = dx(2:end-1);
dy = (wrap(y0(3:end)) - wrap(y0(1:end-2)))/(2*dt); dy = dy(2:end-1);
dz = (wrap(z0(3:end)) - wrap(z0(1:end-2)))/(2*dt); dz = dz(2:end-1);

d2x = (wrap(x0(3:end)) - 2*wrap(x0(2:end-1)) + wrap(x0(1:end-2)))/(dt^2); d2x = d2x(2:end-1);
d2y = (wrap(y0(3:end)) - 2*wrap(y0(2:end-1)) + wrap(y0(1:end-2)))/(dt^2); d2y = d2y(2:end-1);
d2z = (wrap(z0(3:end)) - 2*wrap(z0(2:end-1)) + wrap(z0(1:end-2)))/(dt^2); d2z = d2z(2:end-1);

% speed and curvature (for diagnostics)
v  = sqrt(dx.^2 + dy.^2 + dz.^2);
cr = vecnorm(cross([dx;dy;dz].', [d2x;d2y;d2z].'), 2, 2) ./ (v.^3 + eps);
v = v(:);
% arc-length reparameterization for constant speed
s  = [0; cumsum( 0.5*(v(1:end-1)+v(2:end))*dt )];            % cumulative arc length
L  = s(end);
Ns = 6000;                                                    % final resolution
sU = linspace(0, L, Ns).';                                    % uniform arc-length grid

x = interp1(s, x0(2:end-1).', sU, 'pchip');
y = interp1(s, y0(2:end-1).', sU, 'pchip');
z = interp1(s, z0(2:end-1).', sU, 'pchip');

% re-compute phase along the path for smooth coloring of 3 regions
% use the original parameter as a proxy by mapping sU back to t0
t_of_s = interp1(s, t0(2:end-1).', sU, 'linear','extrap');
phase3 = mod(n*t_of_s + phi, 2*pi);   % [0, 2pi)

% soft region masks (Gaussian windows) for visualization, no hard edges
sigma = 0.45;  % width of window in radians (controls softness)
centers = (0:n-1)*(2*pi/n);
W = zeros(Ns, n);
for k = 1:n
    d = angle(exp(1i*(phase3 - centers(k))));  % shortest angular distance
    W(:,k) = exp(-0.5*(d/sigma).^2);
end
W = W ./ sum(W,2);  % normalize weights

% palette (distinct but soft)
cols = lines(n);

%% Plot the smooth constant-speed path
figure('Color','w'); hold on; axis equal vis3d; grid on;
% faint backbone
plot3(x, y, z, 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

% smoothly colored segments
for k = 1:n
    % draw with per-vertex alpha by splitting into small chunks
    idx = 1:Ns-1;
    for j = 1:50:Ns-1
        jj = j:min(j+49,Ns-1);
        plot3(x(jj), y(jj), z(jj), 'LineWidth', 3, 'Color', cols(k,:)*mean(W(jj,k)) + (1-mean(W(jj,k)))*[0.9 0.9 0.9]);
    end
end

% motion arrow (constant speed)
i0 = round(0.10*Ns); i1 = i0+20;
quiver3(x(i0), y(i0), z(i0), x(i1)-x(i0), y(i1)-y(i0), z(i1)-z(i0), 0, ...
        'k','LineWidth',1.2,'MaxHeadSize',1.5);

xlabel('x'); ylabel('y'); zlabel('z');
title('Smooth 3‑lobed closed path (constant speed, no kinks)'); view(35,24);

% annotate lobe “centroids”
for k = 1:n
    w = W(:,k);
    cx = sum(w.*x)/sum(w); cy = sum(w.*y)/sum(w); cz = sum(w.*z)/sum(w);
    if k==1
        txt = 'Lobe 1  (d-like, ~−1/3 e)';
    else
        txt = sprintf('Lobe %d  (u-part, ~+1/3 e)', k);
    end
    text(cx, cy, cz, ['  ' txt], 'FontSize', 10, 'Color', cols(k,:));
end

%% Curvature diagnostics (optional)
% show that curvature is finite and smooth (no spikes)
figure('Color','w');
plot(linspace(0,1,numel(cr)), cr, 'LineWidth',1.2);
xlabel('normalized parameter'); 
ylabel('\kappa (curvature)');
title('Curvature along the (pre-reparam) path — smooth, bounded'); grid on;

%% Time fractions per lobe with CONSTANT speed
% With equal arc-length parameterization, time ∝ length spent in each soft region.
frac = zeros(1,n);
for k = 1:n
    frac(k) = mean(W(:,k));
end
fprintf('Time fraction per lobe (should be ~1/3 each): [%.4f  %.4f  %.4f]\n', frac);
