% GDM quark-generating closed path (trefoil-like torus knot)
% A single continuous loop split into 3 lobes (sub-tracks).
% One lobe ≈ d (−1/3 e); two lobes together ≈ u (+2/3 e).

clear; 
clc;

% ----- Geometry (adjust these) -----
R = 1.0;        % major radius of the torus (overall size)
r = 0.35;       % minor "twist" amplitude (lobe prominence)
m = 2;          % torus-knot windings around major radius
n = 3;          % number of lobes (kept at 3 for quarks)

% ----- Parametric path (trefoil-like torus knot) -----
N = 6000;                           % resolution
t = linspace(0, 2*pi, N); t(end) = [];   % avoid duplicate endpoint

x = (R + r*cos(n*t)).*cos(m*t);
y = (R + r*cos(n*t)).*sin(m*t);
z =  r*sin(n*t);

% ----- Identify lobe membership (by phase of n*t) -----
phi  = mod(n*t, 2*pi);
bins = 0:(2*pi/n):2*pi;
lobe = discretize(phi, bins);   % 1..n

% ----- Plot -----
figure('Color','w'); hold on; grid on; axis equal
% Faint full curve for continuity
plot3(x, y, z, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);

C = lines(n);
for k = 1:n
    idx = (lobe == k);
    plot3(x(idx), y(idx), z(idx), 'LineWidth', 2.5, 'Color', C(k,:));
end

% Direction arrow
i0 = round(N/12); i1 = i0 + 12;
quiver3(x(i0), y(i0), z(i0), x(i1)-x(i0), y(i1)-y(i0), z(i1)-z(i0), ...
        0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 1.5);

% Labels at lobe centroids
for k = 1:n
    idx = (lobe == k);
    cx = mean(x(idx)); cy = mean(y(idx)); cz = mean(z(idx));
    if k == 1
        text(cx, cy, cz, '  Lobe 1: d-like (−1/3 e)', 'FontSize', 10, 'Color', C(k,:));
    else
        text(cx, cy, cz, sprintf('  Lobe %d: u-part (+1/3 e)', k), 'FontSize', 10, 'Color', C(k,:));
    end
end

view(35, 24);
xlabel('x'); ylabel('y'); zlabel('z');
title('GDM Twisted Charge Path (3-lobed closed loop)');

% ----- Report time fractions per lobe (uniform speed along t) -----
frac = arrayfun(@(k) mean(lobe == k), 1:n);
fprintf('Time fraction per lobe (should be ~1/3 each): [%.3f %.3f %.3f]\n', frac);

% ----- Notes -----
% • Increase r for “sharper” lobes; decrease for a rounder path.
% • Keep n = 3 for quark picture. With uniform circulation,
%   one lobe gives −1/3 e (d-like), the other two together +2/3 e (u-like).
