% GDM quark paths across generations (1,2,3)
% Smooth 3‑lobe closed loops, arc‑length reparam (constant speed), and diagnostics
% Yair Neve‑Oz — Aug 2025

clear; clc;

%% Global knobs (same family; generations change mode + scale)
R0   = 1.00;     % base major radius (gen-1)
a0   = 0.32;     % base 3-lobe corrugation amplitude (gen-1)
b0   = 0.28;     % base out-of-plane amplitude (gen-1)
beta = 0.75;     % radius scale per generation: R_g = R0*beta^(g-1)
eps0 = 0.15;     % phase-modulation depth per mode unit (twist smoothness)
lobes= 3;        % keep 3 lobes for all generations

N0   = 40000;    % dense pre-sampling for smooth derivatives
Ns   = 6000;     % points after arc-length reparam (final uniform speed)

gens = 1:3;      % generations
cols = lines(numel(gens));

% Storage for a small summary table
Summary = struct('Gen',[],'Twist',[],'R',[],'a',[],'b',[],'eps',[], ...
                 'Length',[],'KappaMean',[],'KappaMax',[]);

figure('Color','w','Name','Quark Paths by Generation'); 
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for idx = 1:numel(gens)
    g = gens(idx);                   % generation index (1,2,3)
    twist = g;                       % "degree of twist" (mode index)
    R = R0*beta^(g-1);
    a = a0*beta^(g-1);
    b = b0*beta^(g-1);
    eps = eps0*twist;

    % -------- Smooth parametric path (macro 3-lobe with micro twist) --------
    % Use phase modulation to add twist without creating kinks or straight bits.
    t0 = linspace(0,2*pi,N0+1); t0(end) = [];
    phase = lobes*t0 + eps*sin(twist*t0);      % retains 3-lobe macro pattern
    x0 = (R + a*cos(phase)).*cos(t0);
    y0 = (R + a*cos(phase)).*sin(t0);
    z0 =  b*sin(phase);

    % -------- Derivatives (periodic central difference) --------
    wrap = @(v)[v(end), v, v(1)];
    dt = t0(2)-t0(1);
    dx  = (wrap(x0(3:end)) - wrap(x0(1:end-2)))/(2*dt); dx  = dx(2:end-1);
    dy  = (wrap(y0(3:end)) - wrap(y0(1:end-2)))/(2*dt); dy  = dy(2:end-1);
    dz  = (wrap(z0(3:end)) - wrap(z0(1:end-2)))/(2*dt); dz  = dz(2:end-1);
    d2x = (wrap(x0(3:end)) - 2*wrap(x0(2:end-1)) + wrap(x0(1:end-2)))/(dt^2); d2x = d2x(2:end-1);
    d2y = (wrap(y0(3:end)) - 2*wrap(y0(2:end-1)) + wrap(y0(1:end-2)))/(dt^2); d2y = d2y(2:end-1);
    d2z = (wrap(z0(3:end)) - 2*wrap(z0(2:end-1)) + wrap(z0(1:end-2)))/(dt^2); d2z = d2z(2:end-1);

    v   = sqrt(dx.^2 + dy.^2 + dz.^2);
    cr0 = vecnorm(cross([dx;dy;dz].', [d2x;d2y;d2z].'),2,2) ./ (v.^3 + eps);
    v = v(:);
    % -------- Arc-length reparam (uniform speed) --------
    s   = [0; cumsum( 0.5*(v(1:end-1)+v(2:end))*dt )];
    L   = s(end);
    sU  = linspace(0,L,Ns).';

    x = interp1(s, x0(2:end-1).', sU, 'pchip');
    y = interp1(s, y0(2:end-1).', sU, 'pchip');
    z = interp1(s, z0(2:end-1).', sU, 'pchip');

    % Recompute curvature on uniform-arc path
    ds = sU(2)-sU(1);
    Dx = gradient(x, ds);  Dy = gradient(y, ds);  Dz = gradient(z, ds);
    D2x= gradient(Dx,ds);  D2y= gradient(Dy,ds);  D2z= gradient(Dz,ds);
    vU = sqrt(Dx.^2 + Dy.^2 + Dz.^2);
    cr = vecnorm(cross([Dx Dy Dz], [D2x D2y D2z], 2), 2, 2) ./ (vU.^3 + eps);

    % -------- Plot 3D path --------
    nexttile; hold on; axis equal vis3d; grid on;
    plot3(x, y, z, 'Color', cols(idx,:), 'LineWidth', 2);
    % draw direction arrow
    i0 = round(0.10*Ns); i1 = i0+25;
    quiver3(x(i0),y(i0),z(i0), x(i1)-x(i0),y(i1)-y(i0),z(i1)-z(i0), 0, ...
            'k','LineWidth',1.2,'MaxHeadSize',1.5);
    title(sprintf('Gen %d  (twist=%d)', g, twist));
    xlabel x; ylabel y; zlabel z; view(35,24);

    % -------- Store summary --------
    Summary(idx).Gen        = g;
    Summary(idx).Twist      = twist;
    Summary(idx).R          = R;
    Summary(idx).a          = a;
    Summary(idx).b          = b;
    Summary(idx).eps        = eps;
    Summary(idx).Length     = L;
    Summary(idx).KappaMean  = mean(cr);
    Summary(idx).KappaMax   = max(cr);
    
    % For later overlay plot
    CR{idx}  = cr; %#ok<SAGROW>
end

% -------- Curvature overlay (shows increasing curvature with generation) --------
figure('Color','w','Name','Curvature vs arclength');
hold on; grid on;
for idx = 1:numel(gens)
    s_norm = linspace(0,1,numel(CR{idx}));
    plot(s_norm, CR{idx}, 'LineWidth', 1.5, 'Color', cols(idx,:));
end
xlabel('normalized arclength s/L'); ylabel('\kappa(s)');
title('Curvature along the constant-speed path');
legend(arrayfun(@(g)sprintf('Gen %d',g), gens, 'UniformOutput',false), 'Location','best');

% -------- Print parameter table --------
fprintf('\n%-6s %-6s %-8s %-8s %-8s %-8s %-10s %-12s %-12s\n', ...
        'Gen','Twist','R','a','b','eps','Length','KappaMean','KappaMax');
for i = 1:numel(Summary)
    S = Summary(i);
    fprintf('%-6d %-6d %-8.3f %-8.3f %-8.3f %-8.3f %-10.4f %-12.4f %-12.4f\n', ...
        S.Gen, S.Twist, S.R, S.a, S.b, S.eps, S.Length, S.KappaMean, S.KappaMax);
end

%% Notes
% • All three generations share the same 3‑lobe macro‑geometry (quark charge
%   fractions), while “generation” is encoded as:
%     – Twist index (mode) = 1,2,3 via smooth phase modulation (eps*sin(twist*t))
%     – Geometric compaction: R,a,b scaled by beta^(g-1) (higher gen = tighter)
% • The arc‑length reparam enforces constant speed (no straight bits / no kinks).
% • Tweak beta (size scaling) and eps0 (twist depth) to emphasize differences.
