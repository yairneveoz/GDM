% ricci_like_curve_flow.m
% Curve-shortening flow (a boundary analogue of 2D Ricci flow)
% Starts from a square and evolves smoothly toward a circle.

clear; 
clc; 
close all;

% --------------------------- Parameters -----------------------------------
N        = 400;       % number of curve sample points (more = smoother)
T_final  = 0.30;      % total "flow time"
dt       = 1e-5;      % time step (keep small for stability)
redistrib_every = 20; % reparameterize by arclength every k steps
plot_every      = 10; % plot every k steps

% ------------------------ Initial square curve ----------------------------
L  = 1;                 % half-length of square side
nSide = N/4;            % assumes N divisible by 4
x0 = [linspace(-L, L, nSide),  L*ones(1,nSide), linspace( L,-L, nSide), -L*ones(1,nSide)];
y0 = [ -L*ones(1,nSide), linspace(-L, L, nSide),  L*ones(1,nSide), linspace( L,-L, nSide)];
x = x0(:); 
y = y0(:);

% Slight pre-smoothing to help stability at sharp corners (optional)
[x,y] = smooth_closed_curve(x,y,3);

% -------------------------- Helpers ---------------------------------------
idx = @(k) mod(k-1,N)+1;   % periodic indexing 1..N
t = 0;

% Figure setup
figure('Color','w'); 
axis equal off;
title_handle = title('Curve-shortening flow (square \rightarrow circle)');
colormap(parula);

% -------------------------- Time stepping ---------------------------------
step = 0;
while t < T_final
    step = step + 1;

    % Periodic neighbors
    xm = x(idx((1:N)-1));  xp = x(idx((1:N)+1));
    ym = y(idx((1:N)-1));  yp = y(idx((1:N)+1));

    % Tangents and arclength element
    tx = xp - xm;
    ty = yp - ym;
    ds = hypot(tx,ty);
    % Avoid zeros
    ds(ds==0) = mean(ds(ds>0));

    % Unit tangent
    tx = tx ./ ds;
    ty = ty ./ ds;

    % Central difference of unit tangent along the curve (approx dT/ds)
    Txm = tx(idx((1:N)-1)); Tym = ty(idx((1:N)-1));
    Txp = tx(idx((1:N)+1)); Typ = ty(idx((1:N)+1));

    % Effective ds for central difference (average of neighboring segments)
    ds_f = 0.5*(ds + ds(idx((1:N)+1)));
    ds_b = 0.5*(ds + ds(idx((1:N)-1)));
    dTds_x = (Txp - Txm) ./ (0.5*(ds_f + ds_b));
    dTds_y = (Typ - Tym) ./ (0.5*(ds_f + ds_b));

    % Unit normal = rotate tangent by +90°: n = (-ty, tx)
    nx = -ty; ny = tx;

    % Signed curvature: kappa = dT/ds · n
    kappa = dTds_x .* nx + dTds_y .* ny;

    % Explicit Euler step: x_t = kappa * n  (curve-shortening flow)
    x = x + dt * (kappa .* nx);
    y = y + dt * (kappa .* ny);
    t = t + dt;

    % Optional tangential redistribution to keep points well spaced
    if mod(step, redistrib_every) == 0
        [x,y] = redistribute_arclength(x,y);
        % tiny smoothing to reduce numerical noise
        [x,y] = smooth_closed_curve(x,y,1);
    end

    % Plot
    if mod(step, plot_every) == 0 || t >= T_final
        cla;
        % Color by curvature magnitude (normalized)
        c = rescale(abs(kappa));
        patch('XData',x,'YData',y,'FaceColor','none','EdgeColor','flat','CData',c,'LineWidth',1.5);
        hold on;
        plot([x; x(1)], [y; y(1)], 'k-', 'LineWidth', 1);
        axis equal; axis([-1.3 1.3 -1.3 1.3]);
        set(gca,'XTick',[],'YTick',[],'Box','on');
        title_handle.String = sprintf('Curve-shortening flow  t = %.5f', t);
        drawnow;
    end
end

disp('Done.');

% ======================== Local functions =================================
function [x2,y2] = smooth_closed_curve(x,y,repeat)
    % Very light Laplacian smoothing on a closed polygonal curve
    if nargin<3, repeat=1; end
    N = numel(x); idx = @(k) mod(k-1,N)+1;
    x2 = x; y2 = y;
    w  = 0.1; % smoothing weight (small!)
    for r = 1:repeat
        xm = x2(idx((1:N)-1)); xp = x2(idx((1:N)+1));
        ym = y2(idx((1:N)-1)); yp = y2(idx((1:N)+1));
        x2 = (1-2*w)*x2 + w*(xm+xp);
        y2 = (1-2*w)*y2 + w*(ym+yp);
    end
end

function [x2,y2] = redistribute_arclength(x,y)
    % Reparameterize points to be (approximately) equally spaced by arclength
    N      = numel(x);
    idx    = @(k) mod(k-1,N)+1;

    % Segment lengths (including wrap-around)
    dx     = x(idx((1:N)+1)) - x;
    dy     = y(idx((1:N)+1)) - y;
    seglen = hypot(dx,dy);             % length N, last is wrap segment

    % Cumulative arclength with the 0 at the start and full L at the end
    s      = [0; cumsum(seglen)];      % length N+1
    L      = s(end);

    % Periodic coordinate arrays of length N+1
    xp     = [x; x(1)];
    yp     = [y; y(1)];

    % Uniform arclength targets (N+1 points, last equals L)
    s_uni  = linspace(0, L, N+1).';

    % Interpolate, then drop the duplicated last point
    x2p    = interp1(s, xp, s_uni, 'linear');
    y2p    = interp1(s, yp, s_uni, 'linear');

    x2     = x2p(1:end-1);
    y2     = y2p(1:end-1);
end
