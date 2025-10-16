function proton_quark_paths_GDM()
% PROTON_QUARK_PATHS_GDM
% Builds and plots three smooth closed paths for a proton's uud quarks
% using torus-knot/toroidal-helix parametrizations, then saves the XYZ
% arrays to proton_paths.mat.
%
% Up quarks:   same (p,q), phase-shifted 120°.
% Down quark:  different (p,q) and opposite handedness.
%
% All curves are C^∞ smooth; r/R ~ alpha by default. No sharp corners.

%% ----- Config (edit here) -----
cfg.alphaFS   = 1/137.035999084;          % fine-structure
cfg.N         = 4000;                     % points per curve (smoothness)
cfg.R         = 5.0e-16;                  % major radius [m] (set your scale)
cfg.r         = cfg.alphaFS * cfg.R;      % minor radius [m]
cfg.centering = true;                      % recenter curves to origin

% Winding numbers (choose coprime p,q for closed torus knots)
% You can change these to your preferred GDM mapping.
% Example: up = T(2,3), down = T(1,3) with opposite handedness
cfg.up.p  = 2;  cfg.up.q  = 3;  cfg.up.handed  = +1;
cfg.down.p= 1;  cfg.down.q= 3;  cfg.down.handed= -1;

% Phase offsets to keep the two 'u' loops apart
phi_u1 = 0;                 % radians
phi_u2 = 2*pi/3;            % 120°
phi_d  = pi/3;              % optional offset for down

%% ----- Generate paths -----
[u1.X,u1.Y,u1.Z] = torus_knot(cfg.R,cfg.r,cfg.up.p,cfg.up.q,cfg.N,phi_u1,cfg.up.handed);
[u2.X,u2.Y,u2.Z] = torus_knot(cfg.R,cfg.r,cfg.up.p,cfg.up.q,cfg.N,phi_u2,cfg.up.handed);
[dq.X,dq.Y,dq.Z] = torus_knot(cfg.R,cfg.r,cfg.down.p,cfg.down.q,cfg.N,phi_d ,cfg.down.handed);

if cfg.centering
    % Recenter each curve independently (barycenter → 0)
    [u1.X,u1.Y,u1.Z] = recenter(u1.X,u1.Y,u1.Z);
    [u2.X,u2.Y,u2.Z] = recenter(u2.X,u2.Y,u2.Z);
    [dq.X,dq.Y,dq.Z] = recenter(dq.X,dq.Y,dq.Z);
end

%% ----- Plot -----
figure('Color','w'); hold on; grid on; axis equal
plot3(u1.X,u1.Y,u1.Z,'LineWidth',1.6); 
plot3(u2.X,u2.Y,u2.Z,'LineWidth',1.6);
plot3(dq.X,dq.Y,dq.Z,'LineWidth',1.6);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('Proton (uud) paths   r/R = %.5f,   R = %.2e m', cfg.r/cfg.R, cfg.R));
legend({'u_1 (T(2,3))','u_2 (T(2,3))','d (T(1,3), opposite handedness)'},'Location','best');

% Show the major circle (reference)
tt = linspace(0,2*pi,600);
plot3(cfg.R*cos(tt), cfg.R*sin(tt), 0*tt, '--', 'Color',[0.6 0.6 0.6]);

view(35,22)

%% ----- Save to MAT -----
proton_paths.u1 = u1; proton_paths.u2 = u2; proton_paths.d = dq; proton_paths.cfg = cfg;
save('proton_paths.mat','-struct','proton_paths');
fprintf('Saved paths to proton_paths.mat\n');

end % main

%% ----- Helpers -----
function [X,Y,Z] = torus_knot(R,r,p,q,N,phi,handed)
% TORUS_KNOT  x = (R + r cos(p t)) cos(q t)
%             y = (R + r cos(p t)) sin(q t)
%             z =  r sin(p t)
% p,q coprime integers give a closed knot; phase phi shifts along the torus.
% handed = +1 normal, -1 reverses parameter direction (opposite helicity).

    t = linspace(0,2*pi,N+1)'; t(end) = [];   % column vector, no duplicate end
    if handed<0, t = -t; end
    T = t + phi;

    cp = cos(p*T); sp = sin(p*T);
    cq = cos(q*T); sq = sin(q*T);

    X = (R + r.*cp).*cq;
    Y = (R + r.*cp).*sq;
    Z =  r.*sp;
end

function [Xc,Yc,Zc] = recenter(X,Y,Z)
% recenters by subtracting the barycenter (each as column vec)
    X = X(:); Y = Y(:); Z = Z(:);
    Xc = X - mean(X);  Yc = Y - mean(Y);  Zc = Z - mean(Z);
end
