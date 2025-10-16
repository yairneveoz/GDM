function belt_trick_spinor_demo()
% belt_trick_spinor_demo
% Visualizes the 4π (spin-½) nature using a ribbon on a closed path.
% - Centerline: toroidal helix (closed, smooth).
% - Frame: Bishop (rotation-minimizing) along the curve.
% - Ribbon director uses half-angle: d(θ) = cos(θ/2) N + sin(θ/2) B.
% - The object is rotated by θ about z; d uses θ/2, so after 2π the ribbon flips,
%   after 4π it returns (belt trick).

%% Params
N      = 1500;     % samples along curve (smoothness)
R      = 4e-1;     % major radius (plot units)
alpha  = 1/137;    % just to set a small minor radius visually
r      = alpha*R;  % minor radius
n      = 3;        % poloidal windings of the helix (integer)
width  = 0.6*r;    % ribbon half-width
nFrames= 240;      % animation samples 0 -> 4π

%% Build centerline and Bishop frame
[t, X, Y, Z] = toroidalHelix(R, r, n, N);
C  = [X Y Z];
[Tvec, Nvec, Bvec] = bishopFrame(C);   % tangent, normal, binormal

%% Figure
f = figure('Color','w'); ax = axes('Parent',f); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z'); view(35,25);
title(ax,'Spin-1/2 belt trick: 0 \rightarrow 4\pi');

% Centerline (updated in-place)
hCenter = plot3(ax, C(:,1), C(:,2), C(:,3),'k-','LineWidth',1.5);

% Ribbon (surf with two long edges)
[Xrib, Yrib, Zrib] = ribbonFromDirector(C, Nvec, Bvec, width, 0); % θ=0
hSurf = surf(ax, Xrib, Yrib, Zrib, 'EdgeColor','none', 'FaceAlpha',0.95);
colormap(ax, turbo); % any default colormap is fine

% Two colored edges to show swap at 2π
hEdge1 = plot3(ax, Xrib(:,1), Yrib(:,1), Zrib(:,1), 'r-', 'LineWidth',1.4);
hEdge2 = plot3(ax, Xrib(:,2), Yrib(:,2), Zrib(:,2), 'b-', 'LineWidth',1.4);

% Major reference circle
tt = linspace(0,2*pi,400);
plot3(ax, R*cos(tt), R*sin(tt), 0*tt, '--', 'Color',[0.7 0.7 0.7]);

txt = text(ax, 0, 0, 1.2*max(Z)-0.1, '', 'HorizontalAlignment','center','FontWeight','bold');

%% Animate θ: 0 -> 4π
thetas = linspace(0, 4*pi, nFrames);
for k = 1:nFrames
    th = thetas(k);
    % Rotate centerline + frame by θ about z
    Cth   = rotZ(C, th);
    Tth   = rotZ(Tvec, th);
    Nth   = rotZ(Nvec, th);
    Bth   = rotZ(Bvec, th);

    % Build ribbon using HALF-ANGLE on the director
    [Xrb, Yrb, Zrb] = ribbonFromDirector(Cth, Nth, Bth, width, th);

    % Update graphics
    set(hCenter, 'XData', Cth(:,1), 'YData', Cth(:,2), 'ZData', Cth(:,3));
    set(hSurf,   'XData', Xrb, 'YData', Yrb, 'ZData', Zrb);
    set(hEdge1,  'XData', Xrb(:,1), 'YData', Yrb(:,1), 'ZData', Zrb(:,1));
    set(hEdge2,  'XData', Xrb(:,2), 'YData', Yrb(:,2), 'ZData', Zrb(:,2));

    % Label: show θ and spinor phase θ/2
    set(txt, 'String', sprintf('\\theta = %3.0f^\\circ   (spinor phase = %3.0f^\\circ)', ...
        th*180/pi, (th/2)*180/pi));

    drawnow;
end

% Helpful markers:
disp('Note: At θ = 360° the centerline is back, but red/blue edges are swapped (not identical).');
disp('Only at θ = 720° do the edges return—visual 4π periodicity (spin-½).');

end % main

%% --------- Helpers ---------

function [t, X, Y, Z] = toroidalHelix(R, r, n, N)
% Toroidal helix (closed spiral around a circle)
% x = (R + r cos(n t)) cos t
% y = (R + r cos(n t)) sin t
% z =  r sin(n t)
    t = linspace(0, 2*pi, N+1)'; t(end) = [];
    ct = cos(t); st = sin(t);
    cn = cos(n*t); sn = sin(n*t);
    X = (R + r.*cn).*ct;
    Y = (R + r.*cn).*st;
    Z =  r.*sn;
end

function [T, N, B] = bishopFrame(C)
% Rotation-minimizing (Bishop) frame along curve C (Nx3).
    Npts = size(C,1);
    T = zeros(Npts,3);
    N = zeros(Npts,3);
    B = zeros(Npts,3);

    % Tangents
    dC = gradient(C);
    for i=1:Npts
        T(i,:) = safeNormalize(dC(i,:));
    end

    % Initial normal: pick something not parallel to T(1,:)
    v = [1,0,0];
    if abs(dot(v,T(1,:))) > 0.9, v = [0,1,0]; end
    N(1,:) = safeNormalize(v - dot(v,T(1,:))*T(1,:));
    B(1,:) = cross(T(1,:), N(1,:));

    % Discrete parallel transport
    for i=1:Npts-1
        t1 = T(i,:); t2 = T(i+1,:);
        axis = cross(t1, t2);
        s = norm(axis);
        c = dot(t1, t2);
        if s < 1e-12
            % Nearly same tangent -> keep normal
            N(i+1,:) = N(i,:);
        else
            axis = axis / s;
            ang  = atan2(s, c);
            R = rotFromAxisAngle(axis, ang);
            N(i+1,:) = (R * N(i,:).').';
        end
        % Orthonormalize
        N(i+1,:) = safeNormalize(N(i+1,:) - dot(N(i+1,:),t2)*t2);
        B(i+1,:) = cross(t2, N(i+1,:));
    end
end

function v = safeNormalize(v)
    n = norm(v);
    if n < 1e-15, v = [1,0,0]; else, v = v/n; end
end

function R = rotFromAxisAngle(axis, ang)
% 3x3 rotation from unit axis, angle (Rodrigues)
    x=axis(1); y=axis(2); z=axis(3);
    c=cos(ang); s=sin(ang); C=1-c;
    R = [ x*x*C+c,   x*y*C - z*s, x*z*C + y*s;
          y*x*C + z*s, y*y*C+c,   y*z*C - x*s;
          z*x*C - y*s, z*y*C + x*s, z*z*C + c ];
end

function Pout = rotZ(P, theta)
% Rotate Nx3 (or 1x3) points/vectors about z by angle theta
    Rz = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    if isvector(P)
        Pout = (Rz * P(:)).';
    else
        Pout = (Rz * P.').';
    end
end

function [Xrib, Yrib, Zrib] = ribbonFromDirector(C, Nvec, Bvec, w, theta)
% Build ribbon surface using d(θ) = cos(θ/2) N + sin(θ/2) B
% Returns two-column matrices for surf().
    c = cos(theta/2); s = sin(theta/2);
    D  = c.*Nvec + s.*Bvec;  % Nx3 director at half-angle
    E1 = C + w*D;            % one edge
    E2 = C - w*D;            % opposite edge (flips at 2π)

    Xrib = [E1(:,1) E2(:,1)];
    Yrib = [E1(:,2) E2(:,2)];
    Zrib = [E1(:,3) E2(:,3)];
end
