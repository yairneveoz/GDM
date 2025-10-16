function [r_path] = getKnotLikePathsCLConst(v,loops,R0,r0,n,phase0)

L0 = 2*pi*R0;

cT = 1; % Transverse light speed
cL = 1.6068*cT; % Longitudinal light speed

t0 = 0:0.05:loops*2*pi;
gamma = 1/sqrt(1 - v^2/cT^2);

R = R0/gamma;
r = r0/gamma;


X = R*cos(t0);
Y = R*sin(t0);
Z = L0*(v/cT)*t0;

% R_path = [X',Y',Z'];
% t(R)
x0 = X + r.*cos(n*t0 + 0).*cos(t0 + phase0);
y0 = Y + r.*cos(n*t0 + 0).*sin(t0 + phase0);
z0 = Z + r.*sin(n*t0 + 0);

r0_path = [x0',y0',z0'];

%%%
[diff_path, diff_length] = getDiffPath(r0_path);
mean_diff_length = mean(diff_length);

%%%
end