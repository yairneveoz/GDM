function [r_path] = getKnotLikePaths(v,loops,R0,r0,n,phase0)

L0 = 2*pi*R0;

cT = 1;
% cL = 1.6068*cT;
t = 0:0.02:loops*2*pi;
gamma = 1/sqrt(1 - v^2/cT^2);

R = R0/gamma;
r = r0/gamma;


X = R*cos(t);
Y = R*sin(t);
% Z = R*v*t*gamma;
Z = L0*(v/cT)*t;

% R_path = [X',Y',Z'];

x = X + r.*cos(n*t + 0).*cos(t + phase0);
y = Y + r.*cos(n*t + 0).*sin(t + phase0);
z = Z + r.*sin(n*t + 0);

r_path = [x',y',z'];

end