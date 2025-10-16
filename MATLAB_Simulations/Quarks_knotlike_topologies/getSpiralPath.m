function [R_path, r_path] = getSpiralPath(v,loops)

alpha = 1/137.0366;


r0 = 1.0; 
R0 = r0/alpha; 
L0 = 2*pi*R0;

n = 17.34; % 20; r cycles in R cycles

cT = 1;
cL = 1.6068*cT;
t = 0:0.005:loops*2*pi;
gamma = 1/sqrt(1 - v^2/cT^2);

R = R0/gamma;
r = r0/gamma;


X = R*cos(t);
Y = R*sin(t);
% Z = R*v*t*gamma;
Z = L0*(v/cT)*t;

R_path = [X',Y',Z'];

x = X + r.*cos(n*t).*cos(t);
y = Y + r.*cos(n*t).*sin(t);
z = Z + r.*sin(n*t);

r_path = [x',y',z'];
end