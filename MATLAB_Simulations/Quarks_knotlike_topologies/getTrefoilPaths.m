function [r_path] = getTrefoilPaths(v,loops,R0,nr,nz,phase0)

L0 = 2*pi*R0;

cT = 1;
% cL = 1.6068*cT;
t = 0:0.05:loops*2*pi;
gamma = 1/sqrt(1 - v^2/cT^2);

R = R0/gamma;
% r = r0/gamma;

x = R*(cos(t) - nr*cos(nr*t));
y = R*(sin(t) + nr*sin(nr*t));
z = -nz*R*sin(nz*t);

r_path = [x',y',z'];

figure(3)
subplot(1,2,1)
plot3(x', y', z', 'r-','LineWidth',2)
grid on
axis equal

subplot(1,2,2)
plot3(x', y', z', 'r-','LineWidth',2)
grid on
axis equal
view(2)
end