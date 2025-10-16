% Quarks
clc
clear
%{
up
down
top
bottom
strange
charm

anti_up
anti_down
anti_top
anti_bottom
anti_strange
anti_charm
%}
%% Plot sphere:


%
%% Plot paths:
delta_t = 0.02;
t = 0:delta_t:4;
R = 1;
theta = 2*pi/1; % 1, 
phi = theta/2; % 2, 3
Rt = [];

x = R*cos(theta*t).*sin(phi*t);
y = R*sin(theta*t).*sin(phi*t);
z = R*cos(phi*t);

xq = zeros(size(t));
yq = zeros(size(t));
zq = zeros(size(t));

xqv = x;
yqv = y;
zqv = z;

vx = [0, diff(xqv)/delta_t];
vy = [0, diff(yqv)/delta_t];
vz = [0, diff(zqv)/delta_t];

B = cross([vx' vy' vz'],[xqv' yqv' zqv']);
Bx = B(:,1)';
By = B(:,2)';
Bz = B(:,3)';

%
%% Plot figure:
figure(12)
plot3(x,y,z,'wo-')
hold on
quiver3(xq,yq,zq,xqv,yqv,zqv,0,'w')
quiver3(x,y,z,vx,vy,vz,1,'c')
quiver3(x,y,z,Bx,By,Bz,1,'m')
% sphere
hold off
axis equal
grid on
axis([-R R -R R -R R])
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'Color',0.5*[1 1 1])
%
%{
Q: Density of space
m: Energy of space density
s: Integral of mangetic field over one period
p: m*v

hbar: Relates to angular momentum space
G: 1/10^42
alpha: Ratio of r0 to R
c: translational, longitidual

Feynman - all possible paths
Lagrangian
%}
