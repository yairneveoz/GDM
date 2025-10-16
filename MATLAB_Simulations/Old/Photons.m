% Photons
clc
clear

%{

%}
%% Plot sphere:


%
%% Plot paths:
delta_t = 0.02;
t = 0:delta_t:2;
R = 1;
theta = 0;
phi = 2*pi;
c = 1;

x = R*cos(phi*t);
y = R*sin(phi*t);
z = t;

% vx0 = R*-sin(phi*delta_t);
% vy0 = R*co(phi*delta_t);
% vz0 = v;

xq = zeros(size(t));
yq = zeros(size(t));
zq = t;

xqv = x;
yqv = y;
zqv = zeros(size(t));

vx = [0, diff(x)/delta_t];
vy = [0, diff(y)/delta_t];
vz = [0, diff(z)/delta_t];

B = cross([vx' vy' vz'],[xqv' yqv' zqv']);
Bx = B(:,1)';
By = B(:,2)';
Bz = B(:,3)';
%
%% Plot figure:
figure(13)
clf
plot3(x,y,z,'wo-')
hold on
quiver3(xq,yq,zq,xqv,yqv,zqv,0,'w')
quiver3(x,y,z,vx,vy,vz,1,'c')
quiver3(x,y,z,Bx,By,Bz,1,'m')
hold off
grid on
set(gca,'Color',0.5*[1 1 1])
% axis([-R R -R R 0 max(t)*c])
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
legend({'e^-','R_{e^-}','c','B'})

