% Spirals2
clear
clc

R0 = 137; % 137; % 10
r0 = 1.0; % 1.2
Omega0 = 100;

n = 17.34; % 20; r cycles in R cycles
loops = 1;
t = 0:0.005:loops*2*pi;
c = 1;
v = 0.1*c;
gamma = 1/sqrt(1 - v^2/c^2);

R = R0/gamma;
r = r0/gamma;

X = R*cos(t);
Y = R*sin(t);
Z = R*v*t*gamma;

x = X + r.*cos(n*t).*cos(t);
y = Y + r.*cos(n*t).*sin(t);
z = Z + r.*sin(n*t);

% R^2 + 
% x = (R + r.*cos(n*t)).*cos(t);
% y = (R + r.*cos(n*t)).*sin(t);
% z = (r.*sin(n*t) + R*v*t*gamma);


figure(14)
plot3(x, y, z, 'b-','LineWidth',1)
hold on
plot3(X, Y, Z, '-','LineWidth',3, 'Color',[1.0, 0.0, 0.0, 0.3])
hold off

xlabel('X')
ylabel('Y')
zlabel('Z')

% axis equal
axis(1.2*[-R0 R0 -R0 R0 -1 max(Z)])
grid on
