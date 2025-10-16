% Walter_Robertson
clear
clc

theta = 32;

e1=[cos(theta) -sin(theta) 0];
e2=[sin(theta) cos(theta) 0];
e3=[0 0 1];

% e1=[a1 b1 c1];
% e2=[a2 b2 c2];
% e3=[a3 b3 c3];

n1=[1 0 0];    % Ox
n2=[0 1 0];    % Oy
n3=[0 0 1];    % Oz

R = [e1; e2; e3];

%{
is the rotation matrix already, when we assume, 
that these are the normalized orthogonal vectors 
of the local coordinate system. To convert between 
the two reference systems all you need is R and R.' 
(as long as the translation is ignored).
%}

% A vector v=[x;y;z] in the global reference system 
x = n1;
y = n2;
z = n3;

v = [x; y; z];
v2 = R*v;

v0 = [1 0 0]';
v1 = R*v0;
% in the local system. Then you can convert it back 
% to the global system by:
% A2 = R.'*R*v;

figure(5)
quiver3(0, 0, 0, v0(1),v0(2),v0(3))
hold on
quiver3(0, 0, 0, v1(1),v1(2),v1(3))
hold off
axis equal
grid on



