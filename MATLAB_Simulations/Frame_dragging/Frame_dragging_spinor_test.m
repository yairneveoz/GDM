% Frame_dragging_spinor_test

clc
clear

r_half = 3;
r_max = 20;
r_N = 400;
rr2 = linspace(r_max,r_half,r_N);

phi = linspace(0,2*pi,40);
x_half = 1.5*r_half*cos(phi);
y_half = 1.5*r_half*sin(phi);

sr = 1;
rho_r2 = sr*exp(-abs(rr2));

theta0 = 40*2*pi;
theta2 = theta0*rho_r2;

x02 = cos(theta2);
y02 = sin(theta2);

% sum_theta2 = cumsum(delta_theta2,delta_theta2);

[x2,y2] = pol2cart(theta2,rr2);
% y2 = sin(sum_theta2);

figure(1)
clf
subplot(2,2,1)
plot(rr2, rho_r2, '.-')
hold on
plot(rr2(1), rho_r2(1), 'o')
hold off
axis([0 r_max 0 inf])
xlabel('r')
ylabel('\rho(r)')

subplot(2,2,2)
plot(x02, y02)
hold on
plot(x02(1), y02(1), 'o')
hold off
xlabel('\rho')
ylabel('\theta(r)')
% axis([0 r_max 0 inf])

subplot(2,2,3)
plot(x2, y2)
hold on
plot(x2(1), y2(1), 'o')
plot(x_half, y_half, 'r-')
hold off
xlabel('x')
ylabel('y')
axis equal
% axis([0 r_max 0 inf])


