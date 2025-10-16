% Optical_path

clear
clc

epsilon0 = 1;
epsilon_max = 2;

theta_in = 10; % deg
y_top = 1; % mm
yy = 0:0.2:5;
leyers_y = y_top*exp(-yy);
epsilons_y = linspace(epsilon0,epsilon_max,length(yy));



figure(1)
plot(yy,leyers_y,'.-')
hold on
plot([yy(1) yy(end)],leyers_y(5).*[1 1],'-')
hold off




