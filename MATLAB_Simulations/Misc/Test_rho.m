% Test_rho
clear
clc

rr = -10:0.05:10;
rho_0 = 4*ones(size(rr));
rho_p = 3./(rr.^2); %rho_0.*exp(-(rr.^2)); % 
rho_n = rho_0./(2 - rho_0./rho_p);

q_p = (1 - rho_0./rho_p)/(4*pi);
q_n = (1 - rho_0./rho_n)/(4*pi);


figure(3)
subplot(1,2,1)
plot(rr,rho_p)
hold on
plot(rr,rho_n)
hold off

subplot(1,2,2)
plot(rr,q_p)
hold on
plot(rr,q_n)
hold off

