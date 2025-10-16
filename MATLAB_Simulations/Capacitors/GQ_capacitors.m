% GQ_capacitors

clear
clc

max_r = 20;
surface_r = 6.4;
rr = linspace(surface_r,max_r,1000);

% G:
r_G = 0;
M = 1e5;
m = 1;
G = 0.001;
F_G = G*M*m./rr.^2;
F_G = 1*ones(size(rr));
% Q_m:
W_n = 0.5;
W_p = W_n;
Q_n = 0.05;
Q_p = Q_n;
r_Q_n = surface_r;
r_Q_p = surface_r + 7;
% F_Q_m = -1e-4*Q_m*Q_m./(rr - r_Q_m).^2;
% F_Q_p =  1e-6*Q_p*Q_p./(rr - r_Q_p).^2;
F_Q_n = -2*Q_n*1e1*exp(-abs(rr - r_Q_n)./W_n);
F_Q_p =  1*Q_p*1e1*exp(-abs(rr - r_Q_p)./W_p);

sum_F = F_G + F_Q_n + F_Q_p;
sum_F_n = F_G + F_Q_n;
% prod_F = F_G.*F_Q_m.*F_Q_p;
figure(1)
plot(rr,F_G)
hold on
plot(rr,F_Q_n)
plot(rr,F_Q_p)
plot(rr,sum_F_n)
hold off
% axis([0 inf 0 inf])
legend('F_G','F_{Qn}','F_{Qp}','sum F_n')
