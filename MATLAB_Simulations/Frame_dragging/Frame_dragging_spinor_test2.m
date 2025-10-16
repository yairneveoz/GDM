% Frame_dragging_spinor_test2

clear
clc

N_grid_points = 100;
d0 = 0.5;
delta_x0 = d0*ones(N_grid_points,1);
sx0 = cumsum(delta_x0) - d0/2;

x0 = sx0 - mean(sx0);

lambda = 0.1;
rho_0 = 1;
delta_rho_xp = d0*rho_0*exp(-abs(lambda*x0));
delta_xp = delta_x0 - delta_rho_xp;

%%%
delta_rho_xm = rho_0.*(2 - rho_0./delta_rho_xp);
delta_xm = delta_x0.*(2 - delta_x0./delta_xp);
%%%
sxp = cumsum(delta_xp);
sxm = cumsum(delta_xm);

center_index1 = N_grid_points/2;
center_index2 = N_grid_points/2 + 1;
xp = sxp - mean([sxp(center_index1), sxp(center_index2)]);
xm = sxm - mean([sxm(center_index1), sxm(center_index2)]);

figure(11)
subplot(2,2,1)
plot(x0, delta_rho_xp, '.')

% subplot(2,2,2)
% plot(x0, delta_rho_xm, '.')

subplot(2,2,3)
% plot(sx, ones(size(sx)), '.')
plot(xp, ones(size(xp)), '.')
hold on
plot(x0, 0.5*ones(size(xm)), '.')
plot(xm, 0*ones(size(xm)), '.')
hold off
% xlim([min(x0), max(x0)])



