function plot2D(inputs,outputs)

grid_size = inputs.grid_size;

X0 = outputs.X0;
Y0 = outputs.Y0;
RHO2_p = outputs.RHO2_p;
RHO2_n = outputs.RHO2_n;
Xp = outputs.Xp;
Yp = outputs.Yp;
R2_p = outputs.R2_p;
R2_n = outputs.R2_n;
Xn = outputs.Xn;
Yn = outputs.Yn;
q2p = outputs.q2p;
q2n = outputs.q2n;
% Q1p = outputs.Q1p;
% Q1n = outputs.Q1n;
U2p = outputs.U2p;
U2n = outputs.U2n;
E1p = outputs.E1p;
E1n = outputs.E1n;
% outputs.m2p = m2p;
% outputs.m2n = m2n;
% outputs.half_sum_q20p = half_sum_q20p;
% outputs.half_sum_q20n = half_sum_q20n;
% delta_tX0 = outputs.delta_tX0;
% delta_tY0 = outputs.delta_tY0;
% delta_tXp = outputs.delta_tXp;
% delta_tYp = outputs.delta_tYp;
% delta_tXn = outputs.delta_tXn;
% delta_tYn = outputs.delta_tYn;
%% 2D:
% M = dlmread(filename,delimiter);
rho_cmap = parula64;
length_rho_cmap = length(rho_cmap);
rho_cmap_p = rho_cmap(length_rho_cmap/2+1:length_rho_cmap,:);
rho_cmap_n = rho_cmap(1:length_rho_cmap/2,:);

Nc = inputs.Nc;

reds1 = ones(1,Nc);
greens1 = linspace(1,0,Nc);
blues1 = linspace(1,0,Nc);
white_red_colormap = [reds1', greens1', blues1'];
q_cmap_p = white_red_colormap;
rho_cmap_p = white_red_colormap;

reds2 = linspace(0,1,Nc);
greens2 = linspace(0,1,Nc);
greens3 = ones(1,Nc);
blues2 = ones(1,Nc);
blue_white_colormap = [reds2', greens2', blues2'];
cyan_white_colormap = [reds2', greens3', blues2'];
q_cmap_n = blue_white_colormap;
rho_cmap_n = cyan_white_colormap; %blue_white_colormap;

%%
figure(21)
clf
%% (2,4,1) rho(r0) and q(r0):
ax1 = subplot(2,4,1);
surf(X0, Y0, RHO2_p, 'EdgeAlpha',0.1)
colormap(ax1, rho_cmap_p)
xlabel('X_0')
ylabel('Y_0')
title('\rho_+(X_0,Y_0)')
%

%% (2,4,2)
ax2 = subplot(2,4,2);
surf(X0, Y0, q2p, 'EdgeAlpha',0.1)
colormap(ax2, q_cmap_p)
xlabel('X_0')
ylabel('Y_0')
title('q_+(X_0,Y_0)')
%
%% (2,4,3):

theta = 0:0.1:2*pi;
r0p = 1.2*ones(size(theta));
r0n = 4.0*ones(size(theta));
[x0p,y0p] = pol2cart(theta,r0p);
[x0n,y0n] = pol2cart(theta,r0n);

ax3 = subplot(2,4,3);
surf(Xp, Yp, q2p, 'EdgeAlpha',0.1);
colormap(ax3, q_cmap_p);
hold on
plot3(x0p, y0p, 4*ones(size(x0p)),'-',...
    'Color', 0*[1 1 1], 'LineWidth', 1.0)
hold off
view(2)
axis tight
axis square
xlabel('X')
ylabel('Y')
title('q_+ deformed 2D space')
%
%% (2,4,4)
ax4 = subplot(2,4,4);
surf(Xp, Yp, U2p, 'EdgeAlpha',0.1) % E2p % del2(q2p)
view(2)
colormap(ax4, rho_cmap_p)
% colormap(ax4, 'parula')
xlabel('X')
ylabel('Y')
title('U_+(X,Y)')
axis square
axis tight
% zlim([-4 1])
%
%% (2,4,5)
ax5 = subplot(2,4,5);
surf(X0, Y0, RHO2_n, 'EdgeAlpha',0.1)
colormap(ax5, rho_cmap_n)
xlabel('X_0')
ylabel('Y_0')
title('\rho_-(X_0,Y_0)')
zlim([-4 1])
%
%% (2,4,6)
ax6= subplot(2,4,6);
surf(X0, Y0, q2n, 'EdgeAlpha',0.1)
colormap(ax6, q_cmap_n)
xlabel('X_0')
ylabel('Y_0')
title('q_-(X_0,Y_0)')
%
%% q(r):
% ax3 = subplot(2,2,3);
% plot(Xn(20,:), q2n(20,:), 'b.-');
% hold on
% plot(7.2*[1 1],[min(q1n) 0],'b--')
% plot(-7.2*[1 1],[min(q1n) 0],'b--')
% hold off
%
%% (2,4,7) 
ax7 = subplot(2,4,7);
surf(Xn, Yn, q2n, 'EdgeAlpha',0.1)
colormap(ax7, q_cmap_n);
hold on
plot3(x0n, y0n, 1*ones(size(x0n)),'-',...
    'Color', 1*[1 1 1], 'LineWidth', 1.0)
hold off
view(2)
axis tight
axis square
xlabel('X')
ylabel('Y')
title('q_- deformed 2D space')
%
%% (2,4,8)
ax8 = subplot(2,4,8);
% surf(delta_tXn, delta_tYn, ones(size(delta_tXn)))
surf(Xn, Yn, U2n, 'EdgeAlpha',0.1) % E2n % del2(q2n)
view(2)
colormap(ax8, flipud(rho_cmap_n))
% colormap(ax8, 'parula')
xlabel('X')
ylabel('Y')
title('U_-(X,Y)')
axis square
axis tight
% zlim([-4 1])
%

figure(23)
subplot(2,1,1)
plot(X0(ceil(grid_size/2), ceil(grid_size/2):grid_size), E1p,...
    'r.-','LineWidth',1)
hold on
plot(X0(ceil(grid_size/2), ceil(grid_size/2):grid_size), E1n,...
    'b.-','LineWidth',1)
hold off

subplot(2,1,2)
plot(Xp(ceil(grid_size/2), ceil(grid_size/2):grid_size), E1p,...
    'r.-','LineWidth',1)
hold on
plot(Xn(ceil(grid_size/2), ceil(grid_size/2):grid_size), E1n,...
    'b.-','LineWidth',1)
hold off

%% (2,1,1)
%{
K2p = outputs.K2p;
K2n = outputs.K2n;
figure(24)
ax1 = subplot(2,1,1);
% surf(delta_tXn, delta_tYn, ones(size(delta_tXn)))
surf(Xp, Yp, K2p, 'EdgeAlpha',0.1) % E2n % del2(q2n)
view(2)
colormap(ax1, flipud(rho_cmap_p))
xlabel('X')
ylabel('Y')
title('K_+(X,Y)')
axis square
axis tight
%
%% (2,1,2)
ax2 = subplot(2,1,2);
% surf(delta_tXn, delta_tYn, ones(size(delta_tXn)))
surf(Xn, Yn, K2n, 'EdgeAlpha',0.1) % E2n % del2(q2n)
view(2)
colormap(ax2, flipud(rho_cmap_n))
xlabel('X')
ylabel('Y')
title('K_+(X,Y)')
axis square
axis tight
%
%}
end