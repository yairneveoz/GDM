% Continous_space_density

clear
clc

grid_size = 300; %25; % number of points in each dimention.
rho_0 = 1;
rho_max = 5;
sigma_r = 15;

%
%% 2D:
% Set initial grid:
xx = 1:1:grid_size;
yy = 1:1:grid_size;
[X0, Y0] = meshgrid(xx,yy);
% [theta0,R0] = [];
RHO_0 = ones(size(X0));

% Set coordinates of deformation center:
cx_1 = grid_size*0.25;
cy_1 = grid_size*0.5;

cx_2 = grid_size*0.75;
cy_2 = grid_size*0.5;

X = X0;
Y = Y0;
RHO = RHO_0;
[X1, Y1, RHO_1] = ...
    deformSpacePositive(X0,Y0,RHO_0,cx_1,cy_1,rho_max,sigma_r);
DX1 = X1 - X0 + cx_1;
DY1 = Y1 - Y0 + cy_1;
% [X2,Y2,RHO_2] = ...
%    deformSpaceNegative(X0,Y0,RHO_0,cx_2,cy_2,rho_max,sigma_r);
[X2, Y2, RHO_2] = ...
    deformSpaceNegative(X0,Y0,RHO_0,cx_2,cy_2,rho_max,sigma_r);

DX2 = X2 - X0 + cx_2;
DY2 = Y2 - Y0 + cy_2;

DX12 = DX1 + DX2;
DY12 = DY1 + DY2;

DR1 = sqrt(DX1.^2 + DY1.^2);
DR2 = sqrt(DX2.^2 + DY2.^2);
DR12 = sqrt(DX12.^2 + DY12.^2);

RHO_p = 1./(1 + DR1);
RHO_n = 1./(1 + DR2);

RHO_pn = 1./(1 + DR12);

X12 = X0 + DX12;
Y12 = Y0 + DY12;
%
%%

figure(16)
subplot(3,3,1)
surf(RHO_1, 'EdgeColor','none')
hold on
contour3(RHO_1, 'k')
hold off

%%%
subplot(3,3,2)
surf(RHO_2, 'EdgeColor','none')


%%%
% figure(17)
% plot(X1(:),Y1(:),'.','MarkerSize',1)
% 
% figure(18)
% plot(X2(:),Y2(:),'.','MarkerSize',1)

% figure(19)
% X3 = X0 + DX12;
% Y3 = Y0 + DY12;
% plot(X3(:),Y3(:),'.','MarkerSize',1)
% %%%
subplot(3,3,3)
surf(RHO_pn, 'EdgeColor','none')
% hold on
% contour3(RHO2p_1.*RHO2n_2, 20, 'k')
% hold off

%% %
figure(16)
subplot(3,3,4)
surf(X1, Y1, RHO_1, 'EdgeColor','none')
% % plot(X2p_1, Y2p_1, '.')
% hold on
% contour3(X2p_1, Y2p_1, RHO2p_1, 'k')
% hold off

%% %
subplot(3,3,5)
surf(X2, Y2, RHO_2, 'EdgeColor','none')
% % plot(X2p_1, Y2p_1, '.')
% hold on
% contour3(X2n_2, Y2n_2, RHO2n_2, 'k')
% hold off

%% %
% subplot(2,3,6)
% surf(X2n_2, Y2n_2, RHO2n_2, 'EdgeColor','none')
% % plot(X2p_1, Y2p_1, '.')
% hold on
% contour3(X2n_2, Y2n_2, RHO2n_2, 'k')
% hold off







