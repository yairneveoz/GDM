function plot2D_dynamics_YZ(inputs)

grid_size = 81; %inputs.grid_size;

%% Colormaps:
% M = dlmread(filename,delimiter);
rho_cmap = parula64;
length_rho_cmap = length(rho_cmap);
rho_cmap_p = rho_cmap(length_rho_cmap/2+1:length_rho_cmap,:);
rho_cmap_n = rho_cmap(1:length_rho_cmap/2,:);

Nc = 64; %inputs.Nc;

reds1 = ones(1,Nc);
greens1 = linspace(1,0,Nc);
blues1 = linspace(1,0,Nc);
white_red_colormap = [reds1', greens1', blues1'];
q_cmap_p = white_red_colormap;

reds2 = linspace(0,1,Nc);
greens2 = linspace(0,1,Nc);
blues2 = ones(1,Nc);
blue_white_colormap = [reds2', greens2', blues2'];
q_cmap_n = blue_white_colormap;
%
%% Set parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D transverse wave
% 2D,1D longitudinal wave
% mixed

fsc = 1/137.036; % fine structure constant
cT = 3e10; %cm/sec
cL = cT*1.6068;
r_p = 0.8774e-13; %cm
r_e = (cL/cT)*r_p;
R_p = r_p/fsc;



v = 0.8*cT;
gamma = 1/sqrt(1 - v^2/cT^2);

xx = 1:1:(1.2*gamma)*grid_size;
yy = 1:1:grid_size;
[X0, Y0] = meshgrid(xx,yy);
grid_size_x = length(xx);
grid_size_y = length(yy);

rotation_R0 = 0.5;
rotation_R = rotation_R0/gamma;
rho_0 = 1;
rho_max = 4;
sigma_r_to_grid_size = 1/20;
sigma_r = grid_size*sigma_r_to_grid_size/gamma;

theta = 0:0.1:2*pi;
r0p = 0.5*rotation_R0*grid_size*ones(size(theta));
r0n = 0.5*rotation_R0*grid_size*ones(size(theta));
% [x0p,y0p] = pol2cart(theta,r0p);
% [x0n,y0n] = pol2cart(theta,r0n);
%
%% Plot
figure(24) 
clf
    
t = 0:0.1:2*pi;
for t_ind = 1:length(t)
%     (grid_size_x/10)
%     mu_xt = (grid_size/2)*(1 - rotation_R*cos(t(t_ind)));
    mu_xt = grid_size_x*(1/20 + 0.1*(v/cT)*t(t_ind)); % 
    mu_yt = grid_size_y*(1/2 + 0.5*rotation_R*sin(t(t_ind)));
%     mu_xt2 = mu_x + 8*t(t_ind);
    [theta_space, R0_space] = cart2pol(X0 - mu_xt, Y0 - mu_yt);

    RHO2_p = rho_0*(1 + rho_max*exp(...
        -0.5*(((X0 - mu_xt).^2 + (Y0 - mu_yt).^2)./sigma_r.^2)));
    RHO2_n = 1./(2./rho_0 - 1./RHO2_p);

%     R2_p = R0_space./RHO2_p;
%     R2_n = R0_space./RHO2_n;
    %
    %% Get 2D grids:
%     [Xp, Yp] = pol2cart(theta_space, R2_p);
%     [Xn, Yn] = pol2cart(theta_space, R2_n);
    
    %
    %% rho(r0) and q(r0):
    % q = 1/(4*pi)*(rho - rho_0)./rho;
%     q2p = 1/(4*pi)*(RHO2_p - rho_0)./RHO2_p;
%     q2n = 1/(4*pi)*(RHO2_n - rho_0)./RHO2_n;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Xp,Yp,q2p,Xn,Yn,q2n] = ...
    dynamicPathYZ(X0,Y0,sigma_r,rotation_R,t(t_ind),v/cT);
    %
    %% (2,2,1):
    ax1 = subplot(2,1,1);
    surf(Xp + mu_xt - grid_size/2, Yp + mu_yt - grid_size/2, q2p,...
        'EdgeAlpha',0.1);
    colormap(ax1, q_cmap_p);
%     hold on
%     plot([0 grid_size_x], 0.5*rotation_R0*grid_size*[1 1],'-',...
%         'Color', 0*[1 1 1], 'LineWidth', 1.0)
%     plot([0 grid_size_x], -0.5*rotation_R0*grid_size*[1 1],'-',...
%         'Color', 0*[1 1 1], 'LineWidth', 1.0)   
%     hold off
    view(2)
    
%     axis tight
%     axis square
    xlabel('X')
    ylabel('Y')
    title(['Deformed 2D q_+ space, v = ', num2str(v/cT),'c'])
%     axis([0, 4*grid_size, (grid_size/2)*[-1 1]])
    axis tight
    axis equal   
%
    %% (2,1,2) q(r):
    ax2 = subplot(2,1,2);
    surf(Xn+mu_xt-grid_size/2, Yn+mu_yt-grid_size/2, q2n,...
        'EdgeAlpha',0.1)
    colormap(ax2, q_cmap_n);
%     hold on
%     plot([min(xx) max(xx)], 0.5*rotation_R0*grid_size*[1 1],'-',...
%         'Color', 0*[1 1 1], 'LineWidth', 1.0)
%     plot([0 grid_size_x], -0.5*rotation_R0*grid_size*[1 1],'-',...
%         'Color', 0*[1 1 1], 'LineWidth', 1.0)   
%     hold off
    view(2)
%     axis tight
%     axis square
    xlabel('X')
    ylabel('Y')
    title(['Deformed 2D q_- space, v = ', num2str(v/cT),'c'])
%     axis((grid_size/2)*[-1 1 -1 1])
%     axis((grid_size/2)*[-1 1 -1 1])
    axis tight
    axis equal
    pause(0.05)
    drawnow
%
end
end