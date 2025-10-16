% points_density

clear
clc

grid_size = 25; % number of points in each dimention.
% a0 = inputs.a0; % lattice constant at rest
rho_0 = 1;
rho_max = 5;
sigma_r = 15;
r_0 = 1:1:grid_size;
size_x = 15;
size_y = size_x;
N_circles = 20;
% cT = 3;
% alpha = 1/137.04;
% cL = 1.618*cT;

% rho_p = rho_0*(1 + rho_max*exp(-0.5*(r_0./sigma_r).^2));
% rho_p = rho_0*(1 + rho_max./r_0.^2);
rho_p = rho_0*(1 + rho_max*exp(-r_0./sigma_r));
rho_n = 1./(2./rho_0 - 1./rho_p);

r_p = r_0./rho_p;
r_n = r_0./rho_n;
%
%%
figure(13)
clf
subplot(2,2,1)
plot(r_0,rho_0*ones(size(r_0))','.-')
hold on
plot(r_p,rho_p,'.-')
plot(r_n,rho_n,'.-')
% plot(cumsum(r_p),rho_p,'.-')
% plot(cumsum(r_n),rho_n,'.-')
hold off
legend

subplot(2,2,2)
plot(r_0,1.0*ones(size(r_0)),'.-')
hold on
plot(r_p,1.5*ones(size(r_p)),'.-')
plot(r_n,0.5*ones(size(r_n)),'.-')
hold off
%
%%
figure(14)
clf
subplot(2,2,1)
hold on
diff_r_p = diff(r_p);
X_r_p = [];
Y_r_p = [];

Delta_r_p = zeros(N_circles,1);

for m = 1:N_circles
    
    N_theta_p = round(2*pi*r_p(m)./diff_r_p(m));
    theta_p = linspace(0,2*pi*(1 - 0),N_theta_p);
    delta_r_p = m - N_theta_p/(2*pi); 
    Delta_r_p(m) = delta_r_p; 

    [x_r_p,y_r_p] = pol2cart(theta_p,r_p(m));
    X_r_p = [X_r_p; x_r_p'];
    Y_r_p = [Y_r_p; y_r_p'];
    
    
end

X_r_p = [0; X_r_p];
Y_r_p = [0; Y_r_p];

plot(X_r_p, Y_r_p,'.', 'MarkerSize',2)

hold off
axis([-size_x size_x -size_y size_y])
axis equal

%
%%

[v_p, c_p] = voronoin([X_r_p, Y_r_p]);
% figure(14)
subplot(2,2,3)
voronoi(X_r_p, Y_r_p)
axis([-size_x size_x -size_y size_y])
axis equal

%
%%
% figure(14)
subplot(2,2,2)
hold on
diff_r_n = diff(r_n);
X_r_n = [];
Y_r_n = [];

Delta_r_n = zeros(N_circles,1);

for m = 1:N_circles
    
    delta_theta = diff_r_n(m)./r_n(m);
    N_theta_n = round(2*pi*r_n(m)./diff_r_n(m));
    theta_n = linspace(0,2*pi*(1 - 0),N_theta_n);
    delta_r_n = m - N_theta_n/(2*pi); 
    Delta_r_n(m) = delta_r_n; 

    [x_r_n,y_r_n] = pol2cart(theta_n,r_n(m));
    X_r_n = [X_r_n; x_r_n'];
    Y_r_n = [Y_r_n; y_r_n'];
   
end

X_r_n = [0; X_r_n];
Y_r_n = [0; Y_r_n];

plot(X_r_n, Y_r_n,'.', 'MarkerSize',2)

hold off
axis([-size_x size_x -size_y size_y])
axis equal
%
%%
[v_n, c_n] = voronoin([X_r_n, Y_r_n]);
% figure(15)
subplot(2,2,4)
voronoi(X_r_n, Y_r_n)
% plot(v_n(:,1), v_n(:,2), '-');
% patch(v_n(:,1), v_n(:,2), 'b');
axis([-size_x size_x -size_y size_y])
axis equal
%
%%
figure(13)
subplot(2,2,3)
plot(r_p(1:N_circles), Delta_r_p, '.-');
hold on
plot(r_n(1:N_circles), Delta_r_n, '.-');
hold off
%
%% 2D:
% xx = 1:1:grid_size;
% yy = 1:1:grid_size;
% [X0, Y0] = meshgrid(xx,yy);
% mu_x = grid_size/2.0;
% mu_y = grid_size/2.0;
% [theta, R0] = cart2pol(X0 - mu_x, Y0 - mu_y);

% RHO2p = rho_0*(1 + rho_max*exp(-0.5*(R0./sigma_r).^2));
% RHO2p = rho_0*(1 + rho_max*exp(-0.5*(r_0./sigma_r).^2));
% RHO2n = 1./(2./rho_0 - 1./RHO2p);
% 
% R2p = R0./RHO2p;
% R2n = R0./RHO2n;

% K2p = (pi/12)*(gradient(R2p).^2)/R2p.^2;
% K2n = (pi/12)*(gradient(R2n).^2)/R2n.^2;

% K2p = (4*pi/45)*(gradient(R2p).^2)/R2p.^2;
% K2n = (4*pi/45)*(gradient(R2n).^2)/R2n.^2;
%


% [v, c] = voronoi(x, y);  % compute Voronoi diagram
% for i = 1:N
%     for j = (i+1):N
%         if v(i) == v(j)  % check if cells share an edge
%             % cells i and j are neighbors
%         end
%     end
% end







