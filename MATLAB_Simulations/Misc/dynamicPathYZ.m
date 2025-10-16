function [Xp,Yp,q2p,Xn,Yn,q2n] = ...
    dynamicPathYZ(X0,Y0,sigma_r,rotation_R,t,v)

rho_0 = 1;
rho_max = 4;


grid_size_x = size(X0,2);
grid_size_y = size(X0,1);


% mu_xt = (grid_size_x/2)*(1 - rotation_R*cos(omega_t));
% mu_xt = (grid_size_x/2)*v*t;
mu_xt = grid_size_x*(1/20 + 0.1*v*t); % /gamma
mu_yt = (grid_size_y/2)*(1 + rotation_R*sin(t));

[theta_space, R0_space] = cart2pol(X0 - mu_xt, Y0 - mu_yt);

RHO2_p = rho_0*(1 + rho_max*exp(...
    -0.5*(((X0 - mu_xt).^2 + (Y0 - mu_yt).^2)./sigma_r.^2)));

RHO2_n = 1./(2./rho_0 - 1./RHO2_p);

R2_p = R0_space./RHO2_p;
R2_n = R0_space./RHO2_n;

%% Get deformed 2D grids:
[Xp, Yp] = pol2cart(theta_space, R2_p);
[Xn, Yn] = pol2cart(theta_space, R2_n);
    
%% q(r0):
% q = 1/(4*pi)*(rho - rho_0)./rho;
q2p = 1/(4*pi)*(RHO2_p - rho_0)./RHO2_p;
q2n = 1/(4*pi)*(RHO2_n - rho_0)./RHO2_n;

end