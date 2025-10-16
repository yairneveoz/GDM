function outputs = get2Ddynamics()





%% 2D:
xx = 1:1:grid_size;
yy = 1:1:grid_size;
[X0, Y0] = meshgrid(xx,yy);
mu_x = grid_size/2.0;
mu_y = grid_size/2.0;
[theta, R0] = cart2pol(X0 - mu_x, Y0 - mu_y);

RHO2_p = rho_0*(1 + rho_max*exp(-0.5*(R0./sigma_r).^2));
RHO2_n = 1./(2./rho_0 - 1./RHO2_p);

R2_p = R0./RHO2_p;
R2_n = R0./RHO2_n;

%
%% Get 2D grids:
[Xp, Yp] = pol2cart(theta, R2_p);
[Xn, Yn] = pol2cart(theta, R2_n);
%
%% rho(r0) and q(r0):
% q = 1/(4*pi)*(rho - rho_0)./rho;
q2p = 1/(4*pi)*(RHO2_p - rho_0)./RHO2_p;
q2n = 1/(4*pi)*(RHO2_n - rho_0)./RHO2_n;
end