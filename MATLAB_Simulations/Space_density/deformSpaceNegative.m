function [X_n,Y_n,RHO_n] = ...
    deformSpaceNegative(X0,Y0,RHO_0,cx,cy,rho_max,sigma_r)


[theta, R0] = cart2pol(X0 - cx, Y0 - cy);
% [theta, R0] = cart2pol(X0, Y0);

RHO_p = RHO_0.*(1 + rho_max*exp(-0.5*(R0./sigma_r).^2));
RHO_n = 1./(2./RHO_0 - 1./RHO_p);

R_n = R0./RHO_n;

[X_n, Y_n] = pol2cart(theta, R_n);

end