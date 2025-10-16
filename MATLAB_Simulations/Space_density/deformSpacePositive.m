function [X,Y,RHO] = ...
    deformSpacePositive(X0,Y0,RHO_0,cx,cy,rho_max,sigma_r)


[theta, R0] = cart2pol(X0 - cx, Y0 - cy);

RHO = RHO_0.*(1 + rho_max*exp(-0.5*(R0./sigma_r).^2));

R = R0./RHO;

[X, Y] = pol2cart(theta, R);

end