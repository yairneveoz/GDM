% plot_EM
% Electromagnetism as the Geometrodynamics of Space
%{
The Electric Field of a Moving Charge
Let a point charge Q reside at the origin of a frame k.
The electric field E has the magnitude HQ/r2
, and for a positive charge it is directed radially 
outward.
In the xz plane its components at any point (x, z) are:
%}
c = 3e10; % cm/sec
v = 0.6*c;

Ex = H*Q/r.^2*cos(theta);
Ex = H*Q*x/(x.^2 + z.^2).^(3/2);

Ez = H*Q/r.^2*sin(theta);
Ez = H*Q*z/(x.^2 + z.^2).^(3/2);

gamma = 1./sqrt(1 - v^2/c^2);
