function outputs = get2D(inputs)

grid_size = inputs.grid_size; % number of points in each dimention.
% a0 = inputs.a0; % lattice constant at rest
rho_0 = inputs.rho_0;
rho_max = inputs.rho_max;
sigma_r = inputs.sigma_r;

cT = inputs.cT;
alpha = inputs.alpha;
cL = inputs.cL;

%
%% 2D:
xx = 1:1:grid_size;
yy = 1:1:grid_size;
[X0, Y0] = meshgrid(xx,yy);
mu_x = grid_size/2.0;
mu_y = grid_size/2.0;
[theta, R0] = cart2pol(X0 - mu_x, Y0 - mu_y);

RHO2p = rho_0*(1 + rho_max*exp(-0.5*(R0./sigma_r).^2));
RHO2n = 1./(2./rho_0 - 1./RHO2p);

R2p = R0./RHO2p;
R2n = R0./RHO2n;

K2p = (4*pi/45)*(gradient(R2p).^2)/R2p.^2;
K2n = (4*pi/45)*(gradient(R2n).^2)/R2n.^2;
%
%% Get 2D grids:
[Xp, Yp] = pol2cart(theta, R2p);
[Xn, Yn] = pol2cart(theta, R2n);

%
%% rho(r0) and q(r0):
% q = 1/(4*pi)*(rho - rho_0)./rho;
q2p = 1/(4*pi)*(RHO2p - rho_0)./RHO2p;
q2n = 1/(4*pi)*(RHO2n - rho_0)./RHO2n;
sum_q = q2p + q2n;
% test, sum_q should be 0:
% r_inds = [ceil(grid_size/2),ceil(grid_size/2):grid_size];
% disp(r_inds)
r1p = R2p(ceil(grid_size/2),ceil(grid_size/2):grid_size);
r1n = R2n(ceil(grid_size/2),ceil(grid_size/2):grid_size);

q1p = q2p(ceil(grid_size/2),ceil(grid_size/2):grid_size);
q1n = q2p(ceil(grid_size/2),ceil(grid_size/2):grid_size);

diff0_r1p = diff(r1p);
diff0_r1n = diff(r1n);

diff_r1p = [diff0_r1p,diff0_r1p(end)];
diff_r1n = [diff0_r1n,diff0_r1n(end)];

Q1p = 2*pi*cumsum(r1p.*q1p.*diff_r1p);
Q1n = 2*pi*cumsum(r1n.*q1n.*diff_r1n);

E1p = Q1p./(2*r1p.^2);
E1n = Q1n./(2*r1n.^2);

U2p = q2p.^2./R2p;
U2n = q2n.^2./R2n;
% Q2p = 2*pi*sum(q2p.*R2_p*dr);
%
%% q(r): Plot Schwartzschild rudius:
r0p = 1.2;
r0n = 7.2;
phi = 0:0.1:2*pi;
x0p = r0p*cos(phi);
y0p = r0p*sin(phi);

x0n = r0n*cos(phi);
y0n = r0n*sin(phi);

%% Get rings:

%
%% Energy calculations:
diff_X00 = diff(X0,1,2);
diff_Y00 = diff(Y0,1,1);

diff_X0 = [diff_X00(:,1), diff_X00];
diff_Y0 = [diff_Y00(1,:); diff_Y00];

% delta_tX0 = cT./diff_X0;
% delta_tY0 = cT./diff_Y0;
delta_tX0 = 1./diff_X0;
delta_tY0 = 1./diff_Y0;

diff_Xp0 = diff(Xp,1,2);
diff_Yp0 = diff(Yp,1,1);
diff_Xn0 = diff(Xn,1,2);
diff_Yn0 = diff(Yn,1,1);

diff_Xp = [diff_Xp0(:,1), diff_Xp0];
diff_Yp = [diff_Yp0(1,:); diff_Yp0];
diff_Xn = [diff_Xn0(:,1), diff_Xn0];
diff_Yn = [diff_Yn0(1,:); diff_Yn0];



% delta_tXp = cT./diff_Xp;
% delta_tYp = cT./diff_Yp;
% delta_tXn = cT./diff_Xn;
% delta_tYn = cT./diff_Yn;
delta_tXp = 1./diff_Xp;
delta_tYp = 1./diff_Yp;
delta_tXn = 1./diff_Xn;
delta_tYn = 1./diff_Yn;

k = 1; % Spring constant
% E_Xp = 0.5*k*(diff_Xp - diff_X0).^2;
% E_Yp = 0.5*k*(diff_Yp - diff_Y0).^2;
% E_Xn = 0.5*k*(diff_Xn - diff_X0).^2;
% E_Yn = 0.5*k*(diff_Yn - diff_Y0).^2;
diff_Xq2p0 = 4*pi*diff(q2p,1,2);
diff_Yq2p0 = 4*pi*diff(q2p,1,1);
diff_Xq2n0 = 4*pi*diff(q2n,1,2);
diff_Yq2n0 = 4*pi*diff(q2n,1,1);

diff_Xq2p = [diff_Xq2p0(:,1), diff_Xq2p0];
diff_Yq2p = [diff_Yq2p0(1,:); diff_Yq2p0];
diff_Xq2n = [diff_Xq2n0(:,1), diff_Xq2n0];
diff_Yq2n = [diff_Yq2n0(1,:); diff_Yq2n0];

E_Xp = 0.5*k*(diff_Xq2p).^2;
E_Yp = 0.5*k*(diff_Yq2p).^2;
E_Xn = 0.5*k*(diff_Xq2n).^2;
E_Yn = 0.5*k*(diff_Yq2n).^2;

% E1p = sqrt(E_Xp.^2 + E_Yp.^2);
% E1n = sqrt(E_Xn.^2 + E_Yn.^2);

%
%% Outputs:
outputs.R0 = R0;
outputs.X0 = X0;
outputs.Y0 = Y0;
outputs.RHO2_p = RHO2p;
outputs.RHO2_n = RHO2n;
outputs.R2_p = R2p;
outputs.R2_n = R2n;
outputs.Xp = Xp;
outputs.Yp = Yp;
outputs.Xn = Xn;
outputs.Yn = Yn;
outputs.q2p = q2p;
outputs.q2n = q2n;
outputs.Q1p = Q1p;
outputs.Q1n = Q1n;
outputs.E1p = E1p;
outputs.E1n = E1n;
outputs.U2p = U2p;
outputs.U2n = U2n;
% outputs.m2p = m2p;
% outputs.m2n = m2n;
% outputs.half_sum_q20p = half_sum_q20p;
% outputs.half_sum_q20n = half_sum_q20n;

outputs.diff_X0 = diff_X0;
outputs.diff_Y0 = diff_Y0;
outputs.diff_Xp = diff_Xp;
outputs.diff_Yp = diff_Yp;
outputs.diff_Xn = diff_Xn;
outputs.diff_Yn = diff_Yn;
outputs.delta_tX0 = delta_tX0;
outputs.delta_tY0 = delta_tY0;
outputs.delta_tXp = delta_tXp;
outputs.delta_tYp = delta_tYp;
outputs.delta_tXn = delta_tXn;
outputs.delta_tYn = delta_tYn;
outputs.K2p = K2p;
outputs.K2n = K2n;
%
end


