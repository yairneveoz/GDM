function outputs = get1D(inputs)

grid_size = inputs.grid_size; % number of points in each dimention.
% a0 = inputs.a0; % lattice constant at rest
rho_0 = inputs.rho_0;
rho_max = inputs.rho_max;

cT = inputs.cT;
alpha = inputs.alpha;
cL = inputs.cL;
 
% Gaussians:
r0 = 0:1:floor(grid_size/2);
sigma_r = inputs.sigma_r;

% rho1_p = rho_0*(1 + rho_max*exp(-0.5*(r0./sigma_r).^2));
rho1_p = rho_0*(1 + rho_max*exp(-abs(r0./sigma_r)));
rho1_n = 1./(2./rho_0 - 1./rho1_p);

% Get new r distributions:
delta_1p = 1./rho1_p;
delta_1n = 1./rho1_n;


rp0 = cumsum(delta_1p) - sum(delta_1p(1:floor(grid_size/2)));
rn0 = cumsum(delta_1n) - sum(delta_1n(1:floor(grid_size/2)));
rp = rp0 - min(rp0);
rn = rn0 - min(rn0);

q1p = 1/(4*pi)*(rho1_p - rho_0)./rho1_p;
q1n = 1/(4*pi)*(rho1_n - rho_0)./rho1_n;

%% q(r), Rs
sum_q1p = sum(q1p);
half_sum_q10p = sum_q1p/2;

sum_q1n = sum(q1n);
half_sum_q10n = sum_q1n/2;

%% 2D Q(r):
Q1p = cumsum(q1p);
Q1n = cumsum(q1n);

E1p = (delta_1p - ones(size(delta_1p))).^2;
E1n = (delta_1n - ones(size(delta_1n))).^2;

m1p = E1p/cT^2;
m1n = E1n/cT^2;

%% Outputs:
outputs.r0 = r0;
outputs.rho1_p = rho1_p;
outputs.rho1_n = rho1_n;
outputs.rp = rp;
outputs.rn = rn;
outputs.q1p = q1p;
outputs.q1n = q1n;
outputs.Q1p = Q1p;
outputs.Q1n = Q1n;
outputs.E1p = E1p;
outputs.E1n = E1n;
outputs.m1p = m1p;
outputs.m1n = m1n;
outputs.half_sum_q10p = half_sum_q10p;
outputs.half_sum_q10n = half_sum_q10n;
%
end


