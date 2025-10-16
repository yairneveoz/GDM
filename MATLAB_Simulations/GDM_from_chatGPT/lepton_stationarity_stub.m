function out = lepton_stationarity_stub(theta_geom, theta_medium)
% LEPTON_STATIONARITY_STUB
% Stub for Branch A: returns predicted alpha (r/R) and residuals for a
% stationary closed-loop lepton-at-rest configuration.
%
% INPUTS
%   theta_geom   = [R, r, n]      % geometry: major/minor radius, integer winding
%   theta_medium = struct with fields:
%       K    : effective stiffness / "bulk" elastic scale (>0)
%       chi  : spin-flow to contraction coupling (>0)
%       xi0  : dimensionless bias (placeholder for nonlinear micro-physics)
%
% OUTPUT (struct)
%   alpha_pred  : predicted dimensionless ratio r/R (here simply r/R)
%   EL_res      : Euler–Lagrange "stationarity" residual (toy)
%   Prad_res    : no-radiation proxy residual (toy)
%   mass_pred   : toy mass estimate from R (ħ=c=1 units → m ~ 1/R)
%   info        : diagnostic fields

    R  = theta_geom(1);
    r  = theta_geom(2);
    n  = theta_geom(3);

    K   = theta_medium.K;
    chi = theta_medium.chi;
    xi0 = theta_medium.xi0;

    % --- Dimensionless ratio actually realized by the geometry:
    alpha_pred = r / R;

    % --- Placeholder "preferred" ratio from micro-physics (REPLACE with real GDM):
    % Think of xi_star = f(K,chi,...) as the stationary r/R given your Lagrangian.
    xi_star = chi ./ (4*pi*K) + xi0;   % simple positive function (toy)

    % --- Stationarity residual (should -> 0 at solution):
    EL_res = alpha_pred - xi_star;

    % --- No-radiation residual (toy): small for closed, low-n loop
    Prad_res = 1e-4*(n-1)^2 + 1e-5*(alpha_pred - xi_star).^2;

    % --- Mass proxy (ħ=c=1 → m ~ 1/R); put back units if needed
    mass_pred = 1.0 / R;

    out.alpha_pred = alpha_pred;
    out.EL_res     = EL_res;
    out.Prad_res   = Prad_res;
    out.mass_pred  = mass_pred;
    out.info       = struct('xi_star',xi_star,'R',R,'r',r,'n',n, ...
                            'K',K,'chi',chi,'xi0',xi0);
end
