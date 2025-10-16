% DEMO_ALPHA_G_COUPLING
% Quick check that the stubs are wired correctly.

clear; clc;

% --- geometry guess (electron-like scale just for a number)
R = 3.86e-13;                   % m (reduced Compton scale)
r = R / 140;                    % start near alpha, but not exact
n = 1;

% --- medium parameters (toy)
theta_medium = struct( ...
    'K',    2.0, ...            % arbitrary positive
    'chi',  0.011, ...          % tune to shift xi_star
    'xi0',  0.0, ...            % small bias
    'a_map', 6.6743e-11 ...     % choose so G_pred ~ 6.67e-11 when chi/K ~ 1
);

% --- Branch A: stationarity
A = lepton_stationarity_stub([R,r,n], theta_medium);
fprintf('alpha_pred = r/R = %.9f   (xi_star=%.9f)   EL_res=%.2e\n', ...
    A.alpha_pred, A.info.xi_star, A.EL_res);

% --- Branch B: G from frame-drag
G_pred = G_from_framedrag_stub(theta_medium);
fprintf('G_pred = %.6e (target ~ 6.674e-11)\n', G_pred);
