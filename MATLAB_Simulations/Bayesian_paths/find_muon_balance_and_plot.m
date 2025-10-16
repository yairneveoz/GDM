%% find_muon_balance_and_plot.m
% Searches for a smooth, closed-loop "muon-at-rest" path and plots it.
% Identical structure to the electron script, but with muon constants.
% - Path model: toroidal helix (closed spiral on a circle).
% - Targets: r/R ~ alphaFS, R ~ muon reduced Compton wavelength.

clear; clc;

%% Constants / targets (packed in cfg)
cfg.alphaFS      = 1/137.035999084;      % fine-structure constant
cfg.lambdaC_bar  = 1.867594306e-15;      % muon reduced Compton wavelength [m] ≈ ħ/(m_μ c)
cfg.R_target     = cfg.lambdaC_bar;      % heuristic major-radius target
cfg.cL_over_cT_target = 1.6068;          % not used in geometry here (hook below)

%% Weights & knobs (same as electron script)
cfg.w_ratio   = 1.0;    % (r/R - alphaFS)^2
cfg.w_R       = 0.5;    % ((R - R_target)/R_target)^2
cfg.w_kvar    = 0.1;    % curvature CV^2
cfg.w_self    = 2.0;    % self-intersection penalty
cfg.w_balance = 1.0;    % placeholder "dynamic balance" residual
cfg.minSepFac = 1.5;    % penalize if min separation < minSepFac * r

cfg.Nsamp   = 2000;     % samples on curve
cfg.M_self  = 400;      % downsample for self-intersection check

nCandidates = 1:6;      % integer poloidal windings (closure)

%% Multi-start seeds (scaled to muon size)
R_guess_vec     = cfg.R_target * [0.6, 1.0, 1.6];
r_over_R_guess  = [0.5*cfg.alphaFS, cfg.alphaFS, 2*cfg.alphaFS];

best = struct('J',inf);

fprintf('Searching (muon) over n = %s ...\n', mat2str(nCandidates));
for n = nCandidates
    for R0 = R_guess_vec
        for ratio0 = r_over_R_guess
            r0  = max(1e-18, ratio0*R0);
            x0  = log([R0, r0]);
            obj = @(x) objectiveFromLogVars(x, n, cfg);
            opts = optimset('Display','off','TolX',1e-12,'TolFun',1e-12,...
                            'MaxFunEvals',3e3,'MaxIter',3e3);
            [xopt, J] = fminsearch(obj, x0, opts);

            if J < best.J
                [~, pack] = obj(xopt);
                best = pack; best.J = J; best.n = n;
                fprintf('  n=%d  R=%.3e  r=%.3e  r/R=%.5f  J=%.3e  minSep/r=%.2f\n', ...
                        n, best.R, best.r, best.r/best.R, J, best.minSep/best.r);
            end
        end
    end
end

%% Report best parameters
fprintf('\nBest (muon) parameters:\n');
fprintf('  n   = %d\n', best.n);
fprintf('  R   = %.6e m\n', best.R);
fprintf('  r   = %.6e m\n', best.r);
fprintf('  r/R = %.8f  (target alphaFS = %.8f)\n', best.r/best.R, cfg.alphaFS);
fprintf('  Objective J = %.3e\n', best.J);

%% Path + plots
[t, X, Y, Z] = toroidalHelix(best.R, best.r, best.n, cfg.Nsamp);
[kap, ~, s]  = curvatureAndTorsion(X,Y,Z);

figure('Color','w'); 
plot3(X, Y, Z, 'LineWidth', 1.5); grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('GDM muon-at-rest path (n=%d),  R=%.2e m,  r=%.2e m,  r/R=%.5f', ...
      best.n, best.R, best.r, best.r/best.R));
hold on; tt = linspace(0,2*pi,400);
plot3(best.R*cos(tt), best.R*sin(tt), 0*tt, '--', 'Color',[0.3 0.3 0.3]);

figure('Color','w');
plot(s - s(1), kap, 'LineWidth', 1.2);
xlabel('arclength s [m]'); ylabel('curvature \kappa [1/m]');
title('Curvature along the path'); grid on;

%% ---------------- helpers ----------------
function [J, pack] = objectiveFromLogVars(xlog, n, cfg)
    X = exp(xlog); R = X(1); r = X(2);

    if r >= R || R <= 0 || r <= 0
        J = 1e9 + 1e6*(r/R - cfg.alphaFS)^2; pack = []; return;
    end

    [~, Xc, Yc, Zc] = toroidalHelix(R, r, n, cfg.Nsamp);
    [kap, ~, ~]     = curvatureAndTorsion(Xc,Yc,Zc);
    kap = kap(~isnan(kap) & ~isinf(kap));
    if isempty(kap), J = 1e9; pack = []; return; end

    kmean = mean(kap); kvar = var(kap);
    cov2  = (kvar / max(kmean^2, eps));

    idx = round(linspace(1, cfg.Nsamp, cfg.M_self));
    P   = [Xc(idx), Yc(idx), Zc(idx)];
    minSep = minPairwiseDistance(P, max(3, round(0.01*cfg.M_self)));
    selfPenalty = max(0, (cfg.minSepFac*r - minSep) / (cfg.minSepFac*r) )^2;

    balance = balanceResidual(R, r, n, cfg);  % <<< put real GDM residual here

    JR   = ((R - cfg.R_target)/cfg.R_target)^2;
    Jrat = (r/R - cfg.alphaFS)^2;
    J = cfg.w_ratio*Jrat + cfg.w_R*JR + cfg.w_kvar*cov2 + ...
        cfg.w_self*selfPenalty + cfg.w_balance*balance;

    pack = struct('R',R,'r',r,'n',n,'J',J,'cov2',cov2,'minSep',minSep, ...
                  'JR',JR,'Jrat',Jrat,'selfPenalty',selfPenalty,'balance',balance);
end

function bal = balanceResidual(R, r, n, cfg)
    % Placeholder: replace with your muon-specific GDM stationarity checks
    rn = r/R;
    bal = 1e-3*(n-1)^2 + 1e-4*((rn/cfg.alphaFS - 1))^2;
end

function [t, X, Y, Z] = toroidalHelix(R, r, n, N)
    t = linspace(0, 2*pi, N+1)'; t(end) = [];
    ct = cos(t);  st = sin(t);
    cn = cos(n*t); sn = sin(n*t);
    X = (R + r.*cn).*ct;
    Y = (R + r.*cn).*st;
    Z =  r.*sn;
end

function [kappa, tau, s] = curvatureAndTorsion(X,Y,Z)
    N = numel(X);
    t = linspace(0,1,N)';  dt = mean(diff(t));
    r  = [X Y Z];
    rp = gradient(r, dt);
    rpp= gradient(rp, dt);
    rppp=gradient(rpp, dt);

    cross12 = cross(rp, rpp, 2);
    numK = sqrt(sum(cross12.^2,2));
    denK = (sqrt(sum(rp.^2,2))).^3 + eps;
    kappa = numK ./ denK;

    numT = dot(cross12, rppp, 2);
    denT = sum(cross12.^2,2) + eps;
    tau  = numT ./ denT;

    ds = sqrt(sum(rp.^2,2)) * dt;
    s  = cumsum(ds);
end

function dmin = minPairwiseDistance(P, kskip)
    M = size(P,1); dmin = inf;
    for i = 1:M
        for j = i+kskip:M
            d = P(i,:) - P(j,:);
            dij = d*d.';   % squared distance
            if dij < dmin^2, dmin = sqrt(dij); end
        end
    end
    if isinf(dmin), dmin = 0; end
end
