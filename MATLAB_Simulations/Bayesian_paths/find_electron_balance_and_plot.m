%% find_electron_balance_and_plot.m  (fixed: no 'alpha' name clash)
% Searches for a smooth, closed-loop "electron-at-rest" path and plots it.

clear; clc;

%% Constants / targets (pack into cfg to avoid evalin)
cfg.alphaFS      = 1/137.035999084;     % fine-structure constant (≡ α)
cfg.lambdaC_bar  = 3.8615926796e-13;    % reduced Compton wavelength [m]
cfg.R_target     = cfg.lambdaC_bar;     % heuristic target for major radius
cfg.cL_over_cT_target = 1.6068;         % not used in geometry (hook available)

%% Weights & knobs
cfg.w_ratio   = 1.0;    % weight for (r/R - alphaFS)^2
cfg.w_R       = 0.5;    % weight for (R - R_target)^2 (relative)
cfg.w_kvar    = 0.1;    % weight for curvature CV^2
cfg.w_self    = 2.0;    % weight for self-intersection penalty
cfg.w_balance = 1.0;    % weight for (placeholder) GDM balance residual
cfg.minSepFac = 1.5;    % penalize if min separation < minSepFac * r

cfg.Nsamp   = 2000;     % samples along curve for curvature/smoothness
cfg.M_self  = 400;      % downsample for self-intersection checking

nCandidates = 1:6;      % integer "wobble" counts (toroidal helix frequency)

%% Multi-start seeds
R_guess_vec     = cfg.R_target * [0.6, 1.0, 1.6];
r_over_R_guess  = [0.5*cfg.alphaFS, cfg.alphaFS, 2*cfg.alphaFS];

best = struct('J',inf);

fprintf('Searching over n = %s ...\n', mat2str(nCandidates));
for n = nCandidates
    for R0 = R_guess_vec
        for ratio0 = r_over_R_guess
            r0  = max(1e-16, ratio0*R0);
            x0  = log([R0, r0]);                 % optimize in log-space
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
fprintf('\nBest solution:\n');
fprintf('  n   = %d\n', best.n);
fprintf('  R   = %.6e m\n', best.R);
fprintf('  r   = %.6e m\n', best.r);
fprintf('  r/R = %.8f  (target alphaFS = %.8f)\n', best.r/best.R, cfg.alphaFS);
fprintf('  Objective J = %.3e\n', best.J);

%% Generate full-resolution path and plot
[t, X, Y, Z] = toroidalHelix(best.R, best.r, best.n, cfg.Nsamp);
[kap, ~, s]  = curvatureAndTorsion(X,Y,Z);

figure('Color','w'); 
plot3(X, Y, Z, 'LineWidth', 1.5); grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('GDM electron-at-rest path (n=%d),  R=%.2e m,  r=%.2e m,  r/R=%.5f', ...
      best.n, best.R, best.r, best.r/best.R));
hold on; tt = linspace(0,2*pi,400);
plot3(best.R*cos(tt), best.R*sin(tt), 0*tt, '--', 'Color',[0.3 0.3 0.3]);

figure('Color','w');
plot(s - s(1), kap, 'LineWidth', 1.2);
xlabel('arclength s [m]'); ylabel('curvature \kappa [1/m]');
title('Curvature along the path'); grid on;

%% -------- helpers (no evalin; everything passed via cfg) --------

function [J, pack] = objectiveFromLogVars(xlog, n, cfg)
    X = exp(xlog); R = X(1); r = X(2);

    % validity
    if r >= R || R <= 0 || r <= 0
        J = 1e9 + 1e6*(r/R - cfg.alphaFS)^2; pack = [];
        return;
    end

    % curve + diagnostics
    [~, Xc, Yc, Zc] = toroidalHelix(R, r, n, cfg.Nsamp);
    [kap, ~, ~]     = curvatureAndTorsion(Xc,Yc,Zc);
    kap = kap(~isnan(kap) & ~isinf(kap));
    if isempty(kap), J = 1e9; pack = []; return; end

    % curvature CV^2
    kmean = mean(kap); kvar = var(kap);
    cov2  = (kvar / max(kmean^2, eps));

    % self-intersection penalty
    idx = round(linspace(1, cfg.Nsamp, cfg.M_self));
    P   = [Xc(idx), Yc(idx), Zc(idx)];
    minSep = minPairwiseDistance(P, max(3, round(0.01*cfg.M_self)));
    selfPenalty = max(0, (cfg.minSepFac*r - minSep) / (cfg.minSepFac*r) )^2;

    % placeholder dynamic balance residual (replace with real GDM residual)
    balance = balanceResidual(R, r, n, cfg);

    % objective
    JR   = ((R - cfg.R_target)/cfg.R_target)^2;
    Jrat = (r/R - cfg.alphaFS)^2;
    J = cfg.w_ratio*Jrat + cfg.w_R*JR + cfg.w_kvar*cov2 + ...
        cfg.w_self*selfPenalty + cfg.w_balance*balance;

    pack = struct('R',R,'r',r,'n',n,'J',J,'cov2',cov2,'minSep',minSep, ...
                  'JR',JR,'Jrat',Jrat,'selfPenalty',selfPenalty,'balance',balance);
end

function bal = balanceResidual(R, r, n, cfg)
    % ----- replace with your real GDM stationarity/no-radiation checks -----
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
