# joint_fit_alpha_G.py

import numpy as np, emcee, json, time

# --- Targets (update when you have your own electron loop result)
alpha_obs = 1/137.035999084       # dimensionless
sigma_alpha = 5e-9                 # loose prior width; adjust as you like

G_obs = 6.67430e-11               # SI
sigma_G = 1.5e-14                 # ~2e-4 relative; loosen/tighten as needed

# --- Model (mirror the MATLAB stub formulas)
def xi_star_from_medium(K, chi, xi0):
    return chi/(4*np.pi*K) + xi0  # REPLACE when you have real mapping

def G_from_medium(K, chi, a_map):
    return a_map * (chi / K)      # matches MATLAB stub

# --- Priors: broad log-uniform for positive scales; normal for xi0 near 0
def log_prior(theta):
    # theta = [logK, logChi, logAmap, xi0]
    logK, logChi, logA, xi0 = theta
    K, chi, a_map = np.exp(logK), np.exp(logChi), np.exp(logA)

    # broad bounds to keep sampler sane
    if not (-20 < logK < 20 and -20 < logChi < 20 and -40 < logA < 10 and abs(xi0) < 0.1):
        return -np.inf

    # mild Gaussian prior on xi0 centered at 0
    lp = -0.5*(xi0/0.02)**2
    return lp

def log_likelihood(theta):
    logK, logChi, logA, xi0 = theta
    K, chi, a_map = np.exp(logK), np.exp(logChi), np.exp(logA)

    xi_star = xi_star_from_medium(K, chi, xi0)
    G_pred  = G_from_medium(K, chi, a_map)

    # Likelihoods (independent Gaussians)
    ll_alpha = -0.5*((xi_star - alpha_obs)/sigma_alpha)**2
    ll_G     = -0.5*((G_pred  - G_obs   )/sigma_G   )**2
    return ll_alpha + ll_G

def log_posterior(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp): return -np.inf
    return lp + log_likelihood(theta)

# --- Initialize walkers near a reasonable scale
rng = np.random.default_rng(42)
ndim, nwalkers = 4, 32
p0 = []
for _ in range(nwalkers):
    logK   = rng.normal(0.0, 2.0)      # K ~ e^N(0,2)
    logChi = rng.normal(-2.0, 2.0)     # chi smaller scale
    logA   = np.log(G_obs) - (logChi - logK) + rng.normal(0, 1.0) # make G roughly right
    xi0    = rng.normal(0.0, 0.01)
    p0.append([logK, logChi, logA, xi0])
p0 = np.array(p0)

# --- Run emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior)
t0 = time.time()
sampler.run_mcmc(p0, 2500, progress=True)
print(f"Sampling took {time.time()-t0:.1f}s")

flat = sampler.get_chain(discard=800, thin=5, flat=True)
logp = sampler.get_log_prob(discard=800, thin=5, flat=True)

# Posterior summaries
mean = flat.mean(0); std = flat.std(0)
names = ["logK","logChi","logAmap","xi0"]
print("Posterior means:")
for n,m,s in zip(names, mean, std):
    print(f"  {n:8s} = {m: .4f} ± {s:.4f}")

# Derived posteriors
K, chi, a_map = np.exp(mean[0]), np.exp(mean[1]), np.exp(mean[2])
xi_star = xi_star_from_medium(K, chi, mean[3])
G_pred  = G_from_medium(K, chi, a_map)
print(f"Derived:  xi_star ≈ {xi_star:.12f}  (target α = {alpha_obs:.12f})")
print(f"          G_pred  ≈ {G_pred:.6e}  (target G = {G_obs:.6e})")

# Save
out = {
  "posterior_means": dict(zip(names, map(float,mean))),
  "posterior_stds" : dict(zip(names, map(float,std))),
  "xi_star"        : float(xi_star),
  "G_pred"         : float(G_pred)
}
with open("joint_fit_alpha_G_results.json","w") as f:
    json.dump(out,f,indent=2)
print("Wrote joint_fit_alpha_G_results.json")
