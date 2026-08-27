# -*- coding: utf-8 -*-
"""Is the local volume contraction (sum ln|mu| < 0) physical or a discretization effect?
Rotating-frame Jacobian at the parent for L = 144, 288, 576, 1152 (gamma = tan(pi/L)); report per-step and per-natural-time rates:
  growth exponent of the dominant unstable mode  ln|mu_1| / (2 gamma)   -> should converge (physical)
  sum ln|mu_i| / (2 gamma)                        -> if it goes to 0 as gamma -> 0, contraction is numerical
Also the exact exponential-map integrator at L=144 for comparison.
Output: results/floquet_vs_L.json
"""
import importlib.util, math, json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, "results"); os.makedirs(R, exist_ok=True)
spec = importlib.util.spec_from_file_location("eng", os.path.join(HERE, "run_n_scaling_lowrank_v1_no_sigma_norm.py"))
eng = importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress = lambda m: None

def dense_K(s, m):
    K = np.zeros((m, m))
    for j in range(m):
        e = np.zeros(m); e[j] = 1.0; K[:, j] = s.kmatvec(e)
    return K
def expmap_step(s, z, gamma):
    K = dense_K(s, len(z)); w, V = np.linalg.eigh(1j * K); return V @ (np.exp(-2j * gamma * w) * (V.conj().T @ z))

def floquet(n, L, integrator="cayley", eps=1e-6):
    gamma = math.tan(math.pi / L); eng.GAMMA = gamma
    s = eng.LowRankSystem(n); rng = np.random.default_rng(40260722 + 1000 * n)
    v, res, sig = eng.make_parent(s, rng, iters=1200, tol=1e-12); m = len(v)
    def F(z):
        s.set_theta(np.angle(z)); return s.cayley_step(z, None) if integrator == "cayley" else expmap_step(s, z, gamma)
    Fv = F(v); phase = np.vdot(v, Fv); phase /= abs(phase)
    Rm = lambda z: F(z) / phase
    x0 = np.r_[v.real, v.imag]; J = np.zeros((2 * m, 2 * m))
    for j in range(2 * m):
        dx = np.zeros(2 * m); dx[j] = eps
        zp = (x0 + dx)[:m] + 1j * (x0 + dx)[m:]; zm = (x0 - dx)[:m] + 1j * (x0 - dx)[m:]
        dy = (Rm(zp) - Rm(zm)) / (2 * eps); J[:m, j] = dy.real; J[m:, j] = dy.imag
    ev = np.linalg.eigvals(J); mods = np.abs(ev); ang = np.angle(ev); ds = 2 * gamma
    i1 = int(np.argmax(mods)); dphi = math.atan2(phase.imag, phase.real)
    neutral = sorted({round(float(abs(ang[i])) / dphi, 3) for i in range(len(ev)) if abs(mods[i] - 1) <= 1e-6})
    return {"N": n, "L": L, "integrator": integrator, "gamma": gamma, "ds": ds, "parent_sigma_max": float(sig[0]),
            "dphi": dphi, "dphi_over_ds_sigma": dphi / (ds * float(sig[0])),
            "mu1": float(mods[i1]), "lambda1_per_step": math.log(float(mods[i1])), "lambda1_per_s": math.log(float(mods[i1])) / ds,
            "n_unstable": int(np.sum(mods > 1 + 1e-6)), "sum_log_mod_per_step": float(np.sum(np.log(mods[mods > 0]))),
            "sum_log_mod_per_s": float(np.sum(np.log(mods[mods > 0]))) / ds, "neutral_angle_over_dphi": neutral}

out = []
for n in (5, 8):
    for L in (144, 288, 576, 1152):
        r = floquet(n, L); out.append(r)
        print(f"cayley N={n} L={L:4d}: λ1/s={r['lambda1_per_s']:.5f}  Σln|μ|/s={r['sum_log_mod_per_s']:+.5f}  (per step {r['sum_log_mod_per_step']:+.5f})  n_unst={r['n_unstable']}  Δφ/(ds·σ)={r['dphi_over_ds_sigma']:.5f}  neutral/Δφ={r['neutral_angle_over_dphi']}", flush=True)
for n in (5, 8):
    r = floquet(n, 144, "expmap"); out.append(r)
    print(f"expmap N={n} L= 144: λ1/s={r['lambda1_per_s']:.5f}  Σln|μ|/s={r['sum_log_mod_per_s']:+.5f}  n_unst={r['n_unstable']}  Δφ/(ds·σ)={r['dphi_over_ds_sigma']:.5f}  neutral/Δφ={r['neutral_angle_over_dphi']}", flush=True)
json.dump(out, open(os.path.join(R, "floquet_vs_L.json"), "w"), indent=1); print("DONE")
