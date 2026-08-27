# -*- coding: utf-8 -*-
"""Does anything depend on the clock?  Sweep L (gamma = tan(pi/L)) at equal natural time, and replace the Cayley step
by the exact exponential map exp(2*gamma*K) (same state-dependent K).  Same engine, same seeds.
Readouts: late sigma_eff^2, rigidity, z^2 class structure (N=5), A-B phase offset delta, onset in natural time, late growth rate per natural time.
Output: results/sweep_L_and_expmap.json
"""
import importlib.util, math, json, os, sys
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
    """z <- exp(2*gamma*K) z  (K real antisymmetric -> iK Hermitian)."""
    K = dense_K(s, len(z)); w, V = np.linalg.eigh(1j * K)          # iK = V diag(w) V^H, K = -i V diag(w) V^H
    return V @ (np.exp(-2j * gamma * w) * (V.conj().T @ z))

def classes_N5(Z):
    """cluster z^2 values (10 edges) -> class sizes, and A/B family phase offset delta/pi (mod 1/2)."""
    q = Z * Z; ph = np.angle(q) / (2 * math.pi) % 0.5               # z^2 phase / 2pi, mod 1/2  (sign flips identified)
    r2 = np.abs(Z) ** 2
    vals = sorted(ph); groups = []
    for v in vals:
        if groups and min(abs(v - groups[-1][-1]), 0.5 - abs(v - groups[-1][-1])) < 2e-3: groups[-1].append(v)
        else: groups.append([v])
    sizes = sorted([len(g) for g in groups], reverse=True)
    centers = sorted([float(np.mean(g)) for g in groups])
    delta = None
    if len(centers) == 2:
        d = abs(centers[0] - centers[1]); delta = min(d, 0.5 - d)          # in units of 2pi for z^2 == units of pi for z
    return sizes, [round(c, 6) for c in centers], delta, float(r2.max() / r2.min())

def run(n, L, steps, seed=0, integrator="cayley"):
    gamma = math.tan(math.pi / L); eng.GAMMA = gamma
    s = eng.LowRankSystem(n); rng = np.random.default_rng(40260722 + 1000 * n + seed)
    v, res, sig = eng.make_parent(s, rng, iters=1200, tol=1e-12); Z = v.copy(); wp = rng.normal(size=s.m)
    p = v.real / np.linalg.norm(v.real); qv = v.imag - (v.imag @ p) * p; qv /= np.linalg.norm(qv)
    onset = None; f_series = []
    for t in range(steps):
        hpar = abs(p @ Z) ** 2 + abs(qv @ Z) ** 2; f = max(0.0, 1 - hpar / float(np.vdot(Z, Z).real))
        f_series.append(f)
        if onset is None and f >= 1e-8: onset = t
        s.set_theta(np.angle(Z))
        if integrator == "cayley":
            se, wp = s.sigma_max_power(wp); Z = s.cayley_step(Z, se)
        else:
            Z = expmap_step(s, Z, gamma)
    inc = []
    for t in range(200):
        s.set_theta(np.angle(Z)); Zn = s.cayley_step(Z, None) if integrator == "cayley" else expmap_step(s, Z, gamma)
        inc.append(np.angle(np.vdot(Z, Zn))); Z = Zn
    dphi = float(np.mean(inc))
    sig_eff = (math.tan(dphi / 2) / gamma) if integrator == "cayley" else dphi / (2 * gamma)
    fs = np.array(f_series); mask = (fs > 1e-10) & (fs < 1e-3); idx = np.where(mask)[0]
    rate_step = float(np.polyfit(idx, np.log(fs[idx]), 1)[0]) if len(idx) > 10 else float("nan")
    out = {"N": n, "L": L, "integrator": integrator, "steps": steps, "natural_time_total": 2 * gamma * steps,
           "parent_sigma2": [round(float(x * x), 6) for x in sig[:3]], "late_sigma_eff": sig_eff, "late_sigma_eff2": sig_eff ** 2,
           "onset_step": onset, "onset_natural_time": None if onset is None else 2 * gamma * onset,
           "growth_rate_per_step": rate_step, "growth_rate_per_natural_time": rate_step / (2 * gamma)}
    if n == 5:
        sizes, centers, delta, eq = classes_N5(Z); out.update({"class_sizes": sizes, "class_centers_z2_over_2pi": centers, "delta_over_pi": delta, "equip_max_over_min": eq})
    return out

results = []
for L in (144, 288, 576, 1152):
    for n in (5, 8):
        r = run(n, L, steps=int(5000 * L / 144)); results.append(r)
        extra = f" classes={r['class_sizes']} δ/π={r['delta_over_pi']:.6f}" if n == 5 and r.get("delta_over_pi") is not None else (f" classes={r.get('class_sizes')}" if n == 5 else "")
        print(f"cayley L={L:5d} N={n}: σ_eff²={r['late_sigma_eff2']:.6f} onset s={r['onset_natural_time']:.4f} rate/s={r['growth_rate_per_natural_time']:.5f}{extra}", flush=True)
for L in (144, 288):
    r = run(5, L, steps=int(5000 * L / 144), integrator="expmap"); results.append(r)
    print(f"expmap L={L:5d} N=5: σ_eff²={r['late_sigma_eff2']:.6f} onset s={r['onset_natural_time']:.4f} rate/s={r['growth_rate_per_natural_time']:.5f} classes={r['class_sizes']} δ/π={r['delta_over_pi']:.6f}", flush=True)
json.dump(results, open(os.path.join(R, "sweep_L_and_expmap.json"), "w"), indent=1); print("DONE")
