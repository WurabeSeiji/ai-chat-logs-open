# -*- coding: utf-8 -*-
"""Verification of the sigma = N-1 theorem and its robustness (same engine, no new physics).
 (1) seed dependence of late sigma_eff^2      (2) sigma_eff^2(t) trajectories
 (3) N=3 long run (no vertex closure possible) (4) theorem ingredients on final states
Output: results/verify_theorem.json
"""
import importlib.util, math, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, "results"); os.makedirs(R, exist_ok=True)
spec = importlib.util.spec_from_file_location("eng", os.path.join(HERE, "run_n_scaling_lowrank_v1_no_sigma_norm.py"))
eng = importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress = lambda m: None
G = eng.GAMMA

def evolve(n, seed, steps, record_every=None):
    s = eng.LowRankSystem(n); rng = np.random.default_rng(40260722 + 1000 * n + seed)
    v, res, sig = eng.make_parent(s, rng, iters=1200, tol=1e-12); Z = v.copy(); wp = rng.normal(size=s.m); traj = []
    for t in range(steps):
        s.set_theta(np.angle(Z)); se, wp = s.sigma_max_power(wp); Zn = s.cayley_step(Z, se)
        if record_every and t % record_every == 0:
            d = np.angle(np.vdot(Z, Zn)); traj.append([t, (math.tan(d / 2) / G) ** 2])
        Z = Zn
    inc = []
    for t in range(200):
        s.set_theta(np.angle(Z)); se, wp = s.sigma_max_power(wp); Zn = s.cayley_step(Z, se); inc.append(np.angle(np.vdot(Z, Zn))); Z = Zn
    dphi = float(np.mean(inc))
    return s, Z, (math.tan(dphi / 2) / G) ** 2, [float(x * x) for x in sig[:3]], traj

out = {"gamma": G, "seed_dependence": [], "trajectories": {}, "N3_long": [], "theorem_check": []}
for n in (5, 6, 8, 11):
    for seed in (1, 2):
        _, _, se2, ps, _ = evolve(n, seed, 5000)
        out["seed_dependence"].append({"N": n, "seed": seed, "parent_sigma2": ps, "late_sigma_eff2": se2, "target": (n - 1) ** 2})
        print(f"N={n} seed={seed}: late sigma_eff^2={se2:.6f} target={(n-1)**2}")
for n in (5, 8):
    _, _, _, _, traj = evolve(n, 0, 2001, record_every=100); out["trajectories"][str(n)] = traj
    print(f"N={n} trajectory:", [(t, round(v, 3)) for t, v in traj[::4]])
for steps in (20000, 60000):
    _, _, se2, ps, _ = evolve(3, 0, steps); out["N3_long"].append({"steps": steps, "late_sigma_eff2": se2}); print(f"N=3 steps={steps}: sigma_eff^2={se2:.6f}")
for n in (4, 5, 6, 8, 12):
    s, Z, _, _, _ = evolve(n, 0, 5000)
    ea, eb = s.ea, s.eb; th = np.angle(Z); r = np.abs(Z)
    vc = max(abs(sum(Z[e] ** 2 for e in range(s.m) if ea[e] == i or eb[e] == i)) for i in range(n))
    S2 = []; sigf = []
    for e in range(s.m):
        nb = [f for f in range(s.m) if f != e and (ea[f] in (ea[e], eb[e]) or eb[f] in (ea[e], eb[e]))]
        S2.append(complex(sum(np.exp(2j * (th[f] - th[e])) for f in nb)))
        sigf.append(len(nb) / 2 - 0.5 * sum(np.cos(2 * (th[f] - th[e])) for f in nb))
    S2 = np.array(S2); s.set_theta(th); kv = s.kmatvec(Z)
    mu = float(np.real(np.conj(Z) @ (1j * kv))) / float(np.vdot(Z, Z).real); resid = float(np.linalg.norm(1j * kv - mu * Z) / np.linalg.norm(Z))
    rec = {"N": n, "r_max_over_min": float(r.max() / r.min()), "vertex_closure_max_abs": float(vc),
           "sum_exp2iDelta_re_range": [float(S2.real.min()), float(S2.real.max())], "sum_exp2iDelta_im_max_abs": float(abs(S2.imag).max()),
           "sigma_formula_range": [float(min(sigf)), float(max(sigf))], "eigen_sigma": -mu, "eigenmode_residual": resid, "N_minus_1": n - 1}
    out["theorem_check"].append(rec)
    print(f"N={n}: sum e^(2iΔ) re∈[{rec['sum_exp2iDelta_re_range'][0]:.6f},{rec['sum_exp2iDelta_re_range'][1]:.6f}] | sigma formula∈[{rec['sigma_formula_range'][0]:.6f},{rec['sigma_formula_range'][1]:.6f}] | eigen sigma={-mu:.6f} resid={resid:.1e}")
json.dump(out, open(os.path.join(R, "verify_theorem.json"), "w"), indent=1); print("DONE")
