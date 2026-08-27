# -*- coding: utf-8 -*-
"""Is the slow ordering stage (equipartition + class formation) driven by the discretization?
Run N=5 and N=8 for L = 144, 288, 576, 1152 up to the same natural time s_max, and record in natural time:
  - onset (f >= 1e-8), 99% saturation of H_perp, equipartition error |r2max/r2min - 1| thresholds, sigma_eff^2 thresholds,
  - N=5: within-class spread of z^2 phases (4-class error) thresholds.
Output: results/ordering_vs_L.json
"""
import importlib.util, math, json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, "results"); os.makedirs(R, exist_ok=True)
spec = importlib.util.spec_from_file_location("eng", os.path.join(HERE, "run_n_scaling_lowrank_v1_no_sigma_norm.py"))
eng = importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress = lambda m: None
S_MAX = 900.0
THR_EQ = [1e-2, 1e-4, 1e-6]; THR_CLS = [1e-2, 1e-3, 1e-4]; THR_SIG = [1e-1, 1e-2, 1e-3, 1e-4]

def class_error_N5(Z):
    ph = (np.angle(Z * Z) / (2 * math.pi)) % 0.5
    # fixed class assignment from the paper: A={12,13,45}={0,1,9}, A-={14,15,23}={2,3,4}, B+={24,35}={5,8}, B-={25,34}={6,7}
    groups = [[0, 1, 9], [2, 3, 4], [5, 8], [6, 7]]
    err = 0.0
    for g in groups:
        v = ph[g]; c = np.angle(np.mean(np.exp(2j * math.pi * v / 0.5))) / (2 * math.pi) * 0.5
        d = np.abs(((v - c) + 0.25) % 0.5 - 0.25); err = max(err, float(d.max()))
    return err

def run(n, L):
    gamma = math.tan(math.pi / L); eng.GAMMA = gamma; ds = 2 * gamma
    s = eng.LowRankSystem(n); rng = np.random.default_rng(40260722 + 1000 * n)
    v, res, sig = eng.make_parent(s, rng, iters=1200, tol=1e-12); Z = v.copy(); wp = rng.normal(size=s.m)
    p = v.real / np.linalg.norm(v.real); qv = v.imag - (v.imag @ p) * p; qv /= np.linalg.norm(qv)
    steps = int(S_MAX / ds); target = (n - 1) ** 2
    rec = {"N": n, "L": L, "gamma": gamma, "steps": steps, "s_max": steps * ds, "onset_s": None, "sat99_s": None,
           "equip_s": {str(t): None for t in THR_EQ}, "sigma_s": {str(t): None for t in THR_SIG}, "class_s": {str(t): None for t in THR_CLS},
           "final_sigma_eff2": None, "final_equip_err": None, "final_class_err": None}
    hperp_final_est = None; Zprev = None
    for t in range(steps):
        hpar = abs(p @ Z) ** 2 + abs(qv @ Z) ** 2; ht = float(np.vdot(Z, Z).real); f = max(0.0, 1 - hpar / ht)
        if rec["onset_s"] is None and f >= 1e-8: rec["onset_s"] = t * ds
        s.set_theta(np.angle(Z)); se, wp = s.sigma_max_power(wp); Zn = s.cayley_step(Z, se)
        if t % 20 == 0:
            dphi = float(np.angle(np.vdot(Z, Zn))); sig2 = (math.tan(dphi / 2) / gamma) ** 2
            r2 = np.abs(Z) ** 2; eq = float(r2.max() / r2.min() - 1)
            for th in THR_SIG:
                if rec["sigma_s"][str(th)] is None and abs(sig2 - target) < th * target: rec["sigma_s"][str(th)] = t * ds
            for th in THR_EQ:
                if rec["equip_s"][str(th)] is None and eq < th: rec["equip_s"][str(th)] = t * ds
            if n == 5:
                ce = class_error_N5(Z)
                for th in THR_CLS:
                    if rec["class_s"][str(th)] is None and ce < th: rec["class_s"][str(th)] = t * ds
        Z = Zn
    dphi = float(np.angle(np.vdot(Zprev if Zprev is not None else Z, Z))) if False else None
    s.set_theta(np.angle(Z)); Zn = s.cayley_step(Z, None); dphi = float(np.angle(np.vdot(Z, Zn)))
    rec["final_sigma_eff2"] = (math.tan(dphi / 2) / gamma) ** 2; r2 = np.abs(Z) ** 2; rec["final_equip_err"] = float(r2.max() / r2.min() - 1)
    if n == 5: rec["final_class_err"] = class_error_N5(Z)
    return rec

out = []
for n in (5, 8):
    for L in (144, 288, 576, 1152):
        r = run(n, L); out.append(r)
        cls = f" class<1e-3 at s={r['class_s']['0.001']}" if n == 5 else ""
        print(f"N={n} L={L:4d}: onset s={r['onset_s']:.3f} | σ² within 1e-2 at s={r['sigma_s']['0.01']} 1e-4 at s={r['sigma_s']['0.0001']} | equip<1e-4 at s={r['equip_s']['0.0001']}{cls} | final σ²={r['final_sigma_eff2']:.6f} equip_err={r['final_equip_err']:.1e}", flush=True)
json.dump(out, open(os.path.join(R, "ordering_vs_L.json"), "w"), indent=1); print("DONE")
