# -*- coding: utf-8 -*-
"""N=3..16: parent sigma^2 spectrum, late-time collective rotation (sigma_eff^2), rigidity, final-state K spectrum.
Same engine / seeds / GAMMA as the decompactification package (no new physics)."""
import importlib.util, math, json, sys, time
import numpy as np
HERE = __file__.rsplit("/", 1)[0]
spec = importlib.util.spec_from_file_location("eng", HERE + "/run_n_scaling_lowrank_v1_no_sigma_norm.py")
eng = importlib.util.module_from_spec(spec); spec.loader.exec_module(eng)
eng.progress = lambda msg: None
GAMMA = eng.GAMMA
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
NS = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else list(range(3, 17))
SEED = 0
out = {}
for n in NS:
    t0 = time.time()
    syslr = eng.LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n + SEED)
    v, residual, sig_parent = eng.make_parent(syslr, rng, iters=1200, tol=1e-12)
    Z = v.copy(); wp = rng.normal(size=syslr.m)
    inc_edges = []; inc_global = []
    for t in range(STEPS):
        syslr.set_theta(np.angle(Z))
        sig_est, wp = syslr.sigma_max_power(wp)
        Zn = syslr.cayley_step(Z, sig_est)
        if t >= STEPS - 500:
            inc_edges.append(np.angle(Zn / Z))
            inc_global.append(np.angle(np.vdot(Z, Zn)))
        Z = Zn
    inc_edges = np.array(inc_edges); inc_global = np.array(inc_global)
    edge_mean = inc_edges.mean(axis=0)
    dphi = float(inc_global.mean())
    sig_eff = math.tan(dphi / 2) / GAMMA
    syslr.set_theta(np.angle(Z)); sig_final = syslr.sigma_spectrum()
    a2 = np.abs(Z) ** 2; equip = float(a2.max() / a2.min())
    res = {
        "N": n, "M": syslr.m, "parent_residual": float(residual),
        "parent_sigma2": [round(float(s * s), 6) for s in sig_parent],
        "parent_dphi": float(2 * math.atan(GAMMA * sig_parent[0])),
        "late_dphi": dphi, "late_sigma_eff": sig_eff, "late_sigma_eff2": sig_eff * sig_eff,
        "late_edge_increment_spread": float(edge_mean.max() - edge_mean.min()),
        "late_time_std_max": float(inc_edges.std(axis=0).max()),
        "late_global_time_std": float(inc_global.std()),
        "final_K_sigma2": [round(float(s * s), 6) for s in sig_final[:6]],
        "final_equip_max_over_min": equip,
        "steps": STEPS, "sec": round(time.time() - t0, 1),
    }
    out[n] = res
    print(f"N={n:2d} M={syslr.m:3d} parent σ²={res['parent_sigma2'][:4]} | late Δφ={dphi:.7f} σ_eff²={sig_eff*sig_eff:.6f} | rigid spread={res['late_edge_increment_spread']:.1e} tstd={res['late_time_std_max']:.1e} | final K σ²={res['final_K_sigma2'][:3]} | equip max/min={equip:.4f} | {res['sec']}s", flush=True)
    json.dump(out, open(HERE + "/sweep_sigma2.json", "w"), indent=1)
print("DONE")
