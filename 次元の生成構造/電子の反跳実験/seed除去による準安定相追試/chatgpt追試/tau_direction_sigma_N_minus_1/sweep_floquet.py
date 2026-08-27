# -*- coding: utf-8 -*-
"""Rotating-frame Jacobian at the parent relative equilibrium for N=3..16.
Classifies multipliers: unstable (|mu|>1), neutral (|mu|=1), decaying; reports angles / parent dphi."""
import importlib.util, math, json, sys, time
import numpy as np
HERE = __file__.rsplit("/", 1)[0]
spec = importlib.util.spec_from_file_location("eng", HERE + "/run_n_scaling_lowrank_v1_no_sigma_norm.py")
eng = importlib.util.module_from_spec(spec); spec.loader.exec_module(eng)
eng.progress = lambda msg: None
GAMMA = eng.GAMMA
NS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else list(range(3, 17))
EPS = 1e-6
out = {}
for n in NS:
    t0 = time.time()
    syslr = eng.LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n + 0)
    v, res, sig = eng.make_parent(syslr, rng, iters=1200, tol=1e-12)
    m = len(v)
    syslr.set_theta(np.angle(v)); Fv = syslr.cayley_step(v, None); phase = np.vdot(v, Fv); phase = phase / abs(phase)
    dphi = math.atan2(phase.imag, phase.real)
    def R(z):
        syslr.set_theta(np.angle(z)); return syslr.cayley_step(z, None) / phase
    defect = float(np.linalg.norm(R(v) - v))
    x0 = np.r_[v.real, v.imag]; J = np.zeros((2 * m, 2 * m))
    for j in range(2 * m):
        dx = np.zeros(2 * m); dx[j] = EPS
        zp = (x0 + dx)[:m] + 1j * (x0 + dx)[m:]; zm = (x0 - dx)[:m] + 1j * (x0 - dx)[m:]
        dy = (R(zp) - R(zm)) / (2 * EPS); J[:m, j] = dy.real; J[m:, j] = dy.imag
    ev = np.linalg.eigvals(J); mods = np.abs(ev); ang = np.angle(ev)
    unst = [(float(mods[i]), float(ang[i])) for i in np.argsort(-mods) if mods[i] > 1 + 1e-6]
    neut = sorted({round(float(abs(ang[i])) / dphi, 4) for i in range(len(ev)) if abs(mods[i] - 1) <= 1e-6})
    dec = sorted({(round(float(mods[i]), 5), round(float(abs(ang[i])) / dphi, 3)) for i in range(len(ev)) if mods[i] < 1 - 1e-6})
    sumlog = float(np.sum(np.log(mods[mods > 0])))
    res_n = {"N": n, "M": m, "parent_residual": float(res), "defect": defect, "parent_sigma2": [round(float(s * s), 6) for s in sig[:5]],
             "dphi": dphi, "n_unstable_real_dims": len(unst), "unstable": unst[:8],
             "unstable_angle_over_dphi": [round(a / dphi, 4) for _, a in unst[:8]],
             "neutral_angle_over_dphi": neut, "decaying_mod_angle_over_dphi": dec[:10], "sum_log_mod": sumlog, "sec": round(time.time() - t0, 1)}
    out[n] = res_n
    print(f"N={n:2d} M={m:3d} σ²={res_n['parent_sigma2'][:3]} | unstable dims={len(unst)} |μ|={[round(u[0],6) for u in unst[:4]]} angle/Δφ={res_n['unstable_angle_over_dphi'][:4]} | neutral angle/Δφ={neut} | Σln|μ|={sumlog:+.4f} | defect={defect:.1e} | {res_n['sec']}s", flush=True)
    json.dump(out, open(HERE + "/sweep_floquet.json", "w"), indent=1)
print("DONE")
