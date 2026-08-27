# -*- coding: utf-8 -*-
"""Two independent controls for the 'parent normalization removed' experiment (no re-run of the reference script; its
step functions are re-implemented verbatim here).
 (A) Scale invariance: with the NORMALIZED parent (normalized_engine = published engine, same seed 40260721+1000*5),
     the phase-only baseline must equal the unnormalized trajectory divided by c = |v_unnorm|, and the amplitude-aware
     treatment must equal the unnormalized trajectory with time rescaled by c^2 (K ~ |z|^2).  Also checks the residual
     identity res = sigma*|v|*|1-|v|^2| that explains the reported parent_residual 0.487.
 (B) Engine cross-check: the reference's phase-only exp((2pi/124)K) step is the published low-rank engine's K with the
     exact exponential map; compare onset, growth slope, final sigma1 and PR against the reference summary.
Output: results/scale_and_engine_check.json"""
import os, sys, math, json, importlib.util, numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(HERE, "results"); os.makedirs(R, exist_ok=True)
REF = os.path.join(os.path.dirname(HERE), "N5_linear124_all3fix_seedless_parentnorm_removed_20260828", "data")
def load_engine(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "program", name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); m.progress = lambda s: None; return m
eng_un, eng_n = load_engine("original_engine"), load_engine("normalized_engine")
N, STEPS, SEED, L = 5, 5000, 0, 124

# ---- reference step functions (verbatim logic) ----
def adjacency_mask(s):
    m = s.m; A = np.zeros((m, m))
    for i in range(m):
        share = (s.ea == s.ea[i]) | (s.ea == s.eb[i]) | (s.eb == s.ea[i]) | (s.eb == s.eb[i]); share[i] = False; A[i, share] = 1.0
    return A
def K_phase_only(s, z):
    m = s.m; K = np.zeros((m, m)); th = np.angle(z)
    for i in range(m):
        share = (s.ea == s.ea[i]) | (s.ea == s.eb[i]) | (s.eb == s.ea[i]) | (s.eb == s.eb[i]); share[i] = False; K[i, share] = np.sin(th[share] - th[i])
    return K
def K_amp_aware(s, z):
    K = adjacency_mask(s) * np.imag(np.conj(z)[:, None] * z[None, :]); np.fill_diagonal(K, 0.0); return K
def exp_step(K, z, angle):
    w, V = np.linalg.eigh(1j * K); return V @ (np.exp(-1j * angle * w) * (V.conj().T @ z))
def run(s, v, mode, steps=STEPS, angle=2 * math.pi / L):
    Z = v.copy(); p = v.real / np.linalg.norm(v.real); q = v.imag - (v.imag @ p) * p; q /= np.linalg.norm(q)
    Zs = np.empty((steps + 1, s.m), complex); H = np.empty((steps + 1, 3))
    for t in range(steps + 1):
        Zs[t] = Z; hp = abs(p @ Z) ** 2 + abs(q @ Z) ** 2; ht = float(np.vdot(Z, Z).real); H[t] = (hp, ht - hp, ht)
        if t == steps: break
        K = K_phase_only(s, Z) if mode == "phase" else K_amp_aware(s, Z); Z = exp_step(K, Z, angle)
    return Zs, H
def parents():
    su = eng_un.LowRankSystem(N); ru = np.random.default_rng(40260721 + 1000 * N + SEED); vu, resu, sigu = eng_un.make_parent(su, ru)
    sn = eng_n.LowRankSystem(N); rn = np.random.default_rng(40260721 + 1000 * N + SEED); vn, resn, sign = eng_n.make_parent(sn, rn)
    return (su, vu, resu, sigu), (sn, vn, resn, sign)

out = {}
(su, vu, resu, sigu), (sn, vn, resn, sign) = parents()
c = float(np.linalg.norm(vu)); out["c_norm_unnormalized_parent"] = c; out["c2"] = c * c
out["residual_identity"] = {"reported_residual": float(resu), "sigma*|v|*|1-|v|^2|": float(sigu[0] * c * abs(1 - c * c)), "normalized_parent_residual": float(resn)}
# same direction?  (up to global phase)
ov = abs(np.vdot(vn, vu)) / (np.linalg.norm(vn) * np.linalg.norm(vu)); out["parent_direction_overlap"] = float(ov)
print(f"c=|v_un|={c:.12f}, c^2={c*c:.12f}; residual reported {resu:.6f} vs identity {sigu[0]*c*abs(1-c*c):.6f}; normalized residual {resn:.1e}; direction overlap {ov:.12f}")
# (A1) baseline scale invariance
Zu, Hu = run(su, vu, "phase"); Zn, Hn = run(sn, vn, "phase")
# align global phase of parents (make_parent may return v up to a phase factor); compare |Z| and H
ph = np.vdot(vn, vu) / abs(np.vdot(vn, vu))
dZ = np.abs(Zu - c * ph * Zn).max(); dH = np.abs(Hu - c * c * Hn).max()
out["baseline_scale_check"] = {"max|Z_un - c e^{i a} Z_norm|": float(dZ), "max|H_un - c^2 H_norm|": float(dH), "final_sigma1_un": None}
print(f"(A1) baseline: max|Z_un − c·Z_norm| = {dZ:.2e}, max|H_un − c²H_norm| = {dH:.2e}")
# (A2) treatment time rescaling: unnormalized with angle a  == normalized with angle a*c^2  (K_un = c^2 K_norm at t=0 and stays so)
Zt_u, Ht_u = run(su, vu, "amp"); Zt_n, Ht_n = run(sn, vn, "amp", angle=2 * math.pi / L * c * c)
dZt = np.abs(Zt_u - c * ph * Zt_n).max(); dHt = np.abs(Ht_u - c * c * Ht_n).max()
out["treatment_time_rescale_check"] = {"angle_normalized_run": 2 * math.pi / L * c * c, "max|Z_un(t) - c e^{i a} Z_norm(t; angle*c^2)|": float(dZt), "max|H diff|": float(dHt)}
print(f"(A2) treatment: unnormalized (angle 2π/124) vs normalized with angle·c² : max|ΔZ| = {dZt:.2e}, max|ΔH| = {dHt:.2e}")
# (B) engine cross-check: published low-rank kmatvec vs reference dense K (phase-only), and exp map at gamma = pi/124
sn.set_theta(np.angle(vn)); e = np.eye(sn.m); Kdense_engine = np.column_stack([sn.kmatvec(e[:, j]) for j in range(sn.m)])
out["engine_K_vs_reference_K_maxdiff"] = float(np.abs(Kdense_engine - K_phase_only(sn, vn)).max())
f = Hn[:, 1] / Hn[:, 2]; onset = int(np.argmax(f > 0.05)); mask = (Hn[:, 1] > 1e-10) & (Hn[:, 1] < 1e-3); idx = np.where(mask)[0]
slope = float(np.polyfit(idx, np.log(Hn[idx, 1]), 1)[0]) if len(idx) > 3 else None
sn.set_theta(np.angle(Zn[-1])); sig_final = sn.sigma_spectrum(); a2 = np.abs(Zn[-1]) ** 2; pr = float(a2.sum() ** 2 / (a2 ** 2).sum() / sn.m)
ref = json.load(open(os.path.join(REF, "summary.json")))["baseline"]
out["engine_crosscheck_baseline"] = {"onset_frac_gt_0.05": onset, "ref_onset": ref["onset_Hperp_fraction_gt_0.05"], "growth_slope_per_step": slope,
                                     "ref_growth_slope": ref["growth_fit"]["slope_ln_Hperp_per_step"], "final_sigma1": float(sig_final[0]), "final_PR_over_M": pr,
                                     "ref_final_PR_over_M": ref["final_PR_over_M"], "natural_time_per_step": 2 * math.pi / L, "growth_per_natural_time": slope / (2 * math.pi / L) if slope else None}
print(f"(B) engine K vs reference K max diff {out['engine_K_vs_reference_K_maxdiff']:.1e}; normalized baseline: onset {onset} (ref {ref['onset_Hperp_fraction_gt_0.05']}), slope {slope:.6f} (ref {ref['growth_fit']['slope_ln_Hperp_per_step']:.6f}), final σ1 {sig_final[0]:.6f}, PR/M {pr:.6f} (ref {ref['final_PR_over_M']:.6f})")
json.dump(out, open(os.path.join(R, "scale_and_engine_check.json"), "w"), indent=1); print("DONE")
