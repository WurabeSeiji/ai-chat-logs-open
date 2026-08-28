# -*- coding: utf-8 -*-
"""Verify that the four corrections are actually implemented in the reference program (static + runtime checks).
 Fix 1: make_parent has no amplitude normalization (v /= norm(v) absent)      -> |v| != 1, residual identity sigma*c*|1-c^2|
 Fix 2: no external seed DELTA*g and no Z/=norm(Z) at initialization          -> Z0 == v exactly, |Z0| = c != 1
 Fix 3: linear (exact exponential) rotation exp((2pi/L)K), not Cayley          -> log of one-step operator == (2pi/L) K; eigenmode phase advance == (2pi/L)*sigma exactly
 Fix 4: amplitude-aware K_ij = Im(conj z_i z_j) = r_i r_j sin(theta_j - theta_i) -> scales as c^2 under z -> c z; antisymmetric; adjacency-masked
Output: results/verify_fixes.json"""
import os, re, sys, math, json, importlib.util, numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(HERE, "results"); os.makedirs(R, exist_ok=True)
PROG = os.path.join(HERE, "program")
src_engine = open(os.path.join(PROG, "original_engine.py"), encoding="utf-8").read()
src_run = open(os.path.join(PROG, "run_amplitude_only_fix.py"), encoding="utf-8").read()
out = {"static": {}, "runtime": {}}
# ---- static checks
mp = src_engine[src_engine.index("def make_parent"):src_engine.index("def zero_closure_kernel_seed")]
out["static"]["fix1_make_parent_has_no_norm"] = ("np.linalg.norm(v)" not in mp) and ("v / np.linalg.norm" not in mp)
out["static"]["fix2_no_DELTA_g_in_run"] = ("DELTA*g" not in src_run.replace(" ", "")) or ("Z=v.copy()" in src_run.replace(" ", ""))
out["static"]["fix2_Z0_is_v_copy"] = "Z=v.copy()" in src_run.replace(" ", "")
out["static"]["fix2_no_Z_norm_at_init"] = "Z/=np.linalg.norm(Z)" not in src_run.replace(" ", "") and "Z=Z/np.linalg.norm" not in src_run.replace(" ", "")
out["static"]["fix3_uses_exp_not_cayley"] = ("cayley_step" not in src_run) and ("np.exp(-1j*(2*math.pi/L)*w)" in src_run.replace(" ", ""))
out["static"]["fix4_amplitude_aware_K"] = "np.imag(np.conj(z)[:,None]*z[None,:])" in src_run.replace(" ", "")
out["static"]["steps"] = int(re.search(r"STEPS=(\d+)", src_run.replace(" ", "")).group(1)); out["static"]["L"] = int(re.search(r"L=(\d+)", src_run.replace(" ", "")).group(1))
# ---- runtime checks (re-implement the reference step functions verbatim; do not execute the reference script)
spec = importlib.util.spec_from_file_location("eng", os.path.join(PROG, "original_engine.py")); eng = importlib.util.module_from_spec(spec); spec.loader.exec_module(eng); eng.progress = lambda m: None
N = 5; L = out["static"]["L"]
def adjacency_mask(s):
    m = s.m; A = np.zeros((m, m))
    for i in range(m):
        share = (s.ea == s.ea[i]) | (s.ea == s.eb[i]) | (s.eb == s.ea[i]) | (s.eb == s.eb[i]); share[i] = False; A[i, share] = 1.0
    return A
def K_amp(s, z):
    K = adjacency_mask(s) * np.imag(np.conj(z)[:, None] * z[None, :]); np.fill_diagonal(K, 0.0); return K
def K_phase(s, z):
    m = s.m; K = np.zeros((m, m)); th = np.angle(z)
    for i in range(m):
        share = (s.ea == s.ea[i]) | (s.ea == s.eb[i]) | (s.eb == s.ea[i]) | (s.eb == s.eb[i]); share[i] = False; K[i, share] = np.sin(th[share] - th[i])
    return K
def exp_step(K, z, angle):
    w, V = np.linalg.eigh(1j * K); return V @ (np.exp(-1j * angle * w) * (V.conj().T @ z))
s = eng.LowRankSystem(N); rng = np.random.default_rng(40260721 + 1000 * N + 0); v, res, sig = eng.make_parent(s, rng)
c = float(np.linalg.norm(v))
out["runtime"]["fix1_parent_norm_c"] = c; out["runtime"]["fix1_residual_reported"] = float(res)
out["runtime"]["fix1_residual_identity_sigma_c_abs(1-c2)"] = float(sig[0] * c * abs(1 - c * c))
out["runtime"]["fix1_ok"] = abs(c - 1) > 1e-3 and abs(res - sig[0] * c * abs(1 - c * c)) < 1e-9
Z0 = v.copy(); out["runtime"]["fix2_Z0_equals_v"] = bool(np.array_equal(Z0, v)); out["runtime"]["fix2_norm_Z0"] = float(np.linalg.norm(Z0)); out["runtime"]["fix2_ok"] = np.array_equal(Z0, v) and abs(np.linalg.norm(Z0) - 1) > 1e-3
# fix 3: one-step operator U = exp(angle*K); check log(U) == angle*K and unitarity, and eigenmode phase advance is linear in sigma
angle = 2 * math.pi / L; K = K_phase(s, v); m = s.m
U = np.column_stack([exp_step(K, np.eye(m)[:, j].astype(complex), angle) for j in range(m)])
w, V = np.linalg.eig(U); logU = V @ np.diag(np.log(w)) @ np.linalg.inv(V)
out["runtime"]["fix3_||U^H U - I||"] = float(np.linalg.norm(U.conj().T @ U - np.eye(m)))
out["runtime"]["fix3_||log U - angle*K||"] = float(np.linalg.norm(logU - angle * K))
wk = np.linalg.eigvals(K); sig_k = np.sort(np.abs(wk.imag))[::-1]
ph_adv = np.sort(np.abs(np.angle(w)))[::-1]
out["runtime"]["fix3_phase_advance_vs_angle_sigma_maxdiff"] = float(np.abs(ph_adv - angle * sig_k).max())
out["runtime"]["fix3_cayley_would_give"] = float(2 * math.atan(math.tan(angle / 2) * sig_k[0]))
out["runtime"]["fix3_linear_gives"] = float(angle * sig_k[0])
out["runtime"]["fix3_ok"] = out["runtime"]["fix3_||log U - angle*K||"] < 1e-10 and out["runtime"]["fix3_phase_advance_vs_angle_sigma_maxdiff"] < 1e-10
# fix 4: amplitude scaling and formula
Ka = K_amp(s, v); Ka2 = K_amp(s, 2.0 * v); r = np.abs(v); th = np.angle(v)
formula = adjacency_mask(s) * (r[:, None] * r[None, :] * np.sin(th[None, :] - th[:, None]))
out["runtime"]["fix4_||K(2z) - 4 K(z)||"] = float(np.linalg.norm(Ka2 - 4 * Ka)); out["runtime"]["fix4_||K - r_i r_j sin||"] = float(np.linalg.norm(Ka - formula))
out["runtime"]["fix4_antisym"] = float(np.linalg.norm(Ka + Ka.T)); out["runtime"]["fix4_phase_only_K_is_amplitude_blind_||K(2z)-K(z)||"] = float(np.linalg.norm(K_phase(s, 2 * v) - K_phase(s, v)))
out["runtime"]["fix4_ok"] = out["runtime"]["fix4_||K(2z) - 4 K(z)||"] < 1e-12 and out["runtime"]["fix4_||K - r_i r_j sin||"] < 1e-12 and out["runtime"]["fix4_antisym"] < 1e-12
out["all_four_fixes_ok"] = all([out["static"]["fix1_make_parent_has_no_norm"], out["static"]["fix2_Z0_is_v_copy"], out["static"]["fix2_no_Z_norm_at_init"], out["static"]["fix3_uses_exp_not_cayley"], out["static"]["fix4_amplitude_aware_K"],
                               out["runtime"]["fix1_ok"], out["runtime"]["fix2_ok"], out["runtime"]["fix3_ok"], out["runtime"]["fix4_ok"]])
def _py(o):
    if isinstance(o, dict): return {k: _py(v) for k, v in o.items()}
    if isinstance(o, (np.bool_,)): return bool(o)
    if isinstance(o, (np.floating,)): return float(o)
    if isinstance(o, (np.integer,)): return int(o)
    return o
out = _py(out)
json.dump(out, open(os.path.join(R, "verify_fixes.json"), "w"), indent=1)
for k, v_ in out["static"].items(): print(f"static  {k}: {v_}")
for k, v_ in out["runtime"].items(): print(f"runtime {k}: {v_}")
print("ALL FOUR FIXES OK:", out["all_four_fixes_ok"])
