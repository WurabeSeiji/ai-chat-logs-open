#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス0：干渉保存力学の資格審査（走行前・親生成前に完了必須）。
V0.1 ユニタリ性：乱数状態 1000 step で ‖z‖² 相対ドリフト < 1e-12
V0.2 頂点形恒等式：Hz = conj(z)·(A z²)（read の実装が頂点形と厳密一致）
V0.3 凍結の正当化：連続流 ż=−iH(z)z は ‖z‖ を保存（RK4 で確認）。反線形凍結（M2：S̃ 凍結で ż=−iS̃ z̄）は
     1 step でノルム非保存（双曲的）→ 複素線形凍結 M1=exp(−iΔH) が流れの厳密保存量を保存する唯一の凍結
V0.4 対照テスト：P=Re H を落とすと現行フレーム（v2 補完実験の exp(ΔK)）と厳密一致（シリーズ内完結規約）
V0.5 アンカー定理：hm 親（局所閉塞＋等モジュラー）は H(v)v=−2r²v の厳密固有ベクトル（N=3 は μ=−r²）
V0.6 アンカー上の閉塞不変：hm N=8 を 2000 step 走行して |Σz²|/H が丸め水準に留まる"""
import os, sys, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import edges, adjacency, K_of
from interference_dynamics import hermitian_H, unified_interference_step, current_frame_step, unified_readout, DELTA
import state_provider as sp
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = {}
ok_all = True
def check(name, value, cond, note=""):
    global ok_all
    ok = bool(cond)
    ok_all &= ok
    out[name] = dict(value=value, ok=ok, note=note)
    print(f"{'PASS' if ok else 'FAIL'} {name}: {value:.3e} {note}")

# V0.1 ユニタリ性
drift_max = 0.0
for N in (6, 10, 16):
    A = adjacency(N); M = N*(N-1)//2
    rng = np.random.default_rng(20260831 + N)
    z = rng.standard_normal(M) + 1j*rng.standard_normal(M)
    z *= math.sqrt(M/15.0)/np.linalg.norm(z)
    h0 = float(np.vdot(z, z).real)
    for t in range(1000):
        z = unified_interference_step(z, A)
    drift_max = max(drift_max, abs(float(np.vdot(z, z).real) - h0)/h0)
check("V0.1_unitarity_relative_drift_1000steps", drift_max, drift_max < 1e-12)

# V0.2 頂点形恒等式
dev_max = 0.0
for N in (5, 8, 12):
    A = adjacency(N); M = N*(N-1)//2
    rng = np.random.default_rng(7*N)
    z = rng.standard_normal(M) + 1j*rng.standard_normal(M)
    Hz = hermitian_H(z, A) @ z
    vertex_form = np.conj(z) * (A @ (z*z))
    dev_max = max(dev_max, float(np.linalg.norm(Hz - vertex_form)/np.linalg.norm(Hz)))
check("V0.2_vertex_form_identity", dev_max, dev_max < 1e-14)

# V0.3 流れのノルム保存（RK4）と M2 凍結の非保存
N = 8; A = adjacency(N); M = N*(N-1)//2
rng = np.random.default_rng(3)
z0 = rng.standard_normal(M) + 1j*rng.standard_normal(M)
z0 *= math.sqrt(M/15.0)/np.linalg.norm(z0)
def f(z): return -1j * (hermitian_H(z, A) @ z)
z = z0.copy(); dt = DELTA/50
for t in range(500):
    k1 = f(z); k2 = f(z+dt/2*k1); k3 = f(z+dt/2*k2); k4 = f(z+dt*k3)
    z = z + dt/6*(k1+2*k2+2*k3+k4)
flow_drift = abs(float(np.vdot(z, z).real) - float(np.vdot(z0, z0).real))/float(np.vdot(z0, z0).real)
check("V0.3a_flow_norm_drift_RK4", flow_drift, flow_drift < 1e-10, "(連続流はノルム保存)")
# M2: S̃ を凍結して ż=−i S̃ z̄ を Δ だけ厳密に解く（辺ごと 2×2 実対称生成子 → 双曲的）
s = A @ (z0*z0)
g = 0.0
znew = np.empty(M, complex)
for e in range(M):
    sr, si = s[e].real, s[e].imag
    G2 = np.array([[si, -sr], [-sr, -si]])
    from scipy.linalg import expm
    R2 = expm(DELTA*G2)
    xy = R2 @ np.array([z0[e].real, z0[e].imag])
    znew[e] = xy[0] + 1j*xy[1]
m2_gain = abs(float(np.vdot(znew, znew).real)/float(np.vdot(z0, z0).real) - 1.0)
check("V0.3b_M2_freeze_norm_gain_1step", m2_gain, m2_gain > 1e-8, "(反線形凍結はノルム非保存 → M1 採用の根拠)")

# V0.4 対照テスト：P を落とす → 現行フレームと一致
dev_max = 0.0
for N in (5, 8, 12):
    A = adjacency(N); M = N*(N-1)//2
    rng = np.random.default_rng(11*N)
    z = rng.standard_normal(M) + 1j*rng.standard_normal(M)
    z *= math.sqrt(M/15.0)/np.linalg.norm(z)
    z_cur = current_frame_step(z, A)
    K = K_of(N, z, A)
    w, V = np.linalg.eigh(1j*K)
    z_ref = V @ (np.exp(-1j*DELTA*w) * (V.conj().T @ z))
    dev_max = max(dev_max, float(np.linalg.norm(z_cur - z_ref)/np.linalg.norm(z_ref)))
check("V0.4_contrast_current_frame", dev_max, dev_max < 1e-13)

# V0.5 アンカー定理
worst = 0.0; table = {}
for N in range(3, 17):
    A = adjacency(N)
    v = sp.equimodular(N)
    v *= math.sqrt(len(v)/15.0)/np.linalg.norm(v)
    r2 = float((np.abs(v)**2).mean())
    Hv = hermitian_H(v, A) @ v
    mu = float((np.vdot(v, Hv)/np.vdot(v, v)).real)
    res = float(np.linalg.norm(Hv - mu*v)/np.linalg.norm(v))/r2
    mu_pred = -r2 if N == 3 else -2*r2
    table[N] = dict(mu=mu, mu_pred=mu_pred, res_over_r2=res)
    worst = max(worst, res, abs(mu - mu_pred)/abs(mu_pred))
check("V0.5_anchor_theorem_worst_dev", worst, worst < 1e-12, "(hm: μ=−2r²、N=3 は μ=−r²)")
out["V0.5_table"] = table

# V0.6 アンカー上の閉塞不変（hm N=8、2000 step）
N = 8; A = adjacency(N); E = edges(N)
v = sp.equimodular(N); v *= math.sqrt(len(v)/15.0)/np.linalg.norm(v)
z = v.copy(); cl_max = 0.0
for t in range(2000):
    z = unified_interference_step(z, A)
    ro = unified_readout(z, A, E)
    cl_max = max(cl_max, ro["global_closure"])
check("V0.6_closure_on_anchor_2000steps", cl_max, cl_max < 1e-12)

os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
with open(os.path.join(ROOT, "results", "pass0_qualification.json"), "w") as fjs:
    json.dump(out, fjs, indent=2, ensure_ascii=False, default=float)
if not ok_all:
    raise SystemExit("ABORT: qualification failed")
print("PASS0 OK（資格審査 全項目 PASS）")
