#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40, L=40000, Δτ=2π/N 実験
干渉保存力学で、分母が N に等しい場合のインフレーション観察"""
import os, math, csv, json, platform
import numpy as np
import matplotlib.pyplot as plt

IN=os.path.join(os.path.dirname(__file__), '..', 'hm_mp_free_N3_N40_20260901', 'data')
OUT=os.path.dirname(__file__)
N = 40
STEPS = 40000
DEN = N  # Δτ = 2π/N
KEY_STEPS = sorted(set([0, 10, 50, 100, 500, 1000, 2000, 5000, 10000, 20000, 40000]))

assert np.dtype(np.float64).itemsize==8 and np.dtype(np.complex128).itemsize==16

def edges(N):
    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)

def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    return A

def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)

def one_step(z,A,den):
    H=H_of(z,A); w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)

def plane(v):
    p=v.real.astype(np.float64,copy=True); p/=np.linalg.norm(p); q=v.imag.astype(np.float64,copy=True); q-=np.dot(q,p)*p; q/=np.linalg.norm(q); return p,q

def metrics(z,p,q):
    h=np.vdot(z,z).real; zp=z-p*np.dot(p,z)-q*np.dot(q,z); hp=np.vdot(zp,zp).real; return float(hp/h),float(h),float(abs(z@z)/h)

print(f"N={N}, den={DEN}, STEPS={STEPS}, KEY_STEPS={len(KEY_STEPS)}")

# 入力データ読み込み
npz_file = os.path.join(IN, f'hm_N{N}', 'parent_v.npz')
if not os.path.exists(npz_file):
    print(f'ERROR: {npz_file} not found')
    exit(1)

data = np.load(npz_file)
z0 = np.array(data['v'], dtype=np.complex128, copy=True)
A = adjacency(N)
p, q = plane(z0)

print(f"z0 norm: {np.linalg.norm(z0):.10f}")

# 力学実行
z = z0.copy()
vals = []
key_data = []

for t in range(STEPS + 1):
    h_perp, h_tot, closure = metrics(z, p, q)
    vals.append((t, h_perp, h_tot, closure))

    if t in KEY_STEPS:
        key_data.append((t, h_perp, h_tot, closure))

    if t < STEPS:
        z = one_step(z, A, DEN)

    if (t + 1) % 5000 == 0:
        print(f'  step {t+1:5d}: Hperp/H={h_perp:.3e}', flush=True)

# 判定
vals_arr = np.array(vals)
onset_idx = np.flatnonzero(vals_arr[:, 1] > 0.05)
onset = int(onset_idx[0]) if onset_idx.size else -1

print(f"\n=== 結果 ===")
print(f"Onset (Hperp/H > 0.05): {onset}")
print(f"Initial Hperp/H:        {vals_arr[0, 1]:.3e}")
print(f"Final Hperp/H:          {vals_arr[-1, 1]:.3e}")
print(f"Max Hperp/H:            {vals_arr[:, 1].max():.3e}")

# CSV 出力
with open(os.path.join(OUT, 'timeseries_N40_L40000.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['step', 'Hperp_frac', 'H_total', 'global_closure'])
    w.writerows(vals)

with open(os.path.join(OUT, 'key_steps_N40_L40000.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['step', 'Hperp_frac', 'H_total', 'global_closure'])
    w.writerows(key_data)

# Summary
summary = {
    'N': N,
    'denominator': DEN,
    'steps': STEPS,
    'onset_Hperp_gt_0.05': onset,
    'initial_Hperp_frac': float(vals_arr[0, 1]),
    'final_Hperp_frac': float(vals_arr[-1, 1]),
    'max_Hperp_frac': float(vals_arr[:, 1].max()),
    'dtype_state': 'complex128',
    'dtype_real': 'float64',
    'numpy': np.__version__,
    'python': platform.python_version()
}

with open(os.path.join(OUT, 'summary_N40_L40000.json'), 'w') as f:
    json.dump(summary, f, indent=2)

# 図：対数グラフ
fig, ax = plt.subplots(figsize=(12, 6))
ax.semilogy(vals_arr[:, 0], vals_arr[:, 1], 'b-', linewidth=1.5, label=f'N={N}, Δτ=2π/{DEN}')
ax.axhline(0.05, color='r', linestyle='--', linewidth=1, label='threshold 0.05')
if onset >= 0:
    ax.axvline(onset, color='g', linestyle=':', linewidth=1, alpha=0.7, label=f'onset={onset}')
ax.set_xlabel('step'); ax.set_ylabel('Hperp/H'); ax.set_title(f'N={N}, L={STEPS}, Δτ=2π/N')
ax.grid(alpha=0.3); ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_N40_L40000_dynamics.png'), dpi=180)
plt.close(fig)

# 図：初期段階の詳細（最初 1000 ステップ）
if STEPS >= 1000:
    fig, ax = plt.subplots(figsize=(12, 6))
    mask = vals_arr[:, 0] <= 1000
    ax.semilogy(vals_arr[mask, 0], vals_arr[mask, 1], 'b-', linewidth=1.5)
    ax.set_xlabel('step'); ax.set_ylabel('Hperp/H'); ax.set_title(f'N={N}, L=1000 (初期段階詳細)')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'fig_N40_L40000_initial_1000.png'), dpi=180)
    plt.close(fig)

print(f"\n✓ 実行完了")
print(f"  timeseries: {os.path.join(OUT, 'timeseries_N40_L40000.csv')}")
print(f"  summary: {os.path.join(OUT, 'summary_N40_L40000.json')}")
