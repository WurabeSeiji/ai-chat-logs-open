#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""インフレーションの軌道に、四元数（八元数）として読める組が生じているかを既存データで調べる（2026-10-19）。

入力（読み取りのみ）：N3_N40_stage123_sweep_20260905/results/hm_N40_den_{40,124}_states_500.npz。
  9/18 の細分系列（振幅込み K、den=200L）の states.npz は削除済みのため使えない。
生成子：run_N3_N40_stage123_v1.py の adjacency と one_step 内の K（位相のみ K_ef = A_ef sin(φ_f − φ_e)）を写したもの。

測る量（事前に固定。窓は各段：床・指数成長・遷移・終端）：
  sc_k      : 窓内の複素状態行列 Z（T×M）の特異値（第 1 で規格化）。複素の張りの階層
  sr_k      : 実の張り [Re Z; Im Z]（2T×M）の特異値。倍加の段（2→4→8）に乗るかを見る
  ov_J2     : 第 2 複素モード u₂ と、第 1 モードの四元数の相手 J₂u₁ = conj(J_K u₁) の重なり |⟨J₂u₁, u₂⟩|
              （J_K は窓中央の K の極分解の直交部分。気づきメモ §2 の J₂ = C∘J_K）
  ov_conj   : u₂ と conj(u₁) の重なり
  leak_J2   : 窓の複素の張り（相対特異値 > 1e-6）が J₂ で閉じない割合
注意：J₂ が i と反可換であることは代数的恒等式で、ここでは測らない。測るのはデータの張りが J₂ で閉じるか、
  育つ方向が四元数の相手かどうか。
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, '..', '電子の反跳実験', 'seed除去による準安定相追試', 'chatgpt追試',
                 '干渉保存力学_資格審査とシード無し系列_20260831', 'N3_N40_stage123_sweep_20260905', 'results')
N = 40
WINDOWS = {40: [('床', 10, 40), ('指数成長', 150, 180), ('遷移前半', 300, 330), ('遷移後半', 330, 360),
                ('終端手前', 360, 390), ('終端', 470, 500)],
           124: [('床', 5, 35), ('指数成長', 40, 70), ('遷移前半', 80, 110), ('遷移後半', 110, 140),
                 ('終端手前', 140, 170), ('終端', 470, 500)]}


def adjacency(N):
    # run_N3_N40_stage123_v1.py の adjacency と同じ
    ea, eb = np.triu_indices(N, k=1); M = len(ea); A = np.zeros((M, M), dtype=np.float64)
    for e in range(M):
        share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e]); share[e] = False; A[e, share] = 1.0
    return A


A = adjacency(N)


def K_of(z):
    # one_step の H=H_of(exp(1j*angle(z)),A) の虚部 = K（位相のみ）
    u = np.exp(1j * np.angle(z)); H = A * (np.conj(u)[:, None] * u[None, :]); np.fill_diagonal(H, 0.0)
    return np.imag(H)


def polar_J(K):
    U, s, Vt = np.linalg.svd(K); r = s > 1e-9 * s[0]
    return U[:, r] @ Vt[r, :], int(r.sum())


rows = []
with np.errstate(all='ignore'):   # numpy 2.0 + Accelerate の matmul 偽警告を抑える（値は有限であることを下で確認）
    for den, wins in WINDOWS.items():
        S = np.load(os.path.join(R, f'hm_N{N}_den_{den}_states_500.npz'))['Z'].astype(np.complex128)
        for name, a, b in wins:
            Z = S[a:b]
            U, s, Vh = np.linalg.svd(Z, full_matrices=False); sc = s / s[0]
            sr = np.linalg.svd(np.vstack([Z.real, Z.imag]), compute_uv=False); sr = sr / sr[0]
            d = int(np.sum(sc > 1e-6)); B = Vh[:d].conj().T
            J, rankK = polar_J(K_of(S[(a + b) // 2]))
            J2B = np.conj(J @ B)
            leak = float(np.linalg.norm(J2B - B @ (B.conj().T @ J2B)) / np.linalg.norm(J2B))
            u1, u2 = Vh[0].conj(), Vh[1].conj()
            j2u1 = np.conj(J @ u1); j2u1 /= np.linalg.norm(j2u1)
            ov_J2 = float(abs(np.vdot(j2u1, u2)))
            ov_conj = float(abs(np.vdot(np.conj(u1), u2)))
            assert np.all(np.isfinite([leak, ov_J2, ov_conj]))
            rows.append(dict(den=den, phase=name, window=[a, b], sc=sc[:6].tolist(), sr=sr[:8].tolist(),
                             d_complex_1e6=d, leak_J2=leak, ov_J2=ov_J2, ov_conj=ov_conj, rankK=rankK))

os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
json.dump(rows, open(os.path.join(HERE, 'results', 'quaternion_pairing_N40.json'), 'w'), ensure_ascii=False, indent=1)
with open(os.path.join(HERE, 'results', 'quaternion_pairing_N40.md'), 'w', encoding='utf-8') as f:
    f.write('| den | 段 | 窓 | 複素特異値 s₂/s₁, s₃/s₁, s₄/s₁ | 実特異値（第 3〜6） | 複素の張り (>1e-6) | J₂ で閉じない割合 | ⟨J₂u₁, u₂⟩ | ⟨ū₁, u₂⟩ |\n')
    f.write('|---|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        f.write(f"| {r['den']} | {r['phase']} | {r['window'][0]}–{r['window'][1]} | "
                f"{r['sc'][1]:.1e}, {r['sc'][2]:.1e}, {r['sc'][3]:.1e} | "
                f"{', '.join(f'{x:.1e}' for x in r['sr'][2:6])} | {r['d_complex_1e6']} | {r['leak_J2']:.3f} | "
                f"{r['ov_J2']:.3f} | {r['ov_conj']:.3f} |\n")
print(open(os.path.join(HERE, 'results', 'quaternion_pairing_N40.md'), encoding='utf-8').read())
