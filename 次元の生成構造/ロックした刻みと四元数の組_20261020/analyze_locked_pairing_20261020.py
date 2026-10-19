#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ロックした刻み（q=8, 4）の走行に、四元数として読める組が生じているかを判定する（2026-10-20）。

判定は ../インフレーションの準振動と一周期ロック_20261019/analyze_quaternion_pairing_20261019.py と同じ量：
  sc_k（複素の張りの特異値）、sr_k（実の張り）、ov_J2 = |<J2 u1, u2>|（J2 = C∘J_K）、ov_conj = |<conj u1, u2>|、
  leak_J2（窓の複素の張りが J2 で閉じない割合）。
窓（事前に固定）：長さ 4q ステップ（床の回転ちょうど 4 周）。段は H⊥/H で決める：
  床＝最初の窓、指数成長＝log10(H⊥/H) が −15 を最初に超える点から、遷移＝1e-3 を超える点から、
  終端＝最後の窓。比較対照として den=40（ロックなし）の同じ判定値を前回の結果から並べる。
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, 'run_N40_locked_step_20261020.py')).read().split("if __name__")[0])

N = 40
A = adjacency(N)


def K_of(z):
    return np.imag(H_of(np.exp(1j * np.angle(z)), A))


def polar_J(K):
    U, s, Vt = np.linalg.svd(K); r = s > 1e-9 * s[0]
    return U[:, r] @ Vt[r, :]


def judge(S, a, b):
    Z = S[a:b]
    U, s, Vh = np.linalg.svd(Z, full_matrices=False); sc = s / s[0]
    sr = np.linalg.svd(np.vstack([Z.real, Z.imag]), compute_uv=False); sr = sr / sr[0]
    d = int(np.sum(sc > 1e-6)); B = Vh[:d].conj().T
    J = polar_J(K_of(S[(a + b) // 2]))
    J2B = np.conj(J @ B)
    leak = float(np.linalg.norm(J2B - B @ (B.conj().T @ J2B)) / np.linalg.norm(J2B))
    u1, u2 = Vh[0].conj(), Vh[1].conj()
    j2u1 = np.conj(J @ u1); j2u1 /= np.linalg.norm(j2u1)
    return dict(window=[int(a), int(b)], sc=sc[:6].tolist(), sr=sr[:8].tolist(), d_complex_1e6=d, leak_J2=leak,
                ov_J2=float(abs(np.vdot(j2u1, u2))), ov_conj=float(abs(np.vdot(np.conj(u1), u2))))


rows = []
with np.errstate(all='ignore'):
    for qord in (8, 4):
        f = os.path.join(HERE, 'results', f'locked_N40_q{qord}_states.npz')
        if not os.path.exists(f):
            continue
        S = np.load(f)['Z'].astype(np.complex128)
        p, q = plane(S[0]); h = np.array([metrics(z, p, q)[0] for z in S]); lg = np.log10(np.maximum(h, 1e-40))
        W = 4 * qord; T = len(S)
        marks = [('床', 1)]
        if np.any(lg > -15): marks.append(('指数成長', int(np.argmax(lg > -15))))
        if np.any(lg > -3): marks.append(('遷移', int(np.argmax(lg > -3))))
        marks.append(('終端', T - W))
        for name, a in marks:
            a = min(max(a, 0), T - W)
            r = judge(S, a, a + W); r.update(q=qord, phase=name, log10_Hperp_at_start=float(lg[a]))
            rows.append(r)

json.dump(rows, open(os.path.join(HERE, 'results', 'locked_pairing_N40.json'), 'w'), ensure_ascii=False, indent=1)
with open(os.path.join(HERE, 'results', 'locked_pairing_N40.md'), 'w', encoding='utf-8') as f:
    f.write('| q | 段 | 窓 | log₁₀(H⊥/H) | 複素特異値 s₂/s₁, s₃/s₁, s₄/s₁ | 実特異値（第 3〜6） | 複素の張り | J₂ で閉じない割合 | ⟨J₂u₁, u₂⟩ | ⟨ū₁, u₂⟩ |\n')
    f.write('|---|---|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        f.write(f"| {r['q']} | {r['phase']} | {r['window'][0]}–{r['window'][1]} | {r['log10_Hperp_at_start']:.1f} | "
                f"{r['sc'][1]:.1e}, {r['sc'][2]:.1e}, {r['sc'][3]:.1e} | "
                f"{', '.join(f'{x:.1e}' for x in r['sr'][2:6])} | {r['d_complex_1e6']} | {r['leak_J2']:.3f} | "
                f"{r['ov_J2']:.3f} | {r['ov_conj']:.3f} |\n")
print(open(os.path.join(HERE, 'results', 'locked_pairing_N40.md'), encoding='utf-8').read())
