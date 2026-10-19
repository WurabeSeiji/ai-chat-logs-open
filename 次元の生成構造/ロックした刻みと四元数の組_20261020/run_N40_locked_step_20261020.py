#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40 の床を、系自身の回転にロックした刻みで回す（2026-10-20）。

元プログラム（忠実コピー）：次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/
  干渉保存力学_資格審査とシード無し系列_20260831/N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py
  edges / adjacency / H_of / one_step / plane / metrics を文字どおり写した。
唯一の変更：one_step の刻み Δτ = 2π/den を、one_step_dtau(z, A, dtau) で直接与えられるようにした（本文は同一）。
  ロックした刻み：Δτ = 2π/(qσ)。σ は床（親 Z0）の iK 固有値から測る。床は一歩ごとにちょうど 2π/q 回り U^q = I。
対照：元の one_step(z, A, 40) で 500 step を回し、保存済み hm_N40_den_40_states_500.npz とビット一致を確認する。
総時間：元の den=40 走行と同じ τ = 500·2π/40 = 25π を覆うステップ数にする。
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '電子の反跳実験', 'seed除去による準安定相追試', 'chatgpt追試',
                   '干渉保存力学_資格審査とシード無し系列_20260831', 'N3_N40_stage123_sweep_20260905')
OUT = os.path.join(HERE, 'results')


# ---- 以下、元プログラムからの文字どおりのコピー ----
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
    # 段3の最小変更（唯一の力学変更点）: 位相のみ生成子 Ĥ の虚部だけを取る H=i·K（K=sin(Δθ) 実反対称）。
    # exp(-iΔτ·iK)=exp(Δτ·K) の実直交回転となり、Z^T Z（零閉塞）と ‖Z‖ を厳密保存する。
    H=H_of(np.exp(1j*np.angle(z)),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)
    w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)
def plane(v):
    p=v.real.astype(np.float64,copy=True); p/=np.linalg.norm(p); q=v.imag.astype(np.float64,copy=True); q-=np.dot(q,p)*p; q/=np.linalg.norm(q); return p,q
def metrics(z,p,q):
    h=np.vdot(z,z).real; zp=z-p*np.dot(p,z)-q*np.dot(q,z); hp=np.vdot(zp,zp).real; return float(hp/h),float(h),float(abs(z@z)/h)
# ---- コピーここまで ----


def one_step_dtau(z,A,dtau):
    # one_step と同一。刻みを den ではなく dtau で直接受け取る（唯一の変更）
    H=H_of(np.exp(1j*np.angle(z)),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)
    w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(dtau)*w)
    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)


def run(z0, A, steps, stepper):
    p, q = plane(z0); z = z0.copy()
    states = np.empty((steps + 1, z.size), np.complex128); vals = np.empty(steps + 1)
    for t in range(steps + 1):
        states[t] = z; vals[t] = metrics(z, p, q)[0]
        if t < steps:
            z = stepper(z)
    return states, vals


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    N = 40
    z0 = np.array(np.load(os.path.join(SRC, 'parents', f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],
                  dtype=np.complex128, copy=True)
    A = adjacency(N)
    meta = {}
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'

    if mode in ('control', 'all'):
        states, vals = run(z0, A, 500, lambda z: one_step(z, A, 40))
        ref = np.load(os.path.join(SRC, 'results', 'hm_N40_den_40_states_500.npz'))['Z']
        meta['control_den40_bit_identical'] = bool(np.array_equal(states, ref))
        meta['control_den40_max_abs_diff'] = float(np.max(np.abs(states - ref)))
        print('control bit-identical:', meta['control_den40_bit_identical'], meta['control_den40_max_abs_diff'], flush=True)

    # σ：床 z0 の iK 固有値（レイリー商）
    K = np.imag(H_of(np.exp(1j * np.angle(z0)), A))
    mu = (np.vdot(z0, (1j * K) @ z0) / np.vdot(z0, z0)).real
    resid = float(np.linalg.norm((1j * K) @ z0 - mu * z0) / np.linalg.norm(z0))
    sigma = -mu
    meta.update(sigma=float(sigma), mu=float(mu), eigen_residual=resid)
    print(f'sigma={sigma:.9f} residual={resid:.2e}', flush=True)

    TAU = 500 * 2 * math.pi / 40
    for qord in (8, 4):
        if mode not in (f'q{qord}', 'all'):
            continue
        dtau = 2 * math.pi / (qord * sigma)
        steps = int(round(TAU / dtau))
        states, vals = run(z0, A, steps, lambda z: one_step_dtau(z, A, dtau))
        np.savez_compressed(os.path.join(OUT, f'locked_N40_q{qord}_states.npz'), Z=states, N=np.int64(N),
                            q=np.int64(qord), dtau=np.float64(dtau), sigma=np.float64(sigma), steps=np.int64(steps))
        meta[f'q{qord}'] = dict(dtau=dtau, steps=steps, Hperp_final=float(vals[-1]), Hperp_max=float(vals.max()),
                                onset_gt_1em3=int(np.argmax(vals > 1e-3)) if np.any(vals > 1e-3) else -1)
        print(f'q={qord} dtau={dtau:.6f} steps={steps} final Hperp/H={vals[-1]:.3e} max={vals.max():.3e}', flush=True)

    json.dump(meta, open(os.path.join(OUT, f'run_meta_{mode}.json'), 'w'), indent=1)
