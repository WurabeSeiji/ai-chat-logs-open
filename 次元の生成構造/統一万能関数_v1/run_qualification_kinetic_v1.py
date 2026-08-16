#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""資格審査 Q-K1〜Q-K4: 統一万能運動関数 K v1

恒久ルール「新メンバー追加時は正本との対照テストを資格審査に追加してから
使用」に従い、`unified_kinetic_v1.py` を正本
`万能相互作用多体接続_v1/run_kinetic_dispersion_demo_v3.py` と照合する。

判定（実行前固定）
------------------
  Q-K1  正本 KE と本モジュールの分散適用が **ビット同一**（T=400・全步）
  Q-K2  正本の判定 K2'（並進が予言 (x₀−vt) mod N と ±0.05 セル以内）を再現
  Q-K3  正本の判定 K3'（形不変：整数シフトで PR が初期値へ厳密回帰）を再現
  Q-K4  **恒等性**: ω₁=0 で k_dispersion / k_translate_flat は恒等写像
        （相互作用のみの走行＝CR0 と厳密に一致することの保証）
  Q-K5  **並進が Σψ² を厳密保存**（非零 ω₁・**衝突後の状態**で検定）
        Σψ² は **ノルムではなく複素量** Σ(a+ib)²。実部 Σ(a²−b²) と
        虚部 2Σab（交差項）を**別々に**検定する（絶対値へ潰さない）
        ——真の並進なら位相を変えても何も変わらない。当初 Q-K4 までしか
        無く、非零 ω₁ で `k_translate_flat` を検定する項目が欠けていたため、
        生添字による大域位相のバグ（Σψ² が e^{iNω} 回る）を見逃した

使い方: python3 run_qualification_kinetic_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MB = HERE.parent / "万能相互作用多体接続_v1"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


K = _load("unified_kinetic_v1", HERE / "unified_kinetic_v1.py")
g3 = _load("s3kin_q", MB / "run_genesis_v3_register_local_v1.py")
abl, V2 = g3.abl, g3.V2

N_GRAPH, NREG, DELTA = 12, 16, 3e-2
M_EDGE = N_GRAPH * (N_GRAPH - 1) // 2
OMEGA1 = 2 * np.pi / NREG * 0.05
T = 400


class KE_ref(V2):
    """正本 run_kinetic_dispersion_demo_v3.KE の写し（インライン分散）。"""

    def __init__(self, n, C0, wp, omega1, **kw):
        super().__init__(n, C0, wp, **kw)
        self.disp = np.exp(1j * np.arange(self.nreg) * omega1)[None, :]

    def step(self):
        super().step()
        self.C = self.C * self.disp


class KE_new(V2):
    """統一万能運動関数 K を使う版。"""

    def __init__(self, n, C0, wp, omega1, **kw):
        super().__init__(n, C0, wp, **kw)
        self.omega1 = omega1

    def step(self):
        super().step()
        self.C = K.k_dispersion(self.C, self.omega1, axis=1)


def build_C0():
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    seed_edge = g3.zero_closure_state(M_EDGE, np.random.default_rng(98000))
    odd_ks = [k for k in range(NREG) if k % 2 == 1]
    prof = np.zeros(NREG, complex)
    for k in odd_ks:
        prof[k] = 1.0 / np.sqrt(len(odd_ks))
    C0 = np.zeros((M_EDGE, NREG), complex)
    C0[:, 2] = Z0c
    for k in range(NREG):
        if abs(prof[k]) > 0:
            C0[:, k] += DELTA * prof[k] * seed_edge
    return C0, wp0


def cen(C):
    """正本と同一の位置・PR 計器（巻き数2モーメント）。"""
    N = C.shape[1]
    ks = np.arange(N)
    W = np.fft.ifft(C * ((ks % 2 == 1)[None, :]), axis=1) * N
    P = np.sum(np.abs(W) ** 2, axis=0)
    z2 = np.sum(P * np.exp(2j * np.pi * 2 * ks / N)) / P.sum()
    return (float((np.angle(z2) * N / (4 * np.pi)) % (N / 2)),
            float(P.sum() ** 2 / np.sum(P ** 2) / 2))


def main() -> None:
    t0 = time.time()
    C0, wp0 = build_C0()

    ref = KE_ref(N_GRAPH, C0.copy(), wp0, OMEGA1, vertex_on=True)
    new = KE_new(N_GRAPH, C0.copy(), wp0, OMEGA1, vertex_on=True)

    x0, pr0 = cen(ref.C)
    max_diff = 0.0
    rows = []
    k2 = k3 = True
    frac_prs = {}

    for t in range(T):
        ref.step()
        new.step()
        d = float(np.max(np.abs(ref.C - new.C)))
        max_diff = max(max_diff, d)

        if (t + 1) % 50 == 0:
            x, pr = cen(ref.C)
            s = 0.05 * (t + 1)
            pred = (x0 - s) % 8
            dev = min(abs(x - pred), 8 - abs(x - pred))
            frac = round(s % 1, 3)
            rows.append({"t": t + 1, "x": round(x, 4), "pred": round(pred, 4),
                         "dev": round(dev, 6), "PR": round(pr, 4)})
            k2 &= dev <= 0.05
            if frac == 0.0:
                k3 &= abs(pr - pr0) <= 0.01 * pr0
            else:
                if frac in frac_prs:
                    k3 &= abs(pr - frac_prs[frac]) <= 0.01 * frac_prs[frac]
                else:
                    frac_prs[frac] = pr

    # Q-K4: ω₁=0 の恒等性
    rng = np.random.default_rng(7)
    Ztest = (rng.normal(size=(5, 16)) + 1j * rng.normal(size=(5, 16)))
    id_disp = float(np.max(np.abs(K.k_dispersion(Ztest, 0.0, axis=1) - Ztest)))
    psi = (rng.normal(size=512 * 16) + 1j * rng.normal(size=512 * 16))
    id_tr = float(np.max(np.abs(K.k_translate_flat(psi, 0.0, 512, 16) - psi)))

    # Q-K5: 二体正本の衝突後の状態で Σψ² の保存を検定
    _uni = _load("uni_q", HERE / "unified_interaction_v1.py")
    _cr0 = _load("cr0_q", HERE.parent / "電子の反跳実験" / "run_cr0_control_no_theta_v2.py")
    spb = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    nc, nn = int(spb.chi_grid_n), int(spb.eta_grid_n)
    sl, icp, _ = _cr0.calibrate_shift(spb, nc, nn)
    ta, tb = _cr0.make_pair(spb, _cr0.shift_for_deg(-30.0, sl, icp),
                            _cr0.shift_for_deg(+30.0, sl, icp))
    ta, tb, _ = _uni.collision_step_exact(ta, tb, spb)   # 個別閉包を崩す
    S = lambda x: complex(np.sum(x * x))
    s_before = S(ta)
    rng2 = np.random.default_rng(11)
    tc = ta.copy()
    for _ in range(100):
        tc = K.k_translate_flat(tc, float(rng2.normal(0, 2e-3)), nc, nn)
    s_after = S(tc)
    # **Σψ² はノルムではない。** Σ(a+ib)² の複素量であり、
    # 実部 Σ(a²−b²) と虚部 2Σab（交差項）を別々に監視する。
    # 絶対値へ潰してはならず、ノルムで割って基準を緩めることもしない。
    d_re = float(abs(s_after.real - s_before.real))
    d_im = float(abs(s_after.imag - s_before.imag))
    qk5 = bool(d_re <= 1e-12 and d_im <= 1e-12)

    qk1 = bool(max_diff == 0.0)
    qk2 = bool(k2)
    qk3 = bool(k3)
    qk4 = bool(id_disp == 0.0 and id_tr <= 1e-13)

    print("=" * 70)
    print("資格審査 Q-K: 統一万能運動関数 K v1")
    print("=" * 70)
    print(f"  正本 = 万能相互作用多体接続_v1/run_kinetic_dispersion_demo_v3.py")
    print(f"  N={N_GRAPH} Nreg={NREG} ω₁={OMEGA1:.8f} T={T}")
    print(f"  並進速度 = {K.k_cells_per_step(OMEGA1, NREG):.4f} セル/步（−方向）")
    print()
    for r in rows:
        print(f"    t={r['t']:3d}: 位置={r['x']:.4f} 予言={r['pred']:.4f} "
              f"偏差={r['dev']:.6f} PR={r['PR']:.4f}")
    print()
    print(f"  Q-K1 正本とビット同一（最大差 {max_diff:.3e}）: {qk1}")
    print(f"  Q-K2 並進が予言と ±0.05 セル以内: {qk2}")
    print(f"  Q-K3 形不変（整数シフトで PR 回帰）: {qk3}")
    print(f"  Q-K4 ω₁=0 が恒等（分散 {id_disp:.3e} / 並進 {id_tr:.3e}）: {qk4}")
    print(f"  Q-K5 並進が Σψ² を保存（衝突後・100步）: {qk5}")
    print(f"        実部 Σ(a²−b²): {s_before.real:+.6e} → {s_after.real:+.6e}"
          f"  |Δ|={d_re:.3e}")
    print(f"        虚部 2Σab    : {s_before.imag:+.6e} → {s_after.imag:+.6e}"
          f"  |Δ|={d_im:.3e}")
    print()
    ok = all([qk1, qk2, qk3, qk4, qk5])
    print(f"ALL PASS: {ok}  （所要 {time.time()-t0:.1f}s）")

    json.dump({"qualification": "kinetic_v1", "date": time.strftime("%Y-%m-%d %H:%M:%S"),
               "config": {"N": N_GRAPH, "Nreg": NREG, "omega1": OMEGA1, "T": T,
                          "cells_per_step": K.k_cells_per_step(OMEGA1, NREG)},
               "QK1_bitwise": qk1, "QK1_max_diff": max_diff,
               "QK2_translation": qk2, "QK3_shape": qk3,
               "QK4_identity": qk4, "QK4_disp": id_disp, "QK4_translate": id_tr,
               "QK5_sq_conservation": qk5, "QK5_d_re": d_re, "QK5_d_im": d_im,
               "QK5_re_before": s_before.real, "QK5_im_before": s_before.imag,
               "rows": rows, "all_pass": ok},
              open(HERE / "qualification_kinetic_v1_result.json", "w"), indent=1)
    print("保存: qualification_kinetic_v1_result.json")


if __name__ == "__main__":
    main()
