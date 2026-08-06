#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""分布読出し実験: 熱化した二体の海における質量・電荷・スピン・分類の分布

対象: 二体万能非弾性写像（閉形式厳密解）で点火・熱化させた海（f*≈0.469）。
モードは (k, m) = (χ生bin, η巻き数)。各モードの時系列 A_km(t), B_km(t)
（二チャネル振幅）から、本論文で導出した読出しを適用する：

    パリティ分類: 生bin k 偶数≥4 = フェルミオン的 / 奇数 = ボゾン的（確立読出し）
    質量²: 二チャネルGramの detΓ/T²（§8.1 非コヒーレンス）
    スピン的な量: 同Gramの Blochベクトル (X,Y,Z)/T（§8.1の副産物・探索的読出し）
    電荷: 巻き数 m（§9-iii）。符号均衡 balance=(P₊−P₋)/(P₊+P₋) で
          粒子的/反粒子的/灰色を分類

事前明記のオープン課題: 電荷（巻き数）は ±1 に安定せず複数の値に分布する。
これは既知の未解決課題であり、本実験はその分布の実測を与える。

使い方: python3 run_distribution_readouts_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_dist", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0
SEED_AMP = 0.1
J_THERM = 3000
J_WIN = 200


def measure_epoch(a, b, sp, j_win):
    """窓 j_win 衝突の (k,m) モード時系列から全読出しを計算。"""
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    A = np.zeros((j_win, n, ne), complex)
    B = np.zeros((j_win, n, ne), complex)
    for t in range(j_win):
        a, b, _ = ex.collision_step_exact(a, b, sp)
        fa = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
        fb = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
        A[t] = np.fft.fft(fa, axis=1, norm="ortho")
        B[t] = np.fft.fft(fb, axis=1, norm="ortho")
    P = np.mean(np.abs(A) ** 2 + np.abs(B) ** 2, axis=0) / 2
    Gaa = np.mean(np.abs(A) ** 2, axis=0)
    Gbb = np.mean(np.abs(B) ** 2, axis=0)
    Gab = np.mean(A * np.conj(B), axis=0)
    T = (Gaa + Gbb) / 2
    X = np.real(Gab); Y = np.imag(np.conj(Gab)); Z = (Gaa - Gbb) / 2
    det = np.real(Gaa * Gbb - np.abs(Gab) ** 2)
    Tf = np.maximum(T, 1e-300)
    mass2 = det / Tf ** 2
    sz = Z / Tf
    smag = np.sqrt(X ** 2 + Y ** 2 + Z ** 2) / Tf
    ks = np.arange(n); kk = np.where(ks <= n // 2, ks, ks - n)
    ferm_k = (np.abs(kk) % 2 == 0) & (np.abs(kk) >= 4)
    bos_k = (np.abs(kk) % 2 == 1)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    occ = P > (P.max() * 1e-8)
    F = ferm_k[:, None] & occ
    Bc = bos_k[:, None] & occ

    def whist(mask, val, bins):
        w = (P * mask).ravel(); v_ = np.broadcast_to(val, P.shape).ravel()
        h, edges = np.histogram(v_, bins=bins, weights=w)
        return h.tolist(), edges.tolist()

    mh_f, me_ = whist(F, mass2, np.linspace(0, 1.05, 22))
    mh_b, _ = whist(Bc, mass2, np.linspace(0, 1.05, 22))
    sh_f, se_ = whist(F, smag, np.linspace(0, 1.05, 22))
    sh_b, _ = whist(Bc, smag, np.linspace(0, 1.05, 22))
    zh_f, ze_ = whist(F, sz, np.linspace(-1.05, 1.05, 22))
    zh_b, _ = whist(Bc, sz, np.linspace(-1.05, 1.05, 22))
    charge = {int(m_): float(np.sum(P[F & (mm[None, :] == m_)]))
              for m_ in range(-ne // 2, ne // 2)}
    Pk = np.sum(P, axis=1)
    Pk_pos = np.sum(np.where(mm[None, :] > 0, P, 0.0), axis=1)
    Pk_neg = np.sum(np.where(mm[None, :] < 0, P, 0.0), axis=1)
    balance = (Pk_pos - Pk_neg) / np.maximum(Pk_pos + Pk_neg, 1e-300)
    fk = ferm_k & (Pk > Pk.max() * 1e-8)
    bh, be = np.histogram(balance[fk], bins=np.linspace(-1.05, 1.05, 22), weights=Pk[fk])
    tot = float(np.sum(Pk[fk]))
    fr_p = float(np.sum(Pk[fk & (balance > 0.5)]) / tot)
    fr_a = float(np.sum(Pk[fk & (balance < -0.5)]) / tot)
    return a, b, {
        "P_k": Pk.tolist(), "k_signed": kk.tolist(),
        "ferm_mask": ferm_k.tolist(), "bos_mask": bos_k.tolist(),
        "class_power": {"fermionic": float(np.sum(P[F])), "bosonic": float(np.sum(P[Bc]))},
        "mass_hist": {"ferm": mh_f, "bos": mh_b, "edges": me_},
        "smag_hist": {"ferm": sh_f, "bos": sh_b, "edges": se_},
        "sz_hist": {"ferm": zh_f, "bos": zh_b, "edges": ze_},
        "charge_distribution": charge,
        "balance_hist": {"h": bh.tolist(), "edges": be.tolist()},
        "fractions": {"particle": fr_p, "antiparticle": fr_a, "gray": 1 - fr_p - fr_a},
    }


def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    a = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
    a = a + v1.make_bundle(sp, v1.ODD_KS, "A", scale=SEED_AMP * S)
    b = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)

    # 時期1: 成長期（点火の立ち上がり j=400 から窓100）
    for _ in range(400):
        a, b, _ = ex.collision_step_exact(a, b, sp)
    a, b, early = measure_epoch(a, b, sp, 100)
    print("成長期（j=400–500）測定完了")

    # 時期2: 熱平衡期（j=3000 から窓200）
    for _ in range(3000 - 500):
        a, b, _ = ex.collision_step_exact(a, b, sp)
    a, b, late = measure_epoch(a, b, sp, 200)
    print("熱平衡期（j=3000–3200）測定完了")

    for name, ep in (("成長期", early), ("熱平衡期", late)):
        tc = sorted(ep["charge_distribution"].items(), key=lambda kv: -kv[1])[:5]
        print(f"[{name}] F/B={ep['class_power']['fermionic']:.2f}/{ep['class_power']['bosonic']:.2f} "
              f"粒/反/灰={ep['fractions']['particle']:.2f}/{ep['fractions']['antiparticle']:.2f}/{ep['fractions']['gray']:.2f} "
              f"電荷上位={[(m_, round(p_,3)) for m_, p_ in tc]}")

    # 時期3: 帯電した構造（census型: 単一巻きポンプ＋単一巻き種、窓40）
    def single_winding(v):
        shape = (sp.chi_grid_n, sp.eta_grid_n)
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0
        f[sp.chi_grid_n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    a2 = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b2 = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a2 = a2 + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
    _, _, charged = measure_epoch(a2, b2, sp, 40)
    print("帯電構造（census型・窓40）測定完了")
    tc = sorted(charged["charge_distribution"].items(), key=lambda kv: -kv[1])[:6]
    print(f"[帯電構造] 粒/反/灰={charged['fractions']['particle']:.2f}/"
          f"{charged['fractions']['antiparticle']:.2f}/{charged['fractions']['gray']:.2f} "
          f"電荷上位={[(m_, round(p_,3)) for m_, p_ in tc]}")

    out = {"params": {"S": S, "SEED_AMP": SEED_AMP},
           "early": early, "late": late, "charged": charged,
           "open_problem": "電荷（巻き数）は±1に安定せず複数の値に分布する——既知の未解決課題。本実験はその分布の実測。",
           "runtime_sec": time.time() - t0}
    (HERE / "distribution_readouts_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
