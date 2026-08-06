#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 P2: モード分解番地census——低N領域（第一周期候補）の占有表

仮説（事前記録）: 各SVD平面の回転数 ρ = f_平面/f_時計 は U^n=I の有理番地 m/n に
量子化され、占有番地の集合が N とともに変化する（種の出現閾値）。
走査v1のジグザグ点（ピークN=9,15／谷N=10,13）で番地集合が変わるか。

規約（事前固定）:
  - N=4..16 全点。T=4000、窓[2000,4000]、400サンプル。
  - 各SVD平面 j（上位5平面）: 平面複素読出し c_j = Sp·conj(d_{2j}) の
    支配周波数 f_j（FFTピーク・零詰め4倍で分解能補強）と時計 f_clock の比 ρ_j。
  - 番地写像: ρ_j の連分数近似（分母≤16、特別族{5,25,31,72,124,144}）。
    プラトー資格: FFTピーク幅による誤差 δρ を併記し、|ρ−m/n|<δρ のみ番地候補。
  - 各平面の質量度（対Gram非コヒーレンス）とパワー占有率も記録
    → 質量階層の痕跡（R3: 番地分母と質量度の相関）。

読み（記述的）:
  R1 占有番地集合のN依存（出現閾値）
  R2 ジグザグ点での番地集合の変化
  R3 番地分母×質量度の相関
  R4 回転数がそもそも有理プラトーに乗るか（乗らなければ本計画の反証として記録）

使い方: python3 run_pre_mode_census_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
spec1 = importlib.util.spec_from_file_location(
    "pre1_mc", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl

T_END = 4000
WIN = (2000, 4000)
SAMPLE_EVERY = 5
N_LIST = tuple(range(4, 17))
N_PLANES = 5
SPECIAL_DENOMS = (5, 25, 31, 72, 124, 144)
ZERO_PAD = 4


def address_fit(rho: float, drho: float) -> dict:
    fr = Fraction(rho).limit_denominator(16)
    dev_small = rho - float(fr)
    best_sp = None
    for nd in SPECIAL_DENOMS:
        m = int(round(rho * nd))
        if m < 1:
            continue
        d = rho - m / nd
        if best_sp is None or abs(d) < abs(best_sp[2]):
            best_sp = (m, nd, d)
    return {"small": [fr.numerator, fr.denominator, float(dev_small),
                       bool(abs(dev_small) < drho)],
            "special": ([best_sp[0], best_sp[1], float(best_sp[2]),
                          bool(abs(best_sp[2]) < drho)] if best_sp else None)}


def scan_one(n: int) -> dict:
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    M = sys_lr.m
    samples = []
    for t in range(T_END):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if WIN[0] <= t < WIN[1] and (t % SAMPLE_EVERY == 0):
            samples.append(Z.copy())
    S = np.array(samples)
    ns = S.shape[0]
    Sp = S - np.outer(S @ p, p) - np.outer(S @ q, q)

    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)

    cp = S @ p
    cq = S @ q
    phi = np.unwrap(np.angle(cp + 1j * cq))
    om_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])
    f_clock = om_clock / (2 * np.pi)

    def freq_of(c):
        cc = c - c.mean()
        F = np.abs(np.fft.fft(cc, n=ZERO_PAD * ns))
        pk = int(np.argmax(F))
        L = ZERO_PAD * ns
        f = pk / L if pk <= L // 2 else (pk - L) / L
        return f

    tot_pow = float(np.sum(sv[:2 * N_PLANES] ** 2))
    planes = []
    for j in range(N_PLANES):
        if 2 * j + 1 >= Vt.shape[0]:
            break
        d1, d2 = Vt[2 * j], Vt[2 * j + 1]
        d1c = d1[:M] + 1j * d1[M:]; d1c /= np.linalg.norm(d1c)
        d2c = d2[:M] + 1j * d2[M:]; d2c /= np.linalg.norm(d2c)
        c1 = Sp @ np.conj(d1c)
        c2 = Sp @ np.conj(d2c)
        f1 = freq_of(c1)
        rho = abs(f1) / abs(f_clock) if f_clock != 0 else 0.0
        drho = (1.0 / (ZERO_PAD * ns)) / abs(f_clock)
        # 質量度
        n1 = round(2 * np.pi * f1 / om_clock) if om_clock != 0 else 0
        r1 = c1 * np.exp(-1j * n1 * phi)
        f2 = freq_of(c2)
        n2 = round(2 * np.pi * f2 / om_clock) if om_clock != 0 else 0
        r2 = c2 * np.exp(-1j * n2 * phi)
        G11 = float(np.mean(np.abs(r1) ** 2)); G22 = float(np.mean(np.abs(r2) ** 2))
        G12 = complex(np.mean(r1 * np.conj(r2)))
        T = 0.5 * (G11 + G22)
        mdeg = float((G11 * G22 - abs(G12) ** 2) / T ** 2) if T > 0 else 0.0
        planes.append({"plane": j + 1,
                        "power_share": float((sv[2 * j] ** 2 + sv[2 * j + 1] ** 2) / tot_pow),
                        "rho": float(rho), "drho": float(drho),
                        "address": address_fit(rho, drho),
                        "mass_deg": mdeg})
    return {"N": n, "f_clock_per_sample": f_clock,
            "clock_over_pi72_step": float(om_clock / SAMPLE_EVERY / (np.pi / 72)),
            "planes": planes}


def main() -> None:
    t0 = time.time()
    out = {"T_END": T_END, "WIN": list(WIN), "SAMPLE_EVERY": SAMPLE_EVERY,
           "N_LIST": list(N_LIST), "N_PLANES": N_PLANES, "ZERO_PAD": ZERO_PAD,
           "SPECIAL_DENOMS": list(SPECIAL_DENOMS), "scan": []}
    for n in N_LIST:
        r = scan_one(n)
        out["scan"].append(r)
        print(f"N={n:2d} clock/step比={r['clock_over_pi72_step']:.4f}")
        for pl in r["planes"]:
            a = pl["address"]
            sm = a["small"]
            sp = a["special"]
            tag_s = f"{sm[0]}/{sm[1]}{'✓' if sm[3] else '×'}"
            tag_p = (f"{sp[0]}/{sp[1]}{'✓' if sp[3] else '×'}" if sp else "-")
            print(f"    平面{pl['plane']}: 占有={pl['power_share']:.3f} "
                  f"ρ={pl['rho']:.5f}±{pl['drho']:.5f} 番地小={tag_s} 特={tag_p} "
                  f"質量度={pl['mass_deg']:.2e}")
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_mode_census_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_mode_census_result_v1.json")


if __name__ == "__main__":
    main()
