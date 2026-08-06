#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v2: 長窓census——番地の一意判定と質量=離調二乗則の精密化

v1 の限界（δρ=0.018 で番地一意同定不能）を長窓で突破する。

規約（事前固定）:
  - N=5..16。T=42000、窓[2000,42000]（40000ステップ・8000サンプル）。
  - SVD基底は全窓で1回だけ推定し固定（部分窓のドリフトを力学に帰属させるため）。
  - 各平面: 全窓FFT（零詰め4倍）で ρ=f/f_clock と分解能 δρ を測る。
    第2ピークとのパワー比も記録（線分裂=非定常の指標）。
  - 部分窓検定: 窓を4分割し、固定基底で ρ を再測。最大ドリフトが 2δρ_sub を
    超える平面は「非定常」と記録（番地判定の資格なし）。
  - 番地判定 H_plateau: 定常な平面について |ρ−m/n| < δρ となる最小分母の
    m/n（分母≤24）と特別族{5,25,31,72,124,144,248}を報告。
  - 質量則検定 H_C: 定常平面全点で 質量度 = C·|ρ−1|^p をフィット。
    p の95%CIに 2 が入るか、C の95%CIに sin²(23π/124)=0.302822 が入るかを判定。

使い方: python3 run_census_longwindow_v2.py
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
    "pre1_lw", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl

T_END = 42000
WIN = (2000, 42000)
SAMPLE_EVERY = 5
N_LIST = tuple(range(5, 17))
N_PLANES = 5
ZERO_PAD = 4
SPECIAL_DENOMS = (5, 25, 31, 72, 124, 144, 248)
A23 = float(np.sin(23 * np.pi / 124) ** 2)


def peak_freq(c, pad):
    cc = c - c.mean()
    L = pad * len(cc)
    F = np.abs(np.fft.fft(cc, n=L))
    pk = int(np.argmax(F))
    f1 = pk / L if pk <= L // 2 else (pk - L) / L
    p1 = F[pk]
    # 第2ピーク（第1の±3ビン近傍を除外）
    F2 = F.copy()
    lo = max(0, pk - 3 * pad); hi = min(L, pk + 3 * pad + 1)
    F2[lo:hi] = 0.0
    pk2 = int(np.argmax(F2))
    f2 = pk2 / L if pk2 <= L // 2 else (pk2 - L) / L
    return f1, float(F2[pk2] / max(p1, 1e-300)), f2


def address_fit(rho, drho):
    fr = Fraction(rho).limit_denominator(24)
    dev = rho - float(fr)
    hit_small = bool(abs(dev) < drho)
    best_sp = None
    for nd in SPECIAL_DENOMS:
        m = int(round(rho * nd))
        if m < 1:
            continue
        d = rho - m / nd
        if best_sp is None or abs(d) < abs(best_sp[2]):
            best_sp = (m, nd, d)
    return {"small": [fr.numerator, fr.denominator, float(dev), hit_small],
            "special": [best_sp[0], best_sp[1], float(best_sp[2]),
                         bool(abs(best_sp[2]) < drho)]}


def scan_one(n: int) -> dict:
    t0 = time.time()
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

    tot_pow = float(np.sum(sv[:2 * N_PLANES] ** 2))
    planes = []
    q4 = ns // 4
    for j in range(N_PLANES):
        d1, d2 = Vt[2 * j], Vt[2 * j + 1]
        d1c = d1[:M] + 1j * d1[M:]; d1c /= np.linalg.norm(d1c)
        d2c = d2[:M] + 1j * d2[M:]; d2c /= np.linalg.norm(d2c)
        c1 = Sp @ np.conj(d1c)
        c2 = Sp @ np.conj(d2c)

        f1, split_ratio, f_2nd = peak_freq(c1, ZERO_PAD)
        rho = abs(f1) / abs(f_clock)
        drho = (1.0 / (ZERO_PAD * ns)) / abs(f_clock)

        # 部分窓ドリフト（固定基底）
        rhos = []
        for k in range(4):
            seg = c1[k * q4:(k + 1) * q4]
            fs, _, _ = peak_freq(seg, ZERO_PAD)
            # 時計も同区間で
            ph = phi[k * q4:(k + 1) * q4]
            omc = float(np.polyfit(np.arange(len(ph)), ph, 1)[0])
            rhos.append(abs(fs) / abs(omc / (2 * np.pi)))
        drift = float(max(rhos) - min(rhos))
        drho_sub = (1.0 / (ZERO_PAD * q4)) / abs(f_clock)
        stationary = bool(drift < 2 * drho_sub)

        # 質量度（全窓・復調）
        n1 = round(2 * np.pi * f1 / om_clock) if om_clock != 0 else 0
        r1 = c1 * np.exp(-1j * n1 * phi)
        f2p, _, _ = peak_freq(c2, ZERO_PAD)
        n2 = round(2 * np.pi * f2p / om_clock) if om_clock != 0 else 0
        r2 = c2 * np.exp(-1j * n2 * phi)
        G11 = float(np.mean(np.abs(r1) ** 2)); G22 = float(np.mean(np.abs(r2) ** 2))
        G12 = complex(np.mean(r1 * np.conj(r2)))
        T = 0.5 * (G11 + G22)
        mdeg = float((G11 * G22 - abs(G12) ** 2) / T ** 2) if T > 0 else 0.0

        planes.append({"plane": j + 1,
                        "power_share": float((sv[2 * j] ** 2 + sv[2 * j + 1] ** 2) / tot_pow),
                        "rho": float(rho), "drho": float(drho),
                        "rho_subwins": [float(x) for x in rhos],
                        "drift": drift, "drho_sub": float(drho_sub),
                        "stationary": stationary,
                        "split_ratio": split_ratio,
                        "address": address_fit(rho, drho),
                        "mass_deg": mdeg})
    return {"N": n, "clock_over_pi72_step": float(om_clock / SAMPLE_EVERY / (np.pi / 72)),
            "planes": planes, "runtime_sec": time.time() - t0}


def main() -> None:
    t0 = time.time()
    out = {"T_END": T_END, "WIN": list(WIN), "SAMPLE_EVERY": SAMPLE_EVERY,
           "N_LIST": list(N_LIST), "ZERO_PAD": ZERO_PAD, "A23": A23, "scan": []}
    for n in N_LIST:
        r = scan_one(n)
        out["scan"].append(r)
        print(f"N={n:2d} clock比={r['clock_over_pi72_step']:.5f} [{r['runtime_sec']:.0f}s]")
        for pl in r["planes"]:
            a = pl["address"]
            sm, sp = a["small"], a["special"]
            st = "定常" if pl["stationary"] else f"漂({pl['drift']:.4f})"
            print(f"    平面{pl['plane']}: 占有={pl['power_share']:.3f} "
                  f"ρ={pl['rho']:.6f}±{pl['drho']:.6f} {st} 分裂比={pl['split_ratio']:.2f} "
                  f"小={sm[0]}/{sm[1]}{'✓' if sm[3] else '×'} "
                  f"特={sp[0]}/{sp[1]}{'✓' if sp[3] else '×'} "
                  f"質量度={pl['mass_deg']:.2e}")
        out["runtime_sec"] = time.time() - t0
        (HERE / "census_longwindow_result_v2.json").write_text(
            json.dumps(out, indent=1, ensure_ascii=False))

    # H_C: 定常平面のみで冪則フィット
    xs, ys = [], []
    for r in out["scan"]:
        for pl in r["planes"]:
            dr = abs(pl["rho"] - 1.0)
            if pl["stationary"] and dr > 1e-4 and pl["mass_deg"] > 0:
                xs.append(dr); ys.append(pl["mass_deg"])
    if len(xs) >= 8:
        lx, ly = np.log(np.array(xs)), np.log(np.array(ys))
        A = np.vstack([lx, np.ones_like(lx)]).T
        coef, res, _, _ = np.linalg.lstsq(A, ly, rcond=None)
        pred = A @ coef
        resid = ly - pred
        s2 = float(np.sum(resid ** 2) / (len(lx) - 2))
        cov = s2 * np.linalg.inv(A.T @ A)
        p_se, lc_se = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))
        C = float(np.exp(coef[1]))
        C_lo, C_hi = float(np.exp(coef[1] - 2 * lc_se)), float(np.exp(coef[1] + 2 * lc_se))
        r2 = float(1 - np.sum(resid ** 2) / np.sum((ly - ly.mean()) ** 2))
        hC = bool(C_lo <= A23 <= C_hi)
        hp = bool(coef[0] - 2 * p_se <= 2.0 <= coef[0] + 2 * p_se)
        print(f"\nH_C 冪則（定常{len(xs)}点）: p={coef[0]:.3f}±{2*p_se:.3f}(95%) "
              f"C={C:.4f} [{C_lo:.4f},{C_hi:.4f}] R²={r2:.3f}")
        print(f"  p=2 を含む: {hp}   C=sin²(23π/124)={A23:.6f} を含む: {hC}")
        out["H_C"] = {"n_points": len(xs), "p": float(coef[0]), "p_2se": 2 * p_se,
                       "C": C, "C_95": [C_lo, C_hi], "R2": r2,
                       "contains_p2": hp, "contains_A23": hC}
    out["runtime_sec"] = time.time() - t0
    (HERE / "census_longwindow_result_v2.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → census_longwindow_result_v2.json")


if __name__ == "__main__":
    main()
