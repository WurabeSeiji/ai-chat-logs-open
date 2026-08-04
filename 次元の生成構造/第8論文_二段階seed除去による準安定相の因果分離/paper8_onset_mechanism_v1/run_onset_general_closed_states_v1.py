#!/usr/bin/env python3
"""ONS-3: 一般零閉塞状態の直接注入——単独セクターを超えた代表性の検証

目的（2026-08-05 木原氏合意・開始様式判別論文の機構検証その3）:
    力学は非線形なので、単独倍音セクター42本の即時型から、混合状態や
    閉塞多様体上の任意点の挙動は従わない。二つのアンサンブルを直接注入する。

アンサンブル（各 N につき 10 状態ずつ）:
    A. 直接一般零閉塞状態: M×2 ガウス白色雑音の QR 直交二列から
       Z=(u1+i u2)/√2（rng 94000+j）。閉塞 Z^T Z=0 は恒等成立、‖Z‖=1。
       これは M 次元閉塞多様体のガウス誘導測度によるランダム標本である。
    B. 許容セクター混合: v3 生成器（確定シード N=5:2/N=40:1）の許容
       セクター C[:,k] を複素ガウス係数 α_k（rng 95000+j）で重ね、
       閉塞多様体へ最小変更射影（[X Y] の SVD 特異値均等化）して規格化。
       射影距離を記録する。

判定基準（プロファイル測定と同一・実行前固定）:
    潜伏 t_launch = max{t: f<1e-20}、幾何級数的= 窓 [1e-20,1e-2] で
    6桁以上かつ R²>0.99、即時型= f(2)>1e-6。T=6000、f は毎 step。

予言（実行前固定・事後変更禁止）:
    全 20 状態が即時型。潜伏バースト型・その他は現れない。
    これが破れた場合、その状態と一段残差の照合を追試し、論文の主張を
    実測に合わせて修正する。

再現性: abl・生成器v3 とも read-only import（SHA-256 記録）。
    乱数シードは 94000+j / 95000+j / 注入補助 96000+series。

使い方: python3 run_onset_general_closed_states_v1.py <N>
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = (REPO / "次元の生成構造" / "make_parent_white_managed_v1"
        / "make_parent_white_harmonics_n_only_v3.py")
SEEDS = {5: 2, 40: 1}
N_DRAWS = 10
T_LONG = 6000
FLOOR = 1e-20
FIT_HI = 1e-2
R2_MIN = 0.99
DECADES_MIN = 6.0
IMMEDIATE_F2 = 1e-6
LATENCY_MIN = 100

spec = importlib.util.spec_from_file_location("abl_ons3", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen3_ons3", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def closure(Z):
    return complex(Z @ Z)


def project_closed(Z):
    """[X Y] の SVD 特異値均等化による閉塞多様体への最小変更射影。"""
    B = np.column_stack([Z.real, Z.imag])
    U, s, Vt = np.linalg.svd(B, full_matrices=False)
    s_eq = np.sqrt(np.mean(s ** 2))
    Bp = U @ np.diag([s_eq, s_eq]) @ Vt
    Zp = Bp[:, 0] + 1j * Bp[:, 1]
    dist = float(np.linalg.norm(Bp - B))
    return Zp / np.linalg.norm(Zp), dist


def profile(n, v0, wp):
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(v0))
    p = v0.real / np.linalg.norm(v0.real)
    q = v0.imag - (v0.imag @ p) * p
    q = q / np.linalg.norm(q)

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    Z = v0.copy()
    fs = np.zeros(T_LONG + 1)
    fs[0] = fval(Z)
    for t in range(1, T_LONG + 1):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        fs[t] = fval(Z)
    return fs


def classify(fs):
    crossing = next((t for t, f in enumerate(fs) if f > 0.05), None)
    below = np.nonzero(fs < FLOOR)[0]
    t_launch = int(below.max()) if below.size else 0
    lim = crossing if crossing is not None else len(fs) - 1
    win = [(t, f) for t, f in enumerate(fs[:lim + 1]) if FLOOR <= f <= FIT_HI]
    fit = {"n_points": len(win), "decades": 0.0, "rate_per_step": None, "r2": None}
    geometric = False
    if len(win) >= 5:
        tt = np.array([w[0] for w in win], float)
        lf = np.log(np.array([w[1] for w in win]))
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, lf, rcond=None)
        pred = A @ coef
        ss_res = float(np.sum((lf - pred) ** 2))
        ss_tot = float(np.sum((lf - lf.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        decades = float((lf.max() - lf.min()) / np.log(10))
        fit = {"n_points": len(win), "decades": decades,
               "rate_per_step": float(coef[0]), "r2": float(r2)}
        geometric = decades >= DECADES_MIN and r2 > R2_MIN
    immediate = bool(fs[2] > IMMEDIATE_F2)
    if t_launch > LATENCY_MIN and geometric:
        cls = "潜伏バースト型（幾何級数的）"
    elif immediate:
        cls = "即時型（有限離脱）"
    else:
        cls = "その他（生データ参照）"
    return {"crossing": crossing, "t_launch": t_launch,
            "f2": float(fs[2]), "fit": fit, "geometric": bool(geometric),
            "immediate": immediate, "class": cls}


def main() -> None:
    t0 = time.time()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    seed = SEEDS[n]
    print(f"ONS-3 一般零閉塞状態注入 N={n}（アンサンブル各{N_DRAWS}）")
    print(f"  import: ABL {sha256(ABL)[:16]}…  GEN3 {sha256(GEN3)[:16]}…")
    m = n * (n - 1) // 2

    r = gen3.make_parent(n, seed=seed)
    C = np.fft.fft(r.relation_waves, axis=1) / n
    allowed = [k for k in range(n) if (2 * k) % n != 0]

    series = {}
    idx = 0
    for j in range(N_DRAWS):
        rng = np.random.default_rng(94000 + j)
        A = rng.normal(size=(m, 2))
        Q, _ = np.linalg.qr(A)
        Z = (Q[:, 0] + 1j * Q[:, 1]) / np.sqrt(2)
        cl = closure(Z)
        wp = np.random.default_rng(96000 + idx).normal(size=m)
        fs = profile(n, Z / np.linalg.norm(Z), wp)
        series[f"A{j}"] = {"kind": "direct_qr", "closure_abs": abs(cl),
                            **classify(fs),
                            "f_sampled": [[int(t), float(fs[t])]
                                           for t in range(0, len(fs), 10)]}
        idx += 1
    for j in range(N_DRAWS):
        rng = np.random.default_rng(95000 + j)
        alpha = rng.normal(size=len(allowed)) + 1j * rng.normal(size=len(allowed))
        Zmix = np.zeros(m, complex)
        for a, k in zip(alpha, allowed):
            Zmix += a * C[:, k]
        cl_pre = closure(Zmix / np.linalg.norm(Zmix))
        Z, dist = project_closed(Zmix)
        cl_post = closure(Z)
        wp = np.random.default_rng(96000 + idx).normal(size=m)
        fs = profile(n, Z, wp)
        series[f"B{j}"] = {"kind": "sector_mixture", "closure_abs_pre": abs(cl_pre),
                            "closure_abs_post": abs(cl_post), "repair_dist": dist,
                            **classify(fs),
                            "f_sampled": [[int(t), float(fs[t])]
                                           for t in range(0, len(fs), 10)]}
        idx += 1

    counts = {}
    for s in series.values():
        counts[s["class"]] = counts.get(s["class"], 0) + 1
    print(f"  分類: {counts}")
    for name, s in series.items():
        print(f"    {name} [{s['kind']}]: {s['class']} crossing={s['crossing']} "
              f"f(2)={s['f2']:.2e}")

    out = {"N": n, "seed": seed, "n_draws": N_DRAWS,
           "imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)},
           "criteria": {"FLOOR": FLOOR, "FIT_HI": FIT_HI, "R2_MIN": R2_MIN,
                         "DECADES_MIN": DECADES_MIN, "IMMEDIATE_F2": IMMEDIATE_F2,
                         "LATENCY_MIN": LATENCY_MIN, "T_LONG": T_LONG,
                         "prediction": "all 20 states immediate"},
           "series": series, "class_counts": counts,
           "runtime_sec": time.time() - t0}
    (HERE / f"onset_general_closed_N{n:05d}_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 5))
    for name, s in series.items():
        arr = np.array(s["f_sampled"])
        color = "#4C78A8" if s["kind"] == "direct_qr" else "#E45756"
        ax.semilogy(arr[:, 0], np.clip(arr[:, 1], 1e-34, None), color=color, lw=0.6)
    ax.axhline(0.05, color="red", ls=":", lw=0.8)
    ax.semilogy([], [], color="#4C78A8", label="A: direct QR closed states")
    ax.semilogy([], [], color="#E45756", label="B: sector mixtures (projected)")
    ax.set_xlim(0, T_LONG)
    ax.set_xlabel("step")
    ax.set_ylabel("f (log)")
    ax.set_title(f"N={n} ONS-3: general zero-closure states")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_onset_general_closed_N{n:05d}.png", dpi=130)
    plt.close(fig)
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
