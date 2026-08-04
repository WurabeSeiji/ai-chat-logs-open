#!/usr/bin/env python3
"""E-M9r-P：立ち上がりプロファイル測定 v1——幾何級数的成長か即時有限離脱かの判定

目的（2026-08-04 木原氏指示）:
    白色倍音セクターの立ち上がり（crossing 11〜16 step）が
    幾何級数的（指数）成長なのか、即時の有限角離脱なのかを、
    対照（潜伏付きインフレーション）と同一の測定で判別する。

系列（v3 生成器・シードは確定値 N=5:2 / N=40:1）:
    control（既存 build_init）/ v3 親ベクトル / v3 正当セクター全本

測定: f(t) を毎 step 記録（T_LONG=6000）。
判定基準（実行前固定・事後変更禁止）:
    潜伏 t_launch = max{t: f(t) < 1e-20}（floor滞在の最終時刻。初期からf≥1e-20なら0）
    バースト幅 = crossing − t_launch
    指数成長判定: f ∈ [1e-20, 1e-2] の窓で ln f vs t を最小二乗、
        判定=「窓が6桁以上 かつ R² > 0.99」→ 幾何級数的
    即時型判定: f(2) > 1e-6（2 step 以内の有限離脱）
    分類: 潜伏バースト型 = t_launch > 100 かつ 幾何級数的
          即時型 = f(2) > 1e-6
          （両方満たさない場合は「その他」として生データごと記録）

再現性: v3 生成器・abl とも read-only import・SHA-256 記録。乱数は
    生成器シード（2/1）と wp 用 rng（対照=build_init 内固定、親=91000、
    セクター=92000+k）のみ。全系列の f 時系列を JSON に保存
    （早期 0..100 は毎step、以降は 10 step 間引き）。図も同時出力
    （semilogy、対照=灰 #7F7F7F・親=黒・セクター=青系、lw0.8、dpi130）。

使い方: python3 run_paper8_em9r_profile_v1.py <N>
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
T_LONG = 6000
FLOOR = 1e-20
FIT_HI = 1e-2
R2_MIN = 0.99
DECADES_MIN = 6.0
IMMEDIATE_F2 = 1e-6
LATENCY_MIN = 100

spec = importlib.util.spec_from_file_location("abl_pf", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen3_pf", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


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
    # 指数成長フィット窓
    lim = crossing if crossing is not None else len(fs) - 1
    win = [(t, f) for t, f in enumerate(fs[:lim + 1]) if FLOOR <= f <= FIT_HI]
    fit = {"n_points": len(win), "decades": 0.0, "rate": None, "r2": None}
    geometric = False
    if len(win) >= 5:
        tt = np.array([w[0] for w in win], float)
        lf = np.log(np.array([w[1] for w in win]))
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, res, _, _ = np.linalg.lstsq(A, lf, rcond=None)
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
            "burst_width": (crossing - t_launch) if crossing is not None else None,
            "f2": float(fs[2]), "fit": fit, "geometric": bool(geometric),
            "immediate": immediate, "class": cls}


def pack_series(fs):
    early = [float(x) for x in fs[:101]]
    rest = [[int(t), float(fs[t])] for t in range(110, len(fs), 10)]
    return {"early_0_100_every_step": early, "sampled_every_10": rest}


def main() -> None:
    t0 = time.time()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    seed = SEEDS[n]
    print(f"E-M9r-P プロファイル測定 N={n}（v3 seed={seed}）")
    print(f"  import: ABL {sha256(ABL)[:16]}…  GEN3 {sha256(GEN3)[:16]}…")
    r = gen3.make_parent(n, seed=seed)
    W = r.relation_waves
    m = W.shape[0]
    C = np.fft.fft(W, axis=1) / n

    series = {}
    _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(n, False)
    fs = profile(n, Z0, wp0.copy())
    series["control"] = {"kind": "control", **classify(fs), "f_series": pack_series(fs)}
    print(f"  control: {series['control']['class']} crossing={series['control']['crossing']} "
          f"潜伏={series['control']['t_launch']} "
          f"rate={series['control']['fit'].get('rate_per_step')} "
          f"R²={series['control']['fit'].get('r2')} 桁数={series['control']['fit'].get('decades'):.1f}"
          if series['control']['fit'].get('r2') is not None else "  control: fit不可")

    wp = np.random.default_rng(91000).normal(size=m)
    fs = profile(n, r.parent_vector / np.linalg.norm(r.parent_vector), wp)
    series["parent_vector"] = {"kind": "parent", **classify(fs), "f_series": pack_series(fs)}
    sp = series["parent_vector"]
    print(f"  parent : {sp['class']} crossing={sp['crossing']} 潜伏={sp['t_launch']} f(2)={sp['f2']:.2e}")

    for k in range(n):
        if (2 * k) % n == 0:
            continue
        Zk = C[:, k]
        wp = np.random.default_rng(92000 + k).normal(size=m)
        fs = profile(n, Zk / np.linalg.norm(Zk), wp)
        series[f"k{k}"] = {"kind": "sector", "k": k, **classify(fs),
                            "f_series": pack_series(fs)}
    secs = [s for s in series.values() if s["kind"] == "sector"]
    cls_count = {}
    for s in secs:
        cls_count[s["class"]] = cls_count.get(s["class"], 0) + 1
    f2s = [s["f2"] for s in secs]
    print(f"  sectors({len(secs)}本): 分類={cls_count} f(2)範囲=[{min(f2s):.2e},{max(f2s):.2e}] "
          f"crossing範囲=[{min(s['crossing'] for s in secs)},{max(s['crossing'] for s in secs)}]")

    # 図（semilogy・対照=灰・親=黒・セクター=青系）
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for name, s in series.items():
        fsr = np.array(s["f_series"]["early_0_100_every_step"]
                       + [x[1] for x in s["f_series"]["sampled_every_10"]])
        tsr = np.array(list(range(101)) + [x[0] for x in s["f_series"]["sampled_every_10"]])
        if s["kind"] == "control":
            color, lw, lab = "#7F7F7F", 1.2, "control"
        elif s["kind"] == "parent":
            color, lw, lab = "black", 1.2, "parent_vector (v3)"
        else:
            color, lw, lab = "#4C78A8", 0.6, None
        ax.semilogy(tsr, np.clip(fsr, 1e-34, None), color=color, lw=lw, label=lab)
    ax.axhline(0.05, color="red", ls=":", lw=0.8, label="crossing threshold f=0.05")
    ax.set_xlim(0, T_LONG)
    ax.set_xlabel("step (absolute)")
    ax.set_ylabel("f = 1 - E_P1 (log)")
    ax.set_title(f"N={n} E-M9r-P: rise profiles (control / v3 parent / white sectors, log)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_em9r_profile_N{n:05d}.png", dpi=130)
    plt.close(fig)

    out = {"N": n, "seed": seed,
           "imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)},
           "criteria": {"FLOOR": FLOOR, "FIT_HI": FIT_HI, "R2_MIN": R2_MIN,
                         "DECADES_MIN": DECADES_MIN, "IMMEDIATE_F2": IMMEDIATE_F2,
                         "LATENCY_MIN": LATENCY_MIN, "T_LONG": T_LONG},
           "series": series, "runtime_sec": time.time() - t0}
    (HERE / f"paper8_em9r_profile_N{n:05d}_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s): JSON + fig_em9r_profile_N{n:05d}.png")


if __name__ == "__main__":
    main()
