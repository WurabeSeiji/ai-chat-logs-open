#!/usr/bin/env python3
"""GENESIS v1: 創成通し実験——真空→インフレーション→時空形成→物質誕生の一本走行

これまでの全検証（段階1-3）が支える本番実験。一つの連続シミュレーションで:
    初期状態 = コヒーレント凝縮真空（control親、偶スライス k=2）
             ＋ 極微の奇数種（白色セクター状態 ×δ、奇スライス k=1）
    力学 = v2エンジン（共有O線形部＋媒介非弾性頂点、常時ON・IF文なし）
    を T=4000 走らせ、潜伏（〜288）→ バースト（crossing≈1166）→
    三方向準安定（〜4000）の全史で、物質生成が二因子則
    （ゲート[凝縮]×種²）に従って時代変調されるかを検定する。

事前登録予言（実行前固定・事後変更禁止）:
    P-G1 時代構造: 生成率 g(t)=d(ln P_odd)/dt の窓平均が
         バースト窓 [400,1100] で 準安定窓 [2000,4000] の 0.5倍未満に抑制され
         （時計の破れ＝位相整合の喪失）、準安定窓で回復する。
    P-G2 法則整合: 潜伏窓 [50,250] の生成率が C·f_seed² （C=5.9、段階3実測）
         の [×1/3,×3] に入る。
    P-G3 和則維持: T=250（潜伏末）と T=3000（準安定）で相棒排他性
         P₃/P₄ ≥ 10（全時代でカスケード優位にならない）。
    P-G4 無擾乱: δ=1e-6 では slice-2 の crossing が 1166±5（背反作用なし）。
    探索記録: δ=1e-2 の第二走行で背反作用（物質生成が膨張史を変えるか——
         crossing 遅延・バースト率変化）を記録する（判定は課さない）。

使い方: python3 run_genesis_v1.py
"""

from __future__ import annotations

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
spec3 = importlib.util.spec_from_file_location(
    "s3g", HERE / "run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = s3
spec3.loader.exec_module(s3)
abl = s3.abl
gen3 = s3.gen3
V2 = s3.VertexEngineV2

T_LONG = 4000
C_REF = 5.883        # 段階3 B3 実測
WIN_LAT = (50, 250)
WIN_BURST = (400, 1100)
WIN_META = (2000, 4000)


def window_rate(ts, lnP, lo, hi):
    m_ = (ts >= lo) & (ts < hi)
    tt = ts[m_].astype(float)
    A = np.vstack([tt, np.ones_like(tt)]).T
    coef, _, _, _ = np.linalg.lstsq(A, lnP[m_], rcond=None)
    return float(coef[0])


def run_genesis(delta, tag):
    n, m, nreg = 5, 10, 5
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r2 = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C0 = np.zeros((m, nreg), complex)
    C0[:, 2] = Z0c
    C0[:, 1] = delta * seed_state
    eng = V2(n, C0, wp0, vertex_on=True)
    p2 = C0[:, 2].real / np.linalg.norm(C0[:, 2].real)
    q2 = C0[:, 2].imag - (C0[:, 2].imag @ p2) * p2
    q2 = q2 / np.linalg.norm(q2)

    f2s = np.zeros(T_LONG)
    fseeds = np.zeros(T_LONG)
    P3s = np.zeros(T_LONG)
    P4s = np.zeros(T_LONG)
    rs = np.zeros(T_LONG)
    crossing = None
    Cflat_prev = eng.C.flatten().copy()
    snap = {}
    for t in range(T_LONG):
        eng.step()
        Z2 = eng.C[:, 2]
        Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
        f2 = float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z2) @ Z2))
        f2s[t] = f2
        if crossing is None and f2 > 0.05:
            crossing = t + 1
        d = eng.diagnostics()
        fseeds[t] = d["f_seed"]
        P = np.sum(np.abs(eng.C) ** 2, axis=0)
        P3s[t], P4s[t] = P[3], P[4]
        cf = eng.C.flatten()
        ip = np.vdot(Cflat_prev, cf)
        ph = ip / abs(ip) if abs(ip) > 0 else 1.0
        rs[t] = float(np.linalg.norm(cf - ph * Cflat_prev))
        Cflat_prev = cf.copy()
        if t + 1 in (250, 3000):
            snap[t + 1] = {"P3_P4": float(P[3] / max(P[4], 1e-300))}
    ts = np.arange(1, T_LONG + 1)
    lnP = np.log(np.maximum(fseeds, 1e-300))
    g_lat = window_rate(ts, lnP, *WIN_LAT)
    g_burst = window_rate(ts, lnP, *WIN_BURST)
    g_meta = window_rate(ts, lnP, *WIN_META)
    return {"delta": delta, "crossing": crossing,
            "g_latency": g_lat, "g_burst": g_burst, "g_metastable": g_meta,
            "f_seed0": float(fseeds[0]), "f_seed_final": float(fseeds[-1]),
            "P3_P4_snapshots": snap,
            "series": {"f2": f2s, "fseed": fseeds, "r": rs}}


def main() -> None:
    t0 = time.time()
    print("GENESIS v1: 真空→インフレーション→時空→物質の通し走行（N=5, T=4000）")
    out = {"criteria": {"P_G1": "g_burst < 0.5*g_metastable",
                         "P_G2": "g_latency / (C*f_seed0^2) in [1/3, 3], C=5.883",
                         "P_G3": "P3/P4 >= 10 at T=250 and T=3000",
                         "P_G4": "crossing in [1161,1171] at delta=1e-6"}}

    g = run_genesis(1e-6, "main")
    fs0 = g["f_seed0"]
    law = C_REF * fs0 ** 2
    ratio_law = g["g_latency"] / law
    pg1 = bool(g["g_burst"] < 0.5 * g["g_metastable"])
    pg2 = bool(1 / 3 <= ratio_law <= 3)
    pg3 = bool(all(s["P3_P4"] >= 10 for s in g["P3_P4_snapshots"].values()))
    pg4 = bool(g["crossing"] is not None and 1161 <= g["crossing"] <= 1171)
    print(f"  膨張史: crossing={g['crossing']}（参照1166） → P-G4={pg4}")
    print(f"  生成率: 潜伏={g['g_latency']:.3e} バースト={g['g_burst']:.3e} "
          f"準安定={g['g_metastable']:.3e}")
    print(f"  P-G1 バースト抑制: {g['g_burst']:.3e} < 0.5×{g['g_metastable']:.3e} → {pg1}")
    print(f"  P-G2 法則整合: 潜伏率/法則値={ratio_law:.2f} → {pg2}")
    print(f"  P-G3 和則維持: P3/P4 = "
          f"{ {k: round(v_['P3_P4'],1) for k, v_ in g['P3_P4_snapshots'].items()} } → {pg3}")
    print(f"  f_seed: {fs0:.2e} → {g['f_seed_final']:.2e}"
          f"（{g['f_seed_final']/fs0:.2f}倍）")

    gb = run_genesis(1e-2, "backreaction")
    print(f"  探索（δ=1e-2 背反作用）: crossing={gb['crossing']}（δ=1e-6: {g['crossing']}） "
          f"f_seed {gb['f_seed0']:.1e}→{gb['f_seed_final']:.1e}")

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    ts = np.arange(1, T_LONG + 1)
    axes[0].semilogy(ts, np.clip(g["series"]["f2"], 1e-34, None), color="#7F7F7F", lw=1.0,
                      label="f₂ (inflation of slice-2 vacuum)")
    axes[0].axhline(0.05, color="red", ls=":", lw=0.8)
    axes[0].set_ylabel("f₂ (log)")
    axes[0].legend(fontsize=8)
    axes[1].semilogy(ts, g["series"]["fseed"], color="#4C78A8", lw=1.0,
                      label="f_seed (matter fraction), δ=1e-6")
    axes[1].semilogy(ts, gb["series"]["fseed"], color="#E45756", lw=0.8,
                      label="δ=1e-2 (backreaction run)")
    axes[1].set_ylabel("f_seed (log)")
    axes[1].legend(fontsize=8)
    axes[2].semilogy(ts, np.clip(g["series"]["r"], 1e-17, None), color="black", lw=0.7,
                      label="one-step residual r (clock)")
    axes[2].set_ylabel("r (log)")
    axes[2].set_xlabel("step")
    axes[2].legend(fontsize=8)
    for ax in axes:
        for w, c in ((WIN_LAT, "#2E7D32"), (WIN_BURST, "#B71C1C"), (WIN_META, "#1565C0")):
            ax.axvspan(w[0], w[1], alpha=0.05, color=c)
    axes[0].set_title("GENESIS v1: vacuum → inflation → spacetime → matter (N=5)")
    fig.tight_layout()
    fig.savefig(HERE / "fig_genesis_v1.png", dpi=130)
    plt.close(fig)

    for k_ in ("series",):
        g.pop(k_)
        gb.pop(k_)
    out.update({"main": g, "backreaction": gb,
                 "P_G1": pg1, "P_G2": pg2, "P_G3": pg3, "P_G4": pg4,
                 "all_pass": bool(pg1 and pg2 and pg3 and pg4),
                 "runtime_sec": time.time() - t0})
    print(f"\n判定: {'ALL PASS' if out['all_pass'] else '不成立あり——反証として記録'}")
    (HERE / "genesis_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
