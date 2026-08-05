#!/usr/bin/env python3
"""GATE-0: 相互作用ゲートの計装実験——力学無改変・読出しのみ

目的（2026-08-05 木原氏合意・設計方針_v1.md 段階0）:
    「相互作用はいつ起こり得るか」を設計せず、既存軌道の状態から
    ゲート因子（凝縮＝共有フレームと時計の存在）を毎step読出し、
    自己ゲートの予言を検定する。力学は一切変更しない。

観測量（全て O(M)/step）:
    r(t) = min_φ ‖Z(t+1)−e^{iφ}Z(t)‖   一段残差（時計凝縮の計器）
    O(t) = tr(P(t)P(t−Δ))/2, Δ=10       平面安定度（舞台の計器）
    γ(t) = Σ_v|S_v|²/((N−1)Σ_e|Z_e|²)   頂点コヒーレンス（探索的・合否外）

予言（実行前固定・事後変更禁止）:
    P1 潜伏期（control, t≤t_launch）: 中央値 r ≤ 1e-10 かつ O ≥ 1−1e-6
       ——親凝縮体は完全なゲート開（ただし種ゼロ）
    P2 バースト窓（t_launch<t≤crossing+500）: max r ≥ 1e-3 かつ min O ≤ 0.99
       ——急拡大中はレジスタが滑り、ゲートが閉じる
    P3 準安定末尾（最後の1000 step）: 中央値 r < 0.3×(バースト窓 max r)
       かつ 中央値 O ≥ 0.999 ——安定時空でゲートが再び開く
    P4 白色セクター: r(0)=O(1e-2) から出発し、末尾1000 step では
       P3 と同じ開ゲート署名（中央値 O ≥ 0.999）に合流する
       ——ゲートの開放は出自に依らず、準安定時空の性質である
    P5 白色親 N=5（安定平衡・不活性）: 全期間 中央値 r ≤ 1e-10, O ≥ 1−1e-6
       ——ゲートは開きっぱなしだが種がなく、何も起こらない

系列（v3生成器シード N=5:2 / N=40:1、注入補助は既存実験と同一シード）:
    N=5:  control / white_parent / k1..k4 / general_A0
    N=40: control / white_parent / k1 / general_A0

再現性: abl・gen3 とも read-only import（SHA-256記録）。頂点集約は
    triu_indices(n,1) の辺順序（gen3 と同一、E-M9r 注入実験で整合実証済み）。
    T=6000、r・γ は毎step、O は10step毎。全系列 JSON 保存＋図出力。

使い方: python3 run_stage0_gate_instrumentation_v1.py <N>
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
SERIES_DIR = HERE.parent
PAPER8 = SERIES_DIR / "第8論文_二段階seed除去による準安定相の因果分離"
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = SERIES_DIR / "make_parent_white_managed_v1" / "make_parent_white_harmonics_n_only_v3.py"
SEEDS = {5: 2, 40: 1}
T_LONG = 6000
DELTA_O = 10
TAIL = 1000
BURST_PAD = 500

spec = importlib.util.spec_from_file_location("abl_g0", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen3_g0", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def orth_plane(Z):
    p = Z.real / np.linalg.norm(Z.real)
    q = Z.imag - (Z.imag @ p) * p
    nq = np.linalg.norm(q)
    if nq < 1e-300:
        return np.column_stack([p, np.zeros_like(p)])
    return np.column_stack([p, q / nq])


def run_series(n, Z0, wp, ia, ib, deg):
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(Z0))
    Z = Z0 / np.linalg.norm(Z0)
    p0 = Z.real / np.linalg.norm(Z.real)
    q0 = Z.imag - (Z.imag @ p0) * p0
    q0 = q0 / np.linalg.norm(q0)

    def fval(Zv):
        Zp = Zv - p0 * (p0 @ Zv) - q0 * (q0 @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    def gamma(Zv):
        nverts = n
        S = np.zeros(nverts, complex)
        np.add.at(S, ia, Zv)
        np.add.at(S, ib, Zv)
        return float(np.sum(np.abs(S) ** 2) / (deg * np.sum(np.abs(Zv) ** 2)))

    rs = np.zeros(T_LONG)
    gs = np.zeros(T_LONG + 1)
    fs = np.zeros(T_LONG + 1)
    Os = []
    gs[0] = gamma(Z)
    fs[0] = fval(Z)
    B_prev = orth_plane(Z)
    for t in range(1, T_LONG + 1):
        Znew, wp = abl.evolve(sys_lr, Z, wp)
        ip = np.conj(Z) @ Znew
        ph = ip / abs(ip) if abs(ip) > 0 else 1.0
        rs[t - 1] = float(np.linalg.norm(Znew - ph * Z))
        Z = Znew
        gs[t] = gamma(Z)
        fs[t] = fval(Z)
        if t % DELTA_O == 0:
            B = orth_plane(Z)
            M12 = B_prev.T @ B
            Os.append([t, float(np.sum(M12 ** 2) / 2.0)])
            B_prev = B
    crossing = next((t for t, f in enumerate(fs) if f > 0.05), None)
    below = np.nonzero(fs < 1e-20)[0]
    t_launch = int(below.max()) if below.size else 0
    return {"r": rs, "gamma": gs, "f": fs, "O": np.array(Os),
            "crossing": crossing, "t_launch": t_launch}


def verdicts(res, kind):
    r, O = res["r"], res["O"]
    t_launch, crossing = res["t_launch"], res["crossing"]
    out = {}
    tail_r = float(np.median(r[-TAIL:]))
    tail_O = float(np.median(O[O[:, 0] > T_LONG - TAIL][:, 1])) if len(O) else None
    out["tail_median_r"] = tail_r
    out["tail_median_O"] = tail_O
    if kind == "control":
        lat = r[:max(t_launch, 1)]
        latO = O[O[:, 0] <= max(t_launch, DELTA_O)][:, 1]
        burst_hi = min((crossing or T_LONG) + BURST_PAD, T_LONG)
        burst_r = r[t_launch:burst_hi]
        burst_O = O[(O[:, 0] > t_launch) & (O[:, 0] <= burst_hi)][:, 1]
        out["latency_median_r"] = float(np.median(lat))
        out["latency_min_O"] = float(np.min(latO)) if len(latO) else None
        out["burst_max_r"] = float(np.max(burst_r)) if len(burst_r) else None
        out["burst_min_O"] = float(np.min(burst_O)) if len(burst_O) else None
        out["P1"] = bool(out["latency_median_r"] <= 1e-10
                         and (out["latency_min_O"] or 0) >= 1 - 1e-6)
        out["P2"] = bool((out["burst_max_r"] or 0) >= 1e-3
                         and (out["burst_min_O"] or 1) <= 0.99)
        out["P3"] = bool(tail_r < 0.3 * (out["burst_max_r"] or 1)
                         and (tail_O or 0) >= 0.999)
    elif kind == "sector" or kind == "general":
        out["r0"] = float(r[0])
        out["P4"] = bool((tail_O or 0) >= 0.999)
    elif kind == "inert_parent":
        out["all_median_r"] = float(np.median(r))
        out["all_min_O"] = float(np.min(O[:, 1]))
        out["P5"] = bool(out["all_median_r"] <= 1e-10 and out["all_min_O"] >= 1 - 1e-6)
    elif kind == "burst_parent":
        burst_hi = min((crossing or T_LONG) + BURST_PAD, T_LONG)
        out["burst_max_r"] = float(np.max(r[t_launch:burst_hi]))
        out["P2"] = bool(out["burst_max_r"] >= 1e-3)
        out["P3"] = bool(tail_r < 0.3 * out["burst_max_r"] and (tail_O or 0) >= 0.999)
    return out


def main() -> None:
    t0 = time.time()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    seed = SEEDS[n]
    m = n * (n - 1) // 2
    ia, ib = np.triu_indices(n, k=1)
    deg = n - 1
    print(f"GATE-0 計装 N={n}（v3 seed={seed}）  ABL {sha256(ABL)[:16]}…  GEN3 {sha256(GEN3)[:16]}…")

    r = gen3.make_parent(n, seed=seed)
    C = np.fft.fft(r.relation_waves, axis=1) / n
    vp = r.parent_vector / np.linalg.norm(r.parent_vector)

    series = {}
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    series["control"] = ("control", Z0c, wp0.copy())
    kind_p = "inert_parent" if n == 5 else "burst_parent"
    series["white_parent"] = (kind_p, vp, np.random.default_rng(91000).normal(size=m))
    ks = [k for k in range(n) if (2 * k) % n != 0]
    ks = ks[:4] if n == 5 else ks[:1]
    for k in ks:
        Zk = C[:, k] / np.linalg.norm(C[:, k])
        series[f"k{k}"] = ("sector", Zk, np.random.default_rng(92000 + k).normal(size=m))
    rngA = np.random.default_rng(94000)
    A = rngA.normal(size=(m, 2))
    Q, _ = np.linalg.qr(A)
    Zg = (Q[:, 0] + 1j * Q[:, 1]) / np.sqrt(2)
    series["general_A0"] = ("general", Zg, np.random.default_rng(96000).normal(size=m))

    results = {}
    for name, (kind, Z0, wp) in series.items():
        res = run_series(n, Z0, wp, ia, ib, deg)
        vd = verdicts(res, kind)
        results[name] = {"kind": kind, "crossing": res["crossing"],
                          "t_launch": res["t_launch"], **vd,
                          "r_sampled": [[int(t), float(res["r"][t])]
                                         for t in range(0, T_LONG, 10)],
                          "gamma_sampled": [[int(t), float(res["gamma"][t])]
                                             for t in range(0, T_LONG + 1, 10)],
                          "O_series": [[int(a), float(b)] for a, b in res["O"]]}
        flags = {k: v for k, v in vd.items() if k.startswith("P")}
        print(f"  {name} [{kind}]: crossing={res['crossing']} 判定={flags} "
              f"tail_r={vd['tail_median_r']:.2e} tail_O={vd['tail_median_O']}")
        results[name]["_res"] = res

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    colors = {"control": "#7F7F7F", "white_parent": "black", "general_A0": "#E45756"}
    for name, d in results.items():
        res = d.pop("_res")
        c = colors.get(name, "#4C78A8")
        lw = 1.2 if name in colors else 0.6
        tt = np.arange(T_LONG)
        axes[0].semilogy(tt, np.clip(res["r"], 1e-17, None), color=c, lw=lw,
                          label=name if name in colors or name == "k1" else None)
        axes[1].plot(res["O"][:, 0], res["O"][:, 1], color=c, lw=lw)
        axes[2].plot(np.arange(0, T_LONG + 1), res["gamma"], color=c, lw=lw)
    axes[0].set_ylabel("one-step residual r (log)")
    axes[0].legend(fontsize=8)
    axes[1].set_ylabel("plane stability O")
    axes[1].set_ylim(0.0, 1.02)
    axes[2].set_ylabel("vertex coherence γ")
    axes[2].set_xlabel("step")
    axes[0].set_title(f"N={n} GATE-0: clock (r), stage (O), mediation field (γ)")
    fig.tight_layout()
    fig.savefig(HERE / f"fig_gate0_N{n:05d}.png", dpi=130)
    plt.close(fig)

    out = {"N": n, "seed": seed,
           "imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)},
           "criteria": {"T_LONG": T_LONG, "DELTA_O": DELTA_O, "TAIL": TAIL,
                         "BURST_PAD": BURST_PAD,
                         "predictions": ["P1 latency: med r<=1e-10 & O>=1-1e-6",
                                          "P2 burst: max r>=1e-3 & min O<=0.99",
                                          "P3 tail: med r<0.3*burst max & med O>=0.999",
                                          "P4 sectors/general: tail med O>=0.999",
                                          "P5 inert parent: all med r<=1e-10 & min O>=1-1e-6"]},
           "series": results, "runtime_sec": time.time() - t0}
    (HERE / f"gate0_result_N{n:05d}_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
