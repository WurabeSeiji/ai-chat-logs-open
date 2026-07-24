#!/usr/bin/env python3
"""親の不動点反復の途中経緯トレース。

run_n_scaling_lowrank_v1.py の make_parent（158-190行）のループを忠実に
複製し、指定ピッチごとに星座図 PNG と JSON を出力する。

規約（シリーズ内完結の再現性規約）:
  コピー→対照テスト。実行の最後に原本 make_parent を同一乱数シードで
  呼び、最終状態が一致（<1e-10）しなければ終了コード 1 で失敗させる。

反復上限は原本と同じ「1 リスタートあたり 1200 回・リスタート 8 回」を
既定とする（対照テストを厳密に保つため）。上限内で残差 < tol に達した
時点で安定と判定して停止する（原本と同じ、10 回ごとの判定）。

出力先: lowN_metastable_result_v1/parent_iteration_trace_N{n:05d}_seed{s:03d}/
  restart{r:02d}_iter{it:06d}.png / .json   ピッチごとのスナップショット
  index.json                                残差の全履歴と一覧

使い方:
  python3 trace_parent_iteration_v1.py 6 --seed=0 --pitch=100
"""

import argparse
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from run_n_scaling_lowrank_v1 import (
    LowRankSystem,
    make_parent,
    _eigenmode_residual,
    build_edges,
)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def amplitude_phase_groups(v, ea, eb):
    groups = {}
    for e in range(len(v)):
        key = (round(abs(v[e]), 6), round(float(np.angle(v[e])), 6))
        groups.setdefault(key, []).append((int(ea[e]), int(eb[e])))
    return groups


def save_snapshot(out_dir, restart, it, sys_lr, v, ea, eb):
    sys_lr.set_theta(np.angle(v))
    residual = _eigenmode_residual(sys_lr, v)
    sigmas = sys_lr.sigma_spectrum()
    groups = amplitude_phase_groups(v, ea, eb)

    tag = f"restart{restart:02d}_iter{it:06d}"
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    for r in sorted({k[0] for k in groups}):
        ax.add_patch(plt.Circle((0, 0), r, fill=False, ls="--", color="gray", lw=0.6))
    for (r, th), edges in groups.items():
        x, y = r * np.cos(th), r * np.sin(th)
        ax.plot([0, x], [0, y], color="tab:blue", lw=1.0)
        ax.plot(x, y, "o", ms=9, color="tab:red")
        if len(edges) > 1:
            ax.annotate(
                f"x{len(edges)}", (x, y),
                textcoords="offset points", xytext=(9, 7), fontsize=12,
            )
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(0, color="k", lw=0.5)
    lim = 1.45 * float(np.max(np.abs(v)))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_title(f"restart {restart}  iter {it}  residual {residual:.2e}")
    fig.savefig(os.path.join(out_dir, f"{tag}.png"), dpi=100, bbox_inches="tight")
    plt.close(fig)

    payload = {
        "restart": restart,
        "iter": it,
        "residual": float(residual),
        "sigma_max": float(sigmas[0]) if len(sigmas) else None,
        "distinct_amplitudes": len({k[0] for k in groups}),
        "group_count": len(groups),
        "v_real": [float(x) for x in v.real],
        "v_imag": [float(x) for x in v.imag],
        "amplitude_phase_groups": [
            {"amplitude": r, "phase_rad": th, "count": len(edges), "edges": edges}
            for (r, th), edges in sorted(groups.items())
        ],
    }
    with open(os.path.join(out_dir, f"{tag}.json"), "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    return residual, len({k[0] for k in groups})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--pitch", type=int, default=100)
    ap.add_argument("--max-iters", type=int, default=1200)
    ap.add_argument("--restarts", type=int, default=8)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--beta", type=float, default=0.5)
    args = ap.parse_args()

    sys_lr = LowRankSystem(args.n)
    ea, eb = build_edges(args.n)
    beta_tag = "" if args.beta == 0.5 else f"_beta{args.beta:.10f}".rstrip("0")
    out_dir = os.path.join(
        BASE_DIR,
        "lowN_metastable_result_v1",
        f"parent_iteration_trace_N{args.n:05d}_seed{args.seed:03d}{beta_tag}",
    )
    os.makedirs(out_dir, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    history = []
    best = (None, np.inf)

    # ---- ここから make_parent（158-190行）の忠実な複製。スナップショット行のみ追加 ----
    final_v = None
    for restart in range(args.restarts):
        theta = rng.uniform(0.0, 2.0 * np.pi, sys_lr.m)
        v = None
        converged = False
        for it in range(args.max_iters):
            sys_lr.set_theta(theta)
            ev, EV = np.linalg.eig(sys_lr.J @ sys_lr.G)
            idx = int(np.argmin(ev.imag))
            v = sys_lr.w(EV[:, idx].astype(complex))
            v = v / np.linalg.norm(v)
            theta_new = np.angle(v)
            mix = (1.0 - args.beta) * np.exp(1j * theta) + args.beta * np.exp(1j * theta_new)
            theta = np.angle(mix)
            if it % 10 == 9:
                sys_lr.set_theta(np.angle(v))
                res_now = _eigenmode_residual(sys_lr, v)
                if res_now < args.tol:
                    converged = True
            if it % args.pitch == 0 or converged or it == args.max_iters - 1:
                residual, n_amp = save_snapshot(out_dir, restart, it, sys_lr, v, ea, eb)
                history.append(
                    {"restart": restart, "iter": it, "residual": residual,
                     "distinct_amplitudes": n_amp}
                )
                print(f"リスタート{restart} 反復{it:>5} 残差={residual:.3e} "
                      f"異なる振幅数={n_amp}")
            if converged:
                break
        sys_lr.set_theta(np.angle(v))
        residual = _eigenmode_residual(sys_lr, v)
        if residual < best[1]:
            best = (v, residual)
        if residual < args.tol:
            final_v = v
            break
    v_final, res_final = best
    # ---- 複製ここまで ----

    # 対照テスト: 原本 make_parent を同一シードで実行し最終状態を突き合わせ
    rng_ref = np.random.default_rng(args.seed)
    sys_ref = LowRankSystem(args.n)
    v_ref, res_ref, _ = make_parent(
        sys_ref, rng_ref, iters=args.max_iters, beta=args.beta,
        tol=args.tol, restarts=args.restarts
    )
    diff = float(np.linalg.norm(v_final - v_ref))
    ok = diff < 1e-10

    with open(os.path.join(out_dir, "index.json"), "w", encoding="utf-8") as fh:
        json.dump(
            {
                "n": args.n, "m": sys_lr.m, "seed": args.seed,
                "pitch": args.pitch, "max_iters": args.max_iters,
                "restarts": args.restarts, "tol": args.tol,
                "final_residual": float(res_final),
                "control_test_diff_vs_original": diff,
                "control_test_pass": ok,
                "history": history,
            },
            fh, ensure_ascii=False, indent=2,
        )

    print(f"\n最終残差 = {res_final:.3e}")
    print(f"対照テスト: 原本 make_parent との差 = {diff:.3e} → "
          f"{'合格' if ok else '不合格'}")
    print(f"出力先: {out_dir}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
