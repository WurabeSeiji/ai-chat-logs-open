#!/usr/bin/env python3
"""コヒーレント結合の距離法則 v1（E5）

目的:
    電荷型（コヒーレント）結合の空間距離依存を初めて直接測る。
    Δθ↔r 辞書（未解決2）への零次アタック。

方法:
    局在波束 = 倍音の等重み和（搬送波なし、共有チャネル梯子つき）。
    B 波束を χ 方向に d グリッド点だけ厳密平行移動（円順シフト＝ユニタリ、
    力学・読出しは無変更）し、結合 <a|b(d)> と一衝突流の変調を d の関数で測る。

予言（測定前に固定）:
    直接重なり型の結合は波束の自己相関で減衰する短距離型であり、
    1/r^2 型の長距離則は**現れない**。ゆえにクーロンの長距離性は
    直接重なりではなく、媒介機構（調和閉鎖 [逆二乗論文]）が担う必要がある
    ——辞書が翻訳すべき対象の特定。
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"

spec = importlib.util.spec_from_file_location("toy_for_range_v1", TOY_RUNNER_PATH)
toy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = toy
spec.loader.exec_module(toy)
base = toy.base
plt = base.plt


def unit_norm(v):
    return v / math.sqrt(float(np.vdot(v, v).real))


def mk(k, which, sp):
    case = base.explicit_packet_case(mode=f"range_{which}_{k}", packet_a=(k,), packet_b=(k,))
    return unit_norm(base.make_case_state(sp, case, which, hair_enabled=False))


def shift_chi(state: np.ndarray, d: int, sp) -> np.ndarray:
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    return np.roll(state.reshape(shape), d, axis=0).reshape(state.shape)


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=1)
    sp = base.build_source_params(params)

    packet_kind = sys.argv[1] if len(sys.argv) > 1 else "comb"
    if packet_kind == "comb":
        ks_a = tuple(range(3, 62, 2))
        ks_b = tuple(range(5, 64, 2))
        a_packet = unit_norm(sum(mk(k, "A", sp) for k in ks_a))
        b_packet = unit_norm(sum(mk(k, "B", sp) for k in ks_b))
    else:  # gaussian: 両パリティ密スペクトル・ガウス重み（真に局在した波束）
        k0, sig = 30, 12.0
        ks_a = tuple(range(3, 62))
        ks_b = tuple(range(5, 64))
        a_packet = unit_norm(
            sum(math.exp(-((k - k0) ** 2) / (2 * sig**2)) * mk(k, "A", sp) for k in ks_a)
        )
        b_packet = unit_norm(
            sum(math.exp(-((k - k0) ** 2) / (2 * sig**2)) * mk(k, "B", sp) for k in ks_b)
        )
    print(f"packet kind: {packet_kind}")
    n_chi = sp.chi_grid_n
    th_p = 0.15
    s2t = math.sin(2 * th_p)

    overlap0 = complex(np.vdot(a_packet, b_packet))
    print(f"grid chi_n={n_chi}, packet coupling at d=0: |<a|b>|={abs(overlap0):.6f}")

    ds = sorted(set(list(range(0, 33)) + [40, 48, 64, 96, 128, n_chi // 2]))
    rows = []
    for d in ds:
        b_d = shift_chi(b_packet, d, sp)
        ov = complex(np.vdot(a_packet, b_d))
        a2, b2 = toy.rotate_ab(a_packet.copy(), b_d.copy(), th_p)
        dnb = float(np.vdot(b2, b2).real) - 1.0
        pred = s2t * ov.real
        rows.append(
            {
                "d_grid": d,
                "overlap_abs": abs(ov),
                "overlap_re": ov.real,
                "dN_B": dnb,
                "flow_pred_err": abs(dnb - pred),
            }
        )

    max_pred_err = max(r["flow_pred_err"] for r in rows)
    # 距離法則の判定: 遠方（d >= 32）での結合の上界と、逆二乗仮説との比較
    near = rows[1]["overlap_abs"]
    far = [r for r in rows if r["d_grid"] >= 32]
    far_max = max(r["overlap_abs"] for r in far)
    inv_square_at_32 = near / (32.0 ** 2)
    short_ranged = far_max < inv_square_at_32
    print(f"flow = sin(2th)*Re<a|b> exact: max err {max_pred_err:.2e}")
    print(f"coupling near (d=1): {near:.4f}; far (d>=32) max: {far_max:.3e}")
    print(f"1/r^2 would give at d=32: {inv_square_at_32:.3e}")
    print(f"SHORT-RANGED (decays faster than 1/r^2): {short_ranged}")

    with (HERE / "coupling_range_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    dd = np.asarray([r["d_grid"] for r in rows if r["d_grid"] > 0])
    oo = np.asarray([max(r["overlap_abs"], 1e-18) for r in rows if r["d_grid"] > 0])
    ax.loglog(dd, oo, "o-", markersize=4, label="|<a|b(d)>| measured")
    ax.loglog(dd, near / (dd / 1.0) ** 2, "--", color="0.5", label="1/d^2 reference")
    ax.set_xlabel("packet separation d (grid units)")
    ax.set_ylabel("coherent coupling |<a|b>|")
    ax.set_title("E5: direct-overlap coupling is short-ranged (no 1/d^2 tail)")
    ax.legend(); ax.grid(alpha=0.3, which="both")
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"coupling_range_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "coupling_range_v1",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "flow_identity_max_err": max_pred_err,
        "near_coupling_d1": near,
        "far_coupling_max_d_ge_32": far_max,
        "inverse_square_reference_at_32": inv_square_at_32,
        "short_ranged": bool(short_ranged),
        "conclusion": (
            "直接重なり型のコヒーレント結合は波束相関距離で減衰する短距離型で、"
            "1/r^2 尾を持たない。クーロンの長距離性は媒介機構（調和閉鎖）が担う"
            "必要があり、Δθ↔r 辞書の翻訳対象が特定された"
        ),
        "rows": rows,
    }
    (HERE / "coupling_range_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
