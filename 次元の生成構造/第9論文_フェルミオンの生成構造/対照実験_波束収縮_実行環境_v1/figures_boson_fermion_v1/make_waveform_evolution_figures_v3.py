#!/usr/bin/env python3
"""論文用図版 v3: 原本論文と同一形式の波形発展図

原本 20260713 エンジンの R070_waveform_evolution 図（4×2グリッド、rho_chi/max、
凡例に L と N_eff、collision=0,1,2,3,5,10,20,42）と同じ形式・同じ関数で、

  図F fig_waveform_evolution_fermionic_v3  R = R₁₃₇厳密値（フェルミオン型: 局在が往復移乗）
  図B fig_waveform_evolution_bosonic_v3    R = R₁₃₇, B=基底波（ボゾン型: 両波とも非局在、同一相互作用でも移乗なし）\n  図C fig_waveform_evolution_no_exchange_control_v3  R = 1.0（補助対照: 交換なし条件）

を生成する。状態生成・発展・観測量・描画形式はすべて原本エンジンの
recursive_snapshot_states / 描画ブロック（1329-1347行）の引用であり、
初期条件は両図とも同一（A=基底波 n=1、B=奇数63次カーネル）。
"""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ENGINE = HERE.parent / "20260713" / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py"

spec = importlib.util.spec_from_file_location("engine", ENGINE)
eng = importlib.util.module_from_spec(spec)
sys.modules["engine"] = eng
spec.loader.exec_module(eng)

plt = eng.plt
plt.rcParams["font.family"] = "Hiragino Sans"

R137 = 0.6971778791282474
EVOLUTION_COLLISIONS = [0, 1, 2, 3, 5, 10, 20, 42]   # 原本と同一


def make_figure(r_value: float, r_label: str, suptitle: str, out_stem: str, n_b: int = 63) -> None:
    params = eng.Params()
    chi, _ = eng.make_grids(params)
    x = chi / math.pi
    snapshots = eng.recursive_snapshot_states(params, r_value, EVOLUTION_COLLISIONS, n_b=n_b)
    # ---- 以下、原本エンジンの R070_waveform_evolution 描画ブロックの引用 ----
    fig, axes = plt.subplots(4, 2, figsize=(12, 12), sharex=True, sharey=True, constrained_layout=True)
    for ax, collision_index in zip(axes.flatten(), EVOLUTION_COLLISIONS):
        snap = snapshots[collision_index]
        a_metrics = snap["A_channel"]
        b_metrics = snap["B_channel"]
        rho_a = snap["rho_A"] / np.max(snap["rho_A"])
        rho_b = snap["rho_B"] / np.max(snap["rho_B"])
        ax.plot(x, rho_a, label=f"A L={a_metrics['L']:.3g}, N={a_metrics['N_eff']:.3g}")
        ax.plot(x, rho_b, label=f"B L={b_metrics['L']:.3g}, N={b_metrics['N_eff']:.3g}")
        ax.set_title(f"{r_label}, collision={collision_index}")
        ax.set_ylabel("rho_chi / max")
        ax.legend(fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("chi / pi")
    # ---- 引用ここまで（追加はタイトルのみ） ----
    fig.suptitle(suptitle, fontsize=13)
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"{out_stem}.{ext}", dpi=160)
    plt.close(fig)
    print(out_stem, "saved")


def main() -> None:
    make_figure(
        R137, "R=R_137",
        "フェルミオン型: 局在（L, N_eff）が A・B 間を往復移乗する（A=基底波, B=奇数63次カーネル, R=0.69717788）",
        "fig_waveform_evolution_fermionic_v3")
    make_figure(
        R137, "R=R_137",
        "ボゾン型: 同一の相互作用条件でも、非局在の波どうしには乗り移る局在が無い（A=基底波, B=基底波, R=0.69717788）",
        "fig_waveform_evolution_bosonic_v3", n_b=1)
    make_figure(
        1.0, "R=1.00",
        "補助対照: 交換なし条件（R=1.00）では局在波束を含む対でも全コマ不変",
        "fig_waveform_evolution_no_exchange_control_v3")


if __name__ == "__main__":
    main()
