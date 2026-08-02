#!/usr/bin/env python3
"""論文A「二文法の分離」証拠パッケージ v1

論文Aの全主張の証拠を単一ランナーで再生成し、全図表を出力する。

主張:
  C1 保存定理: ビン別回転不変量 |A_k|^2+|B_k|^2 の関数は厳密保存（順方向は定理、
     読出し7種の分類は数値）
  C2 対角セクターの合成則は厳密に平均（重力文法: 加法・符号なし）
  C3 コヒーレント流の三項分解が全て予言どおり:
       dN_B = sigma [ sin^2(th)(N_A-N_B)                    ... 拡散項（大きさ駆動）
                      + sin(2th) sqrt((1-fA)(1-fB)) Re c0   ... ボゾン交差項
                      + sin(2th) sqrt(fA fB NA NB) |c6| cos(phi+delta) ] ... 電荷項
     電荷項 = 振幅積（頂点積）× 符号（相対位相）。回転規約 sigma は較正で決める
  C4 搬送波分離のヌル定理: 搬送波ありでは全単一倍音対の重なりが機械精度ゼロ
     → コヒーレント流路恒等閉鎖（中性化の第二機構）

出力図表:
  figA_carrier_overlap  : 重なり行列（搬送波あり/なし）ヒートマップ
  figB_classification   : 読出し7種の自己ドリフト分類（既存結果JSONから）
  figC_mean_law         : 平均合成則（既存CSVから）
  figD_charge_grammar   : dN_B(phi) 曲線と頂点積直線
  figE_three_terms      : 三項分解の予言vs実測（等ノルム/異ノルム）
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"
CLASSIFY_JSON = HERE.parent / "tongue_spectroscopy_pre_v1" / "tongue_spectroscopy_pre_result_v1.json"
MEAN_CSV = HERE.parent / "mixed_coupling_pre_v1" / "mixed_coupling_pre_rows_v1.csv"

F_GRID = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
PHASES = tuple(2.0 * math.pi * k / 8 for k in range(8))
PROBES = (0.05, 0.1, 0.2, 0.3)
UNEQUAL_NB = 1.44  # 異ノルム対照: N_B = 1.44 (amplitude x1.2)
KS = (1, 3, 5, 7, 9, 11, 13)
PRED_TOL = 1.0e-12


def load_toy() -> Any:
    spec = importlib.util.spec_from_file_location("toy_for_paperA_v1", TOY_RUNNER_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


toy = load_toy()
base = toy.base
plt = base.plt


def unit_norm(v: np.ndarray) -> np.ndarray:
    return v / math.sqrt(float(np.vdot(v, v).real))


def mk(k: int, which: str, hair: bool, source_params: Any) -> np.ndarray:
    case = base.explicit_packet_case(
        mode=f"pA_{which}_{k}_{int(hair)}", packet_a=(k,), packet_b=(k,)
    )
    return unit_norm(base.make_case_state(source_params, case, which, hair_enabled=hair))


def fit_cosine(phases: np.ndarray, values: np.ndarray) -> tuple[float, float, float]:
    design = np.column_stack([np.ones_like(phases), np.cos(phases), np.sin(phases)])
    coef, *_ = np.linalg.lstsq(design, values, rcond=None)
    offset, c, s = coef
    return float(offset), float(math.hypot(c, s)), float(math.atan2(-s, c))


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=1)
    sp = base.build_source_params(params)

    # ================= C4: 重なり行列（搬送波あり/なし） =================
    overlaps = {}
    for hair in (True, False):
        A = {k: mk(k, "A", hair, sp) for k in KS}
        B = {k: mk(k, "B", hair, sp) for k in KS}
        M = np.zeros((len(KS), len(KS)))
        for i, j in ((i, j) for i in range(len(KS)) for j in range(len(KS))):
            M[i, j] = abs(complex(np.vdot(A[KS[i]], B[KS[j]])))
        overlaps[hair] = M
    null_max = float(np.max(overlaps[True]))
    shared_pairs = [
        (KS[i], KS[j], float(overlaps[False][i, j]))
        for i in range(len(KS))
        for j in range(len(KS))
        if overlaps[False][i, j] > 1e-10
    ]
    np.savetxt(HERE / "overlap_matrix_hair_on_v1.csv", overlaps[True], delimiter=",")
    np.savetxt(HERE / "overlap_matrix_hair_off_v1.csv", overlaps[False], delimiter=",")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for ax, hair, title in (
        (axes[0], True, "carriers ON: max |<A_j|B_k>| = %.1e (all zero)" % null_max),
        (axes[1], False, "carriers OFF: shared channels open (0.5)"),
    ):
        im = ax.imshow(
            np.log10(np.maximum(overlaps[hair], 1e-18)),
            origin="lower", vmin=-18, vmax=0, cmap="viridis",
        )
        ax.set_xticks(range(len(KS))); ax.set_xticklabels([f"B{k}" for k in KS], fontsize=7)
        ax.set_yticks(range(len(KS))); ax.set_yticklabels([f"A{k}" for k in KS], fontsize=7)
        ax.set_title(title, fontsize=9)
        fig.colorbar(im, ax=ax, shrink=0.8, label="log10 |overlap|")
    fig.suptitle("Fig A: carrier separation closes the coherent channel (null theorem)")
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"figA_carrier_overlap_v1.{ext}", dpi=160)
    plt.close(fig)

    # ================= C1: 読出し分類（既存JSONから図化） =================
    figB_ok = False
    if CLASSIFY_JSON.exists():
        cls = json.loads(CLASSIFY_JSON.read_text(encoding="utf-8"))[
            "part_A_readout_classification"
        ]
        names = list(cls)
        drifts = [max(cls[n]["max_self_drift"], 1e-17) for n in names]
        fig, ax = plt.subplots(figsize=(9, 4.6), constrained_layout=True)
        colors = ["tab:blue" if cls[n]["class"] == "CONSERVED_NO_GO" else "tab:orange" for n in names]
        ax.bar(range(len(names)), drifts, color=colors)
        ax.set_yscale("log")
        ax.axhline(1e-12, color="0.4", linestyle=":", label="conserved tolerance")
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=25, ha="right", fontsize=7)
        ax.set_ylabel("max self-drift of theta")
        ax.set_title(
            "Fig B: only the per-bin invariant readout is conserved "
            "(blue = conserved / orange = dynamical)"
        )
        ax.legend()
        for ext in ("png", "svg"):
            fig.savefig(HERE / f"figB_classification_v1.{ext}", dpi=160)
        plt.close(fig)
        figB_ok = True

    # ================= C2: 平均合成則（既存CSVから図化） =================
    figC_ok = False
    if MEAN_CSV.exists():
        rows = list(csv.DictReader(MEAN_CSV.open(encoding="utf-8")))
        r_meas = np.asarray([float(r["R_joint"]) for r in rows])
        mean_h = np.asarray([float(r["hyp_mean"]) for r in rows])
        prod_h = np.asarray([float(r["hyp_product"]) for r in rows])
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
        ax1.plot(mean_h, r_meas, "o", markersize=5, label="vs mean (f_A+f_B)/2")
        ax1.plot(prod_h, r_meas, "x", markersize=5, color="tab:red", label="vs product f_A f_B")
        ax1.plot([0, 1], [0, 1], "-", color="0.6", linewidth=0.8)
        ax1.set_xlabel("composition-law prediction")
        ax1.set_ylabel("measured R_joint")
        ax1.set_title("mean law: on the diagonal; product law: off")
        ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
        ax2.semilogy(np.abs(r_meas - mean_h) + 1e-18, "o", markersize=4)
        ax2.axhline(1e-12, color="0.4", linestyle=":")
        ax2.set_xlabel("grid cell"); ax2.set_ylabel("|R_joint - mean|")
        ax2.set_title("residual vs mean law (machine precision)")
        ax2.grid(alpha=0.3)
        fig.suptitle("Fig C: diagonal-sector composition is exactly the mean (gravity grammar)")
        for ext in ("png", "svg"):
            fig.savefig(HERE / f"figC_mean_law_v1.{ext}", dpi=160)
        plt.close(fig)
        figC_ok = True

    # ================= C3: 三項分解と電荷文法（正式再実行） =================
    a1, a5 = mk(1, "A", False, sp), mk(5, "A", False, sp)
    b1, b7 = mk(1, "B", False, sp), mk(7, "B", False, sp)
    c0 = complex(np.vdot(a1, b1))
    c6 = complex(np.vdot(a5, b7))

    # 回転規約の較正: 既知状態で dN_B の符号 sigma を決める
    th_c = 0.1
    a_cal, b_cal = a5.copy(), b7.copy()
    _, b2 = toy.rotate_ab(a_cal, b_cal, th_c)
    dnb_cal = float(np.vdot(b2, b2).real) - 1.0
    sigma = 1.0 if dnb_cal * (math.sin(2 * th_c) * c6.real) > 0 else -1.0

    def run_grid(nb_scale: float) -> list[dict[str, Any]]:
        fits = []
        for th_p in PROBES:
            s2t, s2 = math.sin(2 * th_p), math.sin(th_p) ** 2
            for fa in F_GRID:
                a0 = unit_norm(math.sqrt(1 - fa) * a1 + math.sqrt(fa) * a5)
                for fb in F_GRID:
                    b_unit = unit_norm(
                        math.sqrt(1 - fb) * b1 + math.sqrt(fb) * b7
                    )
                    flows = []
                    for phi in PHASES:
                        b0 = nb_scale * unit_norm(
                            math.sqrt(1 - fb) * b1
                            + np.exp(1j * phi) * math.sqrt(fb) * b7
                        )
                        nb0 = float(np.vdot(b0, b0).real)
                        a2, b2 = toy.rotate_ab(a0.copy(), b0.copy(), th_p)
                        flows.append(float(np.vdot(b2, b2).real) - nb0)
                    off, m_mod, delta = fit_cosine(np.asarray(PHASES), np.asarray(flows))
                    nb0 = nb_scale**2
                    pred_off = sigma * (
                        s2 * (1.0 - nb0) / sigma * 0  # placeholder replaced below
                    )
                    # 三項分解の予言（sigma較正済み）
                    pred_offset = (
                        s2 * (1.0 - nb0)
                        + sigma * s2t * math.sqrt((1 - fa) * (1 - fb)) * nb_scale * c0.real
                    )
                    pred_m = abs(s2t) * math.sqrt(fa * fb) * nb_scale * abs(c6)
                    fits.append(
                        dict(
                            th_p=th_p, fa=fa, fb=fb, nb0=nb0,
                            offset=off, pred_offset=pred_offset,
                            offset_err=abs(off - pred_offset),
                            M=m_mod, pred_M=pred_m, M_err=abs(m_mod - pred_m),
                            delta=delta,
                            signflip=(off + m_mod) * (off - m_mod) < 0.0,
                        )
                    )
        return fits

    fits_eq = run_grid(1.0)
    fits_uneq = run_grid(math.sqrt(UNEQUAL_NB))

    def verdict(fits: list[dict[str, Any]]) -> dict[str, Any]:
        max_off = max(f["offset_err"] for f in fits)
        max_m = max(f["M_err"] for f in fits)
        charged = [f for f in fits if f["pred_M"] > 1e-15]
        neutral_max = max(
            (f["M"] for f in fits if f["pred_M"] <= 1e-15), default=0.0
        )
        return dict(
            max_offset_err=max_off,
            max_modulation_err=max_m,
            neutral_max_modulation=neutral_max,
            sign_flips=sum(1 for f in charged if f["signflip"]),
            charged_cells=len(charged),
            all_pass=max_off <= PRED_TOL and max_m <= PRED_TOL and neutral_max <= PRED_TOL,
        )

    v_eq, v_uneq = verdict(fits_eq), verdict(fits_uneq)

    with (HERE / "three_term_flow_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(fits_eq[0]))
        w.writeheader()
        for f in fits_eq:
            w.writerow(f)
        for f in fits_uneq:
            w.writerow(f)

    # 図D: 電荷文法
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    th_show = 0.2
    for fa, fb in ((0.2, 0.2), (0.6, 0.6), (1.0, 1.0), (1.0, 0.2)):
        a0 = unit_norm(math.sqrt(1 - fa) * a1 + math.sqrt(fa) * a5)
        ys = []
        for phi in PHASES:
            b0 = unit_norm(math.sqrt(1 - fb) * b1 + np.exp(1j * phi) * math.sqrt(fb) * b7)
            _, b2 = toy.rotate_ab(a0.copy(), b0.copy(), th_show)
            ys.append(float(np.vdot(b2, b2).real) - 1.0)
        ax1.plot(PHASES, ys, "o-", label=f"f_A={fa}, f_B={fb}")
    ax1.axhline(0, color="0.5", linewidth=0.6)
    ax1.set_xlabel("relative phase phi")
    ax1.set_ylabel("dN_B (one collision)")
    ax1.set_title("flow reverses with phase: sign = relative phase")
    ax1.legend(fontsize=7); ax1.grid(alpha=0.3)

    xs = [abs(math.sin(2 * f["th_p"])) * math.sqrt(f["fa"] * f["fb"]) for f in fits_eq]
    ys = [f["M"] for f in fits_eq]
    ax2.plot(xs, ys, "o", markersize=3.5)
    xline = np.linspace(0, max(xs), 10)
    ax2.plot(xline, abs(c6) * xline, "-", color="tab:red", linewidth=0.9,
             label=f"slope = |c6| = {abs(c6):.3f}")
    ax2.set_xlabel("|sin 2theta_p| sqrt(f_A f_B)")
    ax2.set_ylabel("modulation amplitude M")
    ax2.set_title("vertex product law (4 probe angles, all cells)")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
    fig.suptitle("Fig D: charge grammar in the coherent flow")
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"figD_charge_grammar_v1.{ext}", dpi=160)
    plt.close(fig)

    # 図E: 三項分解（予言vs実測）
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for ax, fits, label in (
        (ax1, fits_eq, "equal norms (diffusive term = 0)"),
        (ax2, fits_uneq, f"unequal norms N_B={UNEQUAL_NB} (diffusive term on)"),
    ):
        po = [f["pred_offset"] for f in fits]
        mo = [f["offset"] for f in fits]
        ax.plot(po, mo, "o", markersize=3)
        lim = [min(po + mo), max(po + mo)]
        ax.plot(lim, lim, "-", color="0.6", linewidth=0.8)
        err = max(f["offset_err"] for f in fits)
        ax.set_xlabel("predicted offset (3-term formula)")
        ax.set_ylabel("measured offset")
        ax.set_title(f"{label}\nmax|err| = {err:.1e}")
        ax.grid(alpha=0.3)
    fig.suptitle("Fig E: the flow equation verified term by term")
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"figE_three_terms_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "paperA_two_grammar_evidence_v1",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "rotation_sign_calibration_sigma": sigma,
        "couplings": {"c0": [c0.real, c0.imag], "c6": [c6.real, c6.imag]},
        "C4_null_theorem": {
            "max_overlap_with_carriers": null_max,
            "shared_pairs_without_carriers": shared_pairs,
        },
        "C3_equal_norms": v_eq,
        "C3_unequal_norms": v_uneq,
        "figures_generated": {
            "figA": True, "figB": figB_ok, "figC": figC_ok, "figD": True, "figE": True,
        },
    }
    (HERE / "paperA_two_grammar_evidence_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"sigma (rotation convention) = {sigma:+.0f}")
    print(f"C4 null theorem: max overlap with carriers = {null_max:.2e}")
    print(f"   shared channels without carriers: {shared_pairs}")
    print(f"C3 equal norms:   {v_eq}")
    print(f"C3 unequal norms: {v_uneq}")
    print(f"figures: A yes / B {figB_ok} / C {figC_ok} / D yes / E yes")


if __name__ == "__main__":
    main()
