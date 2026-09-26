#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate revised Paper 5 figures and tables from archived audit data only.

No physics simulation is performed. The program reads stored CSV/JSON audit outputs,
creates the revised 2->5->6 state-closure figures/tables, and writes SHA256 manifests.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import font_manager, rcParams


def configure_matplotlib() -> None:
    names = {f.name for f in font_manager.fontManager.ttflist}
    if "Noto Sans CJK JP" in names:
        rcParams["font.family"] = "Noto Sans CJK JP"
    rcParams["axes.unicode_minus"] = False
    rcParams["svg.fonttype"] = "path"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save_figure(fig: plt.Figure, svg: Path, png_preview: Path | None = None) -> None:
    fig.savefig(svg, bbox_inches="tight")
    if png_preview is not None:
        fig.savefig(png_preview, dpi=180, bbox_inches="tight")
    plt.close(fig)


def box(ax, x, y, w, h, title, body, linestyle="-"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015", fill=False, linewidth=1.5, linestyle=linestyle)
    ax.add_patch(p)
    ax.text(x + w/2, y + h*0.68, title, ha="center", va="center", fontsize=12, fontweight="bold")
    ax.text(x + w/2, y + h*0.30, body, ha="center", va="center", fontsize=10)


def arrow(ax, x1, y1, x2, y2, text=None):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="->", mutation_scale=14, linewidth=1.4)
    ax.add_patch(a)
    if text:
        ax.text((x1+x2)/2, (y1+y2)/2 + 0.025, text, ha="center", va="bottom", fontsize=9)


def figure01_state_audit(out: Path, preview: Path):
    fig, ax = plt.subplots(figsize=(11, 5.7))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_title("第五思考実験・改訂版：状態監査の経路（2 → 5 → 6）", fontsize=15, pad=16)

    box(ax, 0.04, 0.56, 0.24, 0.25,
        "Paper 4 の見かけの状態",
        "$X=(a,b)$\n二成分の派生表現")
    box(ax, 0.38, 0.56, 0.26, 0.25,
        "第1監査：5永続状態",
        "$(\\chi,p,e,\\phi,t)$\n$\\Downarrow$\n$Z_5=(U,P,E,H,Q)$")
    box(ax, 0.73, 0.56, 0.23, 0.25,
        "追加監査：第6状態",
        "外部 $\\nu$ を発見\n$\\mathcal{N}=\\nu$ を内部化")
    arrow(ax, 0.28, 0.685, 0.38, 0.685, "上流情報を追跡")
    arrow(ax, 0.64, 0.685, 0.73, 0.685, "外部引数を監査")

    box(ax, 0.27, 0.18, 0.46, 0.24,
        "改訂後の永続論理状態",
        r"$Z_6=(U,P,E,H,Q,\mathcal{N})$" + "\n" +
        "6は現在表現の永続状態数。数学的最小性は未証明。")
    arrow(ax, 0.845, 0.56, 0.70, 0.42)
    ax.text(0.05, 0.025,
            "注：RK4 work state と 11相 one-hot 位相 q は明示 machine state として別に存在する。\n"
            "したがって『全実スカラー自由度が6』という主張ではない。",
            fontsize=9.5, va="bottom")
    save_figure(fig, out / "figure01_state_audit_2_to_5_to_6_v2.svg", preview / "figure01.png")


def figure02_interaction_scope(out: Path, preview: Path):
    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_title("改訂版の状態閉包と同期更新：本論文で実証する範囲", fontsize=15, pad=16)

    box(ax, 0.04, 0.58, 0.28, 0.24, "旧 full state を一度に読む",
        r"永続状態 $Z_6$" + "\n" + r"+ work state + one-hot $q$")
    box(ax, 0.39, 0.58, 0.25, 0.24, "同一 microstep 写像",
        r"$\Psi_{m+1}=\sum_j q_{j,m}F_j(\Psi_m)$" + "\n" + "更新後値を同一stepで再利用しない")
    box(ax, 0.71, 0.58, 0.25, 0.24, "次 full state",
        r"$Z_6'=(U',P',E',H',Q',\mathcal{N}')$" + "\n" + r"$\mathcal{N}'=\mathcal{N}$ も同じ写像内")
    arrow(ax, 0.32, 0.70, 0.39, 0.70)
    arrow(ax, 0.64, 0.70, 0.71, 0.70)

    box(ax, 0.13, 0.18, 0.27, 0.20, "読み出し（非帰還）",
        r"$x,y,r,t$" + "\n" + r"$t=\tau_0\log|Q|$" + "\n状態更新へ戻さない")
    arrow(ax, 0.24, 0.58, 0.27, 0.38)

    box(ax, 0.56, 0.16, 0.36, 0.24, "主張しないこと", 
        "6状態の数学的最小性\n11相機械の自然界での必然性\n既存B方式全体の『完全乗法内部演算』適合",
        linestyle="--")
    ax.text(0.04, 0.06,
            "この図は状態閉包・同期更新・読み出し非帰還・ν内部化を示す。\n"
            "後から定めた最厳格R4の完全監査までを合格済みとする図ではない。",
            fontsize=9.5)
    save_figure(fig, out / "figure02_state_closure_sync_scope_v2.svg", preview / "figure02.png")


def figure03_methodB_heatmap(src: Path, out: Path, preview: Path):
    df = pd.read_csv(src / "figure05_methodB_strict_cases.csv")
    states = ["U","P","E","H","Q"]
    data = df[states].to_numpy(float)
    floor = 1e-18
    z = np.log10(np.maximum(np.abs(data), floor))
    labels = [f"p={r.p0:g}, e={r.e0:g}, q={r.q:g}" for r in df.itertuples()]
    fig, ax = plt.subplots(figsize=(9.8, 7.0))
    im = ax.imshow(z, aspect="auto")
    ax.set_xticks(np.arange(len(states)), states)
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.set_title("11相 Method B：12条件の状態差（log10 |Δ|、0は表示下限）")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10 |Δ|")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            txt = "0" if data[i,j] == 0 else f"{data[i,j]:.1e}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7)
    save_figure(fig, out / "figure03_methodB_strict_12case_heatmap_v2.svg", preview / "figure03.png")


def figure04_angular(src: Path, out: Path, preview: Path):
    df = pd.read_csv(src / "figure06_angular_order_convergence_source.csv")
    fig, ax = plt.subplots(figsize=(8.7, 5.7))
    for (p,e), g in df.groupby(["p0","e0"], sort=True):
        g = g.sort_values("N")
        y = np.maximum(g["phi_abs_error"].to_numpy(float), 1e-16)
        ax.plot(g["N"], y, marker="o", label=f"p={p:g}, e={e:g}")
    ax.set_yscale("log")
    ax.set_xlabel("角度級数打切り次数 N")
    ax.set_ylabel("3軌道後 |Δφ|")
    ax.set_title("角度級数 N=12,16,20,24 の収束（9条件）")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncol=3)
    save_figure(fig, out / "figure04_angular_order_convergence_v2.svg", preview / "figure04.png")


def figure05_integration(src: Path, out: Path, preview: Path):
    df = pd.read_csv(src / "figure07_integration_order_convergence_source.csv")
    # Baseline representative case used in the paper text.
    d = df[(df.p0==24.3) & (df.e0==0.45)].copy()
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    for order, g in d.groupby("order", sort=True):
        g = g.sort_values("steps_per_orbit")
        ax.loglog(g["steps_per_orbit"], g["combined"], marker="o", label=f"設計次数 {int(order)}")
    ax.set_xlabel("steps / orbit")
    ax.set_ylabel("combined error")
    ax.set_title("同じ永続状態での4・6・8次収束（p=24.3, e=0.45）")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    save_figure(fig, out / "figure05_integration_order_convergence_v2.svg", preview / "figure05.png")


def figure06_nine_case(src: Path, out: Path, preview: Path):
    df = pd.read_csv(src / "figure08_integration_order_9case_finest.csv")
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    for col in ["4","6","8"]:
        ax.plot(x, df[col], marker="o", label=f"設計次数 {col}")
    ax.set_yscale("log")
    ax.set_xticks(x, df["case"], rotation=35, ha="right")
    ax.set_ylabel("64 steps/orbit での combined error")
    ax.set_title("9条件における高次化の移植性（最細分解能）")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    save_figure(fig, out / "figure06_integration_order_9case_summary_v2.svg", preview / "figure06.png")


def write_tables(src: Path, out: Path):
    # Table 1: revised audit inventory.
    rows1 = [
        ["X=(a,b)", "Paper 4で表示された二成分", "派生表現", "r,phiから構成され、単独では次状態を閉じない", "Paper 4→第1監査"],
        ["U=e^{i chi}", "処理位相", "永続論理状態", "chiの周期情報を保持", "5状態監査"],
        ["P=p", "軌道要素", "永続論理状態", "放射反作用で更新", "5状態監査"],
        ["E=e", "軌道要素", "永続論理状態", "放射反作用で更新", "5状態監査"],
        ["H=e^{i phi}", "軌道方位位相", "永続論理状態", "角度進行を乗法的に保持", "5状態監査"],
        ["Q=e^{t/tau0}", "時計状態", "永続論理状態", "tは読み出し時のみlogで取得", "5状態監査"],
        ["N=nu", "対称質量比", "永続論理状態", "旧実装では外部引数。追加監査で状態漏れと判定し内部化", "追加nu監査"],
        ["q[0..10]", "11相 one-hot 処理位相", "machine state", "外部loop counterの代替。macrostep後に巡回", "11相同期監査"],
        ["RK work registers", "k1..k4等", "work/machine state", "更新途中値を暗黙保持しないため明示", "11相同期監査"],
    ]
    with (out/"table01_state_audit_summary_v2.csv").open("w", newline="", encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["量","役割","分類","監査理由","監査段階"]); w.writerows(rows1)

    rows2 = [
        ["散逸物理モデル", "leading 2.5PN radiation reaction", "高次PN項は未導入", "物理モデル近似次数", "6永続状態とは別問題"],
        ["角度級数", "N=12 baseline; tested 12,16,20,24", "O(u^(N+1))", "級数打切り", "状態数を増やさず高次化"],
        ["Paper 4時間積分", "RK4", "local O(h^5), global O(h^4)", "数値積分次数", "11相展開はこの処理順の明示化"],
        ["B 11相展開", "RK4の代数的再配置", "新たな打切り誤差なし", "machine-state深度", "物理的11段階性の主張ではない"],
        ["Taylor次数検証", "4,6,8", "実測収束を別途確認", "数値近似深度", "永続状態数を増やさず高次化"],
        ["nu内部化", "N=nu identity state", "近似ではない", "状態閉包の修正", "5→6は精度向上ではなく情報漏れ除去"],
    ]
    with (out/"table02_approximation_and_state_scope_v2.csv").open("w", newline="", encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["層","現在/検証次数","誤差・未保持項","意味","状態数との関係"]); w.writerows(rows2)

    nu = json.load((src/"nu_internal_state_control_summary.json").open(encoding="utf-8"))
    saved = json.load((src/"saved_C1_check.json").open(encoding="utf-8"))
    angular = pd.read_csv(src/"figure06_angular_order_convergence_source.csv")
    hard = angular[(angular.p0==24.3) & (angular.e0==0.55)].sort_values("N")
    strict_global = pd.read_csv(src/"figure05b_methodB_strict_global.csv")
    strict_map = dict(zip(strict_global["state"], strict_global["global_max_abs"]))
    rows3 = [
        ["strict Method B", "12条件 global max |ΔH|", strict_map.get("H", math.nan), "11相状態機械 vs reference"],
        ["strict Method B", "12条件 global max |ΔQ|", strict_map.get("Q", math.nan), "11相状態機械 vs reference"],
        ["nu内部化", "12条件×4000 old5/new6 最大差", max(nu["global_max_abs_old5_vs_new6_first5"].values()), "保存された全比較量で0"],
        ["nu内部化", "12条件 raw core bitwise mismatch", nu["global_raw_bitwise_core_mismatch_count"], "0"],
        ["nu内部化", "N bitwise mismatch", nu["global_raw_N_bitwise_mismatch_count"], "0"],
        ["nu内部化", "3条件×12000 N micro hold failure", nu["long_run_total_N_micro_hold_fail_count"], "0"],
        ["nu内部化", "3条件×12000 phase failure", nu["long_run_total_phase_fail_count"], "0"],
        ["nu内部化", "saved C1 old/new error profile identical", saved.get("error_dicts_exactly_equal"), "旧B系の既存微小誤差まで同一"],
    ]
    for r in hard.itertuples():
        rows3.append(["角度次数", f"p=24.3,e=0.55,N={int(r.N)} の |Δphi|", r.phi_abs_error, "3軌道後 exact angular factor との比較"])
    with (out/"table03_key_numerical_results_v2.csv").open("w", newline="", encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["実験","指標","値","解釈"]); w.writerows(rows3)

    # One markdown file for direct insertion/checking.
    md = []
    for fn, title in [
        ("table01_state_audit_summary_v2.csv", "表1 状態監査"),
        ("table02_approximation_and_state_scope_v2.csv", "表2 近似次数と状態閉包の区別"),
        ("table03_key_numerical_results_v2.csv", "表3 主要数値結果"),
    ]:
        df = pd.read_csv(out/fn)
        md.append(f"## {title}\n\n" + df.to_markdown(index=False) + "\n")
    (out/"tables_for_paper_v2_ja.md").write_text("\n".join(md), encoding="utf-8")


def write_indices(root: Path, src: Path, figs: Path, tables: Path):
    fig_index = """# 第五思考実験・改訂版 Figure Index v2\n\n| Figure | File | Role | Source |\n|---|---|---|---|\n| 1 | figure01_state_audit_2_to_5_to_6_v2.svg | 2→5→6の状態監査経路 | Paper 4監査 + nu内部化監査 |\n| 2 | figure02_state_closure_sync_scope_v2.svg | 状態閉包・同期更新・主張範囲 | strict audit rules / nu derivation |\n| 3 | figure03_methodB_strict_12case_heatmap_v2.svg | 11相Method Bの12条件照合 | figure05_methodB_strict_cases.csv |\n| 4 | figure04_angular_order_convergence_v2.svg | 角度級数Nの収束 | figure06_angular_order_convergence_source.csv |\n| 5 | figure05_integration_order_convergence_v2.svg | 4/6/8次数値収束 | figure07_integration_order_convergence_source.csv |\n| 6 | figure06_integration_order_9case_summary_v2.svg | 9条件への高次化移植性 | figure08_integration_order_9case_finest.csv |\n"""
    (root/"FIGURE_INDEX_v2_ja.md").write_text(fig_index, encoding="utf-8")
    exp_index = [
        ["00","00_reference_from_paper4","Paper4 baseline and C1 raw data","reference"],
        ["01","01_五状態同期写像の再構成","identify persistent five-state map","verified"],
        ["02","02_乗法状態表現と同期インライン化の経緯","multiplicative encoding history","historical"],
        ["03","03_A局所exp-log_vs_B内部状態展開","compare A/B","compared"],
        ["04","04_B方式_11相同期状態機械_厳密監査","old-state-only 11-phase state machine","verified"],
        ["05","05_演算ブロック完全展開_Laurent_N12","explicit algebra blocks","verified"],
        ["06","06_近似次数と状態_S不増加検証","raise approximation order without persistent-state increase","verified"],
        ["08","08_ν第6状態内部化_対照実験_20260925","external nu -> internal identity state N","verified"],
    ]
    with (root/"EXPERIMENT_INDEX_v2.csv").open("w", newline="", encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["experiment","folder","purpose","status"]); w.writerows(exp_index)


def write_manifest(root: Path, provenance: dict[str, str]):
    with (root/"SOURCE_PROVENANCE_v2.json").open("w", encoding="utf-8") as f:
        json.dump(provenance, f, ensure_ascii=False, indent=2)
    files = [p for p in root.rglob("*") if p.is_file() and p.name != "SHA256SUMS_v2.txt"]
    lines=[]
    for p in sorted(files):
        lines.append(f"{sha256(p)}  {p.relative_to(root)}")
    (root/"SHA256SUMS_v2.txt").write_text("\n".join(lines)+"\n", encoding="utf-8")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args=ap.parse_args()
    root=args.root.resolve(); src=root/"source_data"; figs=root/"figures"; tables=root/"tables"; preview=root/"previews"
    for d in [src,figs,tables,preview]: d.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()
    figure01_state_audit(figs, preview)
    figure02_interaction_scope(figs, preview)
    figure03_methodB_heatmap(src, figs, preview)
    figure04_angular(src, figs, preview)
    figure05_integration(src, figs, preview)
    figure06_nine_case(src, figs, preview)
    write_tables(src, tables)
    write_indices(root, src, figs, tables)
    provenance = {
        "figure05_methodB_strict_cases.csv": "Google Drive Paper5/plot_data; Drive file id 1KatO5Xzs9vP-eI8OAqSxG9dkMEiE0-qH",
        "figure05b_methodB_strict_global.csv": "Google Drive Paper5/plot_data; Drive file id 1XUDCbNpN6U0WoJOVCXxpV5qUBiOZTRgJ",
        "figure06_angular_order_convergence_source.csv": "Google Drive Paper5/plot_data; Drive file id 1PyLx92F1rn84fkvYT4A2Yv2qjqfi7BP5",
        "figure07_integration_order_convergence_source.csv": "Google Drive Paper5/plot_data; Drive file id 1XnJPUzCQnau41QTVcWKdPmibnnTN898G",
        "figure08_integration_order_9case_finest.csv": "Google Drive Paper5/plot_data; Drive file id 1TTmxEPXPGQ8By8dPPw5UnUKUY7HZVEEN",
        "nu_internal_state_control_summary.json": "Google Drive Paper5/08_nu; Drive file id 1MutSKARNhu0LS1jidsHWLMOdr6xCUplO",
        "nu_internal_state_control_case_summary.csv": "Google Drive Paper5/08_nu; Drive file id 1STIF4pA0RCIYP3T28Fg2IF_Lxo2UYs00",
        "long_run_checks.json": "Google Drive Paper5/08_nu; Drive file id 16sxoNd1Q37pnZ9rIyQXK6Q9rtTKW5D4a",
        "saved_C1_check.json": "Google Drive Paper5/08_nu; Drive file id 10EdAPw4YIAvBqrL66jzgZnzRdftn6Zd4"
    }
    write_manifest(root, provenance)
    print(f"Generated revised Paper 5 figures/tables at {root}")

if __name__ == "__main__":
    main()
