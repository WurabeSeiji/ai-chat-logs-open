#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""基礎論文（体ゲージ）用の新規概念説明図 F1〜F6 を生成する。

出所: 運動量レジスタの検討.md §18.2 の図表計画（2026-08-19）。
既存実験事実の図は先行論文から引用する（§18.1）。本スクリプトが作るのは
新規概念の説明図のみで、実験データは一切使わない（決定的・乱数は固定シード）。

出力: 本フォルダに fig_f1〜fig_f6 の .png と .svg（各図2形式）。

  F1 体ゲージの構成要素（並進・回転・間隔を一つの枠が担う）      … §13.1
  F2 絶対中心と主軸枠（材料は配置自身のみ）                      … §13.3b/c
  F3 正準極分解（動径=観測不能側 / 接線=可視側）                 … §13.3b
  F4 三種の観測時計（同一位相データへの三種の復調）              … §13.6・§16
  F5 対スカラー族→ベクトル（距離幾何による復元）                … §5
  F6 二層アーキテクチャ（F=状態のみ / D=体ゲージ / G=対読出し）  … §3・§14.2

使い方: python3 make_foundation_figures_v1.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent

plt.rcParams["font.family"] = ["Hiragino Sans", "Hiragino Kaku Gothic ProN",
                               "Yu Gothic", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

RNG_SEED = 42  # 決定的。全図共通


def _save(fig, stem: str) -> None:
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  {stem}.png / .svg")


# ----------------------------------------------------------------------
def fig_f1_body_gauge() -> None:
    """F1: 体ゲージ Φ_A = (θ_A, O_A, (Ĝ_A, ρ_A)) の構成要素。"""
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.set_aspect("equal")

    # 広域枠（絶対中心と主軸）
    ax.plot(0, 0, marker="*", ms=18, color="black", zorder=5)
    ax.annotate("絶対中心（広域枠の原点）", (0, 0), (0.15, -0.55), fontsize=10)
    for dx, dy, lab in [(1.0, 0.0, "X（第1主軸）"), (0.0, 1.0, "Y（第2主軸）")]:
        ax.add_patch(FancyArrowPatch((0, 0), (dx * 1.5, dy * 1.5),
                                     arrowstyle="-|>", mutation_scale=16,
                                     color="black", lw=1.4))
        ax.annotate(lab, (dx * 1.55, dy * 1.55), fontsize=10)

    # 体 A の位置（並進オフセット θ_A）
    ax_pos = np.array([2.6, 1.7])
    ax.add_patch(FancyArrowPatch((0, 0), tuple(ax_pos), arrowstyle="-|>",
                                 mutation_scale=14, color="tab:blue",
                                 lw=1.6, linestyle="--"))
    ax.annotate(r"並進オフセット $\theta_A$（率＝速度）",
                (0.55, 0.05), color="tab:blue", fontsize=10, rotation=0)

    # 体 A の局所枠（回転 O_A）
    ang = np.deg2rad(28.0)
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    for v, lab in [(np.array([0.9, 0.0]), "x'"), (np.array([0.0, 0.9]), "y'")]:
        w = rot @ v
        ax.add_patch(FancyArrowPatch(tuple(ax_pos), tuple(ax_pos + w),
                                     arrowstyle="-|>", mutation_scale=13,
                                     color="tab:red", lw=1.5))
        ax.annotate(lab, tuple(ax_pos + w * 1.12), color="tab:red", fontsize=10)
    ax.annotate(r"回転 $O_A$（率＝角速度）", (ax_pos[0] + 0.15, ax_pos[1] - 0.62),
                color="tab:red", fontsize=10)

    # 間隔（局所計量 Ĝ_A・スケール ρ_A）＝楕円
    ell = Ellipse(tuple(ax_pos), width=2.3, height=1.5,
                  angle=np.rad2deg(ang), fill=False, color="tab:green",
                  lw=1.8, linestyle="-")
    ax.add_patch(ell)
    ax.annotate(r"間隔 $(\widehat G_A,\ \rho_A)$（歪み率＝加速度）",
                (ax_pos[0] - 0.3, ax_pos[1] + 1.0), color="tab:green", fontsize=10)

    ax.plot(*ax_pos, marker="o", ms=9, color="tab:blue", zorder=6)
    ax.annotate("体 A", (ax_pos[0] + 0.1, ax_pos[1] + 0.12), fontsize=11)

    ax.set_title("F1  体ゲージ $\\Phi_A$ ——並進・回転・間隔を一つの枠が担う（Cartan 動標構型）",
                 fontsize=11)
    ax.set_xlim(-1.2, 4.6)
    ax.set_ylim(-1.3, 3.5)
    ax.axis("off")
    _save(fig, "fig_f1_body_gauge_v1")


# ----------------------------------------------------------------------
def fig_f2_absolute_frame() -> None:
    """F2: 絶対中心と主軸枠。材料は配置（関係の距離幾何復元）自身のみ。"""
    rng = np.random.default_rng(RNG_SEED)
    # 異方的な配置（N=16）
    n = 16
    base = rng.normal(size=(n, 2)) @ np.array([[1.9, 0.55], [0.0, 0.85]])
    base -= base.mean(axis=0)  # 正準中心化

    t = base.T @ base / n
    lam, vec = np.linalg.eigh(t)
    order = np.argsort(lam)[::-1]
    lam, vec = lam[order], vec[:, order]

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.set_aspect("equal")
    ax.scatter(base[:, 0], base[:, 1], s=45, color="tab:blue", zorder=4,
               label="頂点（体）N=16")

    # 中心
    ax.plot(0, 0, marker="*", ms=20, color="black", zorder=6)
    ax.annotate("絶対中心＝重心\n（二重中心化が関係だけから毎瞬定める）",
                (0.08, -0.15), fontsize=9, va="top")

    # 主軸（長さ順）
    for k, (colr, lab) in enumerate([("tab:red", "X＝第1主軸（最長）"),
                                     ("tab:green", "Y＝第2主軸")]):
        v = vec[:, k] * np.sqrt(lam[k]) * 2.2
        ax.add_patch(FancyArrowPatch((0, 0), tuple(v), arrowstyle="-|>",
                                     mutation_scale=16, color=colr, lw=2.0))
        ax.add_patch(FancyArrowPatch((0, 0), tuple(-v), arrowstyle="-",
                                     color=colr, lw=1.2, linestyle=":"))
        ax.annotate(lab, tuple(v * 1.12), color=colr, fontsize=10)

    ax.annotate("符号（±）と掌性は規約\n縮退点は逐次整列で通過（§13.3c）",
                (0.98, 0.02), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="gray"))

    ax.set_title("F2  絶対ゲージ空間——材料は配置自身のみ（外部参照なし＝無名性を破らない）",
                 fontsize=11)
    ax.legend(loc="upper left", fontsize=9)
    ax.axis("off")
    _save(fig, "fig_f2_absolute_frame_v1")


# ----------------------------------------------------------------------
def fig_f3_polar_decomposition() -> None:
    """F3: 正準極分解。動径=観測不能側(t,R,Q) / 接線=可視 xyz 運動量。"""
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    ax.set_aspect("equal")

    ell = Ellipse((0, 0), width=5.4, height=4.0, fill=False, color="gray", lw=1.6)
    ax.add_patch(ell)
    ax.annotate("読出し面（$C$ 固定の楕円体）", (-2.65, 1.95), fontsize=10,
                color="gray")

    ax.plot(0, 0, marker="*", ms=18, color="black", zorder=5)
    ax.annotate("絶対中心\n（中心方向＝観測不能側）", (-1.75, -0.55), fontsize=9)

    # 体 A（楕円上）
    th = np.deg2rad(35.0)
    pos = np.array([2.7 * np.cos(th), 2.0 * np.sin(th)])
    ax.plot(*pos, marker="o", ms=10, color="tab:blue", zorder=6)
    ax.annotate("体 A", (pos[0] + 0.12, pos[1] + 0.12), fontsize=11)

    # 動径（中心へ）
    ax.add_patch(FancyArrowPatch(tuple(pos), (0.12, 0.09), arrowstyle="-|>",
                                 mutation_scale=15, color="tab:purple", lw=1.8,
                                 linestyle="--"))
    ax.annotate("動径ドリフト $\\dot r_A$（観測不能方向）\n＝ $t, R, Q$ スカラー・スケール $\\rho$\n枠不要の絶対量",
                (0.35, -1.55), color="tab:purple", fontsize=9)

    # 接線
    tang = np.array([-2.7 * np.sin(th) / 1.0, 2.0 * np.cos(th)])
    tang = tang / np.linalg.norm(tang) * 1.5
    ax.add_patch(FancyArrowPatch(tuple(pos), tuple(pos + tang), arrowstyle="-|>",
                                 mutation_scale=15, color="tab:orange", lw=2.0))
    ax.annotate("接線ドリフト（可視3方向）\n＝ xyz 運動量。向きにのみ接続が要る",
                (pos[0] - 0.4, pos[1] + 1.15), color="tab:orange", fontsize=9)

    ax.set_title("F3  正準極分解——間隔は動径（絶対）、方向は接線（ゲージ）", fontsize=11)
    ax.set_xlim(-3.4, 3.9)
    ax.set_ylim(-2.6, 3.4)
    ax.axis("off")
    _save(fig, "fig_f3_polar_decomposition_v1")


# ----------------------------------------------------------------------
def fig_f4_three_clocks() -> None:
    """F4: 三種の観測時計＝同一位相データへの三種の復調。"""
    fig, axes = plt.subplots(3, 1, figsize=(7.6, 7.6))
    n = np.linspace(0, 288, 1000)

    # (1) 集団時計（大域搬送波）
    ax = axes[0]
    ax.plot(n, np.cos(np.pi / 72 * n), color="tab:blue", lw=1.5)
    ax.set_title("集団時計（連続・大域）ω=π/72 ——空間 x,y はこの時計での復調（quadrature）",
                 fontsize=10)
    ax.set_ylabel("位相基準")
    ax.set_xticks([0, 144, 288])
    ax.set_xticklabels(["0", "144（一周）", "288"])

    # (2) 局所時計場（場所ごとの速さ）
    ax = axes[1]
    x = np.linspace(0, 1, 400)
    omega_x = 0.6 + 0.8 * np.exp(-((x - 0.55) ** 2) / 0.02)
    ax.plot(x, omega_x, color="tab:red", lw=1.8)
    ax.fill_between(x, 0.55, omega_x, where=omega_x > 0.62, alpha=0.15,
                    color="tab:red")
    ax.annotate("物質の台（質量 m=⟨ω⟩）", (0.82, 1.22), ha="left",
                color="tab:red", fontsize=9)
    ax.annotate("真空は時を刻まない（G2）", (0.02, 0.72), fontsize=9, color="gray")
    ax.set_title("局所時計場 ω(x)（連続・局所）——質量・寿命・分裂 $t_{split}=π/Δω$",
                 fontsize=10)
    ax.set_ylabel("ω(x)")
    ax.set_xlabel("位置 x")

    # (3) 巡回時計（mod 3）
    ax = axes[2]
    ax.set_aspect("equal")
    circ = Circle((0, 0), 1.0, fill=False, color="tab:green", lw=1.8)
    ax.add_patch(circ)
    for k, lab in enumerate(["0", "+1", "−1"]):
        a = np.pi / 2 - 2 * np.pi * k / 3
        ax.plot([0, np.cos(a)], [0, np.sin(a)], color="tab:green", lw=1.0,
                linestyle=":")
        ax.annotate(lab, (1.22 * np.cos(a), 1.22 * np.sin(a)), ha="center",
                    va="center", fontsize=12, color="tab:green")
    ladder = [2, 4, 8, 16, 32]
    txt = "対生成の倍化梯子:  " + "  ".join(
        f"{m}≡{'+1' if m % 3 == 1 else '−1'}" for m in ladder) + "  (mod 3)"
    ax.annotate(txt, (1.7, 0.45), fontsize=10)
    ax.annotate("3で読む時計だけが全荷電内容を大きさ1に整流（実測100%）\n"
                "3で割り切れない巻きは一価に読めない＝閉じ込め",
                (1.7, -0.35), fontsize=9)
    ax.set_xlim(-1.6, 7.2)
    ax.set_ylim(-1.5, 1.6)
    ax.axis("off")
    ax.set_title("巡回時計「3で読む時計」（離散・$\\mathbb{Z}_3$）——電荷の整流。統計は二重被覆 $\\mathbb{Z}_2$ に住む",
                 fontsize=10)

    fig.suptitle("F4  三種の観測時計＝同一位相データへの三種の復調", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    _save(fig, "fig_f4_three_clocks_v1")


# ----------------------------------------------------------------------
def fig_f5_pair_to_vector() -> None:
    """F5: 対スカラー族 {v_AB} → 距離幾何でベクトル復元（rank≤3, Σp=0）。"""
    rng = np.random.default_rng(RNG_SEED)
    pos = np.array([[0.0, 0.0], [2.4, 0.4], [1.4, 2.0], [-0.6, 1.6]])
    pos -= pos.mean(axis=0)
    vel = rng.normal(size=(4, 2)) * 0.55
    vel -= vel.mean(axis=0)  # Σp = 0（中心枠で恒等）

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 5.2))
    names = ["A", "B", "C", "D"]
    xlim, ylim = (-2.3, 2.6), (-1.9, 1.9)

    # ラベルを線から離すオフセット（対ごとに手で指定・重なり回避）
    label_off = {(0, 1): (0.0, -0.16), (0, 2): (-0.55, 0.10),
                 (0, 3): (-0.30, 0.0), (1, 2): (0.30, 0.0),
                 (1, 3): (0.42, -0.22), (2, 3): (0.0, 0.14)}

    # 左: 対スカラー（枠なしの一次データ）
    ax1.set_aspect("equal")
    ax1.set_xlim(*xlim)
    ax1.set_ylim(*ylim)
    for i in range(4):
        ax1.plot(*pos[i], marker="o", ms=10, color="tab:blue", zorder=5)
        ax1.annotate(names[i], pos[i] + np.array([0.08, 0.10]), fontsize=12)
    for i in range(4):
        for j in range(i + 1, 4):
            mid = (pos[i] + pos[j]) / 2 + np.array(label_off[(i, j)])
            d = pos[j] - pos[i]
            u = d / np.linalg.norm(d)
            v_ab = float(u @ (vel[i] - vel[j]))
            ax1.plot(*np.array([pos[i], pos[j]]).T, color="gray", lw=1.0,
                     linestyle="--")
            ax1.annotate(f"$v_{{{names[i]}{names[j]}}}$={v_ab:+.2f}", mid,
                         fontsize=8.5, color="dimgray", ha="center")
    ax1.set_title("対スカラー族 $\\{v_{AB}\\}$（枠なしで定義できる一次データ）\n"
                  "$v_{AB}=-v_{BA}$（反対称）", fontsize=10)
    ax1.axis("off")

    # 右: 復元されたベクトル（左と同一の座標系・矢印は2倍表示）
    ax2.set_aspect("equal")
    ax2.set_xlim(*xlim)
    ax2.set_ylim(*ylim)
    ax2.plot(0, 0, marker="*", ms=16, color="black", zorder=6)
    for i in range(4):
        ax2.plot(*pos[i], marker="o", ms=10, color="tab:blue", zorder=5)
        ax2.annotate(names[i], pos[i] + np.array([0.08, 0.10]), fontsize=12)
        ax2.add_patch(FancyArrowPatch(tuple(pos[i]), tuple(pos[i] + vel[i] * 1.2),
                                      arrowstyle="-|>", mutation_scale=14,
                                      color="tab:orange", lw=2.0))
    ax2.set_title("ベクトル＝対スカラー族のランク3パターンの名前\n"
                  "（矢印 $\\vec p_A$ は1.2倍表示）", fontsize=10)
    ax2.axis("off")
    fig.text(0.5, 0.015,
             "距離幾何（Schoenberg）で配置と $\\vec p_A$ を復元："
             "埋め込み rank ≤ 3、中心枠で $\\sum\\vec p_A=0$ は恒等",
             ha="center", fontsize=9,
             bbox=dict(boxstyle="round", fc="white", ec="gray"))

    fig.suptitle("F5  空間運動量はどこにレジストされるか——関係対に分散した登録", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _save(fig, "fig_f5_pair_to_vector_v1")


# ----------------------------------------------------------------------
def fig_f6_architecture() -> None:
    """F6: 二層アーキテクチャ。F=状態のみ、D=体ゲージ保持、G=対読出し。"""
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    ax.axis("off")

    def box(x, y, w, h, text, fc):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04",
                                    fc=fc, ec="black", lw=1.2))
        ax.annotate(text, (x + w / 2, y + h / 2), ha="center", va="center",
                    fontsize=10)

    # 状態層
    box(0.06, 0.62, 0.36, 0.25,
        "状態層（無名・静的な関係配置）\n$\\{z_e\\}$,  $M=N(N-1)/2$", "#dce9f7")
    box(0.58, 0.62, 0.36, 0.25,
        "F 万能相互作用\n（衝突のみ・共変・ゲージを参照しない）\n媒介頂点 O(M) 縮約", "#dce9f7")
    ax.add_patch(FancyArrowPatch((0.42, 0.745), (0.58, 0.745), arrowstyle="<|-|>",
                                 mutation_scale=15, lw=1.5,
                                 transform=ax.transAxes))

    # 読出し層
    box(0.06, 0.16, 0.36, 0.30,
        "D 万能次元読出し\n体ゲージ $\\Phi_A$ を保持・更新\n（絶対中心・主軸枠・接続）", "#e8f5dc")
    box(0.58, 0.16, 0.36, 0.30,
        "G 万能読出し\nゲージ対を通した対読出し\n（$r_A$・対スカラー族・三時計）", "#e8f5dc")

    ax.add_patch(FancyArrowPatch((0.24, 0.62), (0.24, 0.46), arrowstyle="-|>",
                                 mutation_scale=15, lw=1.5,
                                 transform=ax.transAxes))
    ax.annotate("読出しのみ（状態へ書き戻さない）", (0.25, 0.53),
                xycoords="axes fraction", fontsize=9)
    ax.add_patch(FancyArrowPatch((0.76, 0.62), (0.76, 0.46), arrowstyle="-|>",
                                 mutation_scale=15, lw=1.5,
                                 transform=ax.transAxes))
    ax.add_patch(FancyArrowPatch((0.42, 0.31), (0.58, 0.31), arrowstyle="-|>",
                                 mutation_scale=15, lw=1.5,
                                 transform=ax.transAxes))

    ax.annotate("運動学（位置・速度・時間）はすべて読出し層＝ゲージ側に住む。\n"
                "状態側の並進は空虚（CR6/7 で全量凍結）かつ無名性違反（片側並進の実測）——二重に棄却済み。",
                (0.5, 0.045), xycoords="axes fraction", ha="center", fontsize=9,
                bbox=dict(boxstyle="round", fc="#fff6dc", ec="gray"))

    ax.set_title("F6  二層アーキテクチャ——状態は動かさず、ゲージが運動学を担う", fontsize=12)
    _save(fig, "fig_f6_architecture_v1")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    print("基礎論文 新規概念図 F1〜F6 を生成:")
    fig_f1_body_gauge()
    fig_f2_absolute_frame()
    fig_f3_polar_decomposition()
    fig_f4_three_clocks()
    fig_f5_pair_to_vector()
    fig_f6_architecture()
    print("完了")
