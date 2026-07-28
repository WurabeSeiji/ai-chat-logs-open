#!/usr/bin/env python3
"""論文用図版: フェルミオン型（局在移乗）とボゾン型（移乗する局在なし）の衝突比較 v1

データ源（すべて取得済みダンプ、新走行なし）:
  フェルミオン型: 03_oddN_B63_keyR  R=0.6971778791282474（B=フル奇数カーネル、強局在）
  ボゾン型:      03_oddN_B1_keyR   同R（両者とも基底波、局在なし）

図1 fig1_schematic     概念図: 左右から接近→相互作用→反跳。局在が乗り移る/乗り移るものが無い
図2 fig2_waveforms     実データ: 衝突前(v=0)/相互作用中(v≈1/2)/衝突後(v=1) の χ空間密度
図3 fig3_transfer_map  実データ: 局在の往復（チャネル別 χ密度の衝突発展ヒートマップ）

注記: 本モデルの「衝突」は離散交換イベントであり、空間的な接近・反跳は概念図（図1）の
抽象化である。図2・図3は実データの再構成波形。
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
PROD = HERE.parent / "production_dump_v1"

# 検証済みパレット（dataviz 既定、順序固定: slot1=青, slot2=橙）
C_A = "#2a78d6"      # A チャネル
C_B = "#eb6834"      # B チャネル
INK = "#33322e"
MUTED = "#8a887d"
BLUE_RAMP = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#104281"]

plt.rcParams.update({
    "font.family": "Hiragino Sans",
    "text.color": INK, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})

R137_TAG = "R0p6971778791"


def load_case(run: str):
    fn = [p for p in glob.glob(str(PROD / run / "output/harmonic_dump_v1/*.npz")) if R137_TAG in p][0]
    z = np.load(fn)
    return z["coeffs"], z["harmonics"], json.loads(str(z["meta"]))


def density_chi(coeffs_k_ch: np.ndarray, harms: np.ndarray, chi_n: int = 512) -> np.ndarray:
    """chi 周辺密度 ρ(χ) = Σ_η |ψ(χ,η)|² を再構成（合計1に規格化済み状態）。"""
    full = np.zeros((chi_n, coeffs_k_ch.shape[1]), dtype=complex)
    for i, n in enumerate(harms):
        full[int(n) % chi_n, :] = coeffs_k_ch[i, :]
    field = np.fft.ifft(full, axis=0, norm="ortho")
    return np.sum(np.abs(field) ** 2, axis=1)


def transfer_phase(coeffs: np.ndarray) -> np.ndarray:
    a0, b0 = coeffs[0, 0], coeffs[0, 1]
    basis = np.column_stack([(a0 @ a0.conj().T).real.ravel(), (b0 @ b0.conj().T).real.ravel()])
    v = np.zeros(coeffs.shape[0])
    for k in range(coeffs.shape[0]):
        G = (coeffs[k, 0] @ coeffs[k, 0].conj().T).real.ravel()
        sol, *_ = np.linalg.lstsq(basis, G, rcond=None)
        v[k] = sol[1]
    return v


def centered(rho: np.ndarray, shift: int) -> np.ndarray:
    return np.roll(rho, shift)


def main() -> None:
    cf, hf, _ = load_case("03_oddN_B63_keyR")   # フェルミオン型
    cb, hb, _ = load_case("03_oddN_B1_keyR")    # ボゾン型
    chi_n = 512
    chi = np.linspace(-np.pi, np.pi, chi_n, endpoint=False)

    v = transfer_phase(cf)
    k0 = 0
    k_mid = int(np.argmin(np.abs(v[:20] - 0.5)))
    k_full = int(np.argmax(v[:20] > 0.999)) or int(np.argmax(v[:20]))

    # 局在ピークを χ=0 に表示するためのシフト（B63 の B チャネル初期ピーク基準）
    rho_ref = density_chi(cf[0, 1], hf)
    shift = chi_n // 2 - int(np.argmax(rho_ref))

    def panel_data(coeffs, harms, k):
        return (centered(density_chi(coeffs[k, 0], harms), shift),
                centered(density_chi(coeffs[k, 1], harms), shift))

    # ---------- 図2: 実データ波形（2行×3列） ----------
    times = [(k0, "衝突前  (v = 0)"), (k_mid, f"相互作用中  (v ≈ 0.5)"), (k_full, "衝突後  (v ≈ 1)")]
    rows = [("フェルミオン型（B = 奇数63次フルカーネル）: 局在が乗り移る", cf, hf),
            ("ボゾン型（B = 基底波のみ）: 乗り移る局在が無い", cb, hb)]
    ymax = 0.0
    for _, coeffs, harms in rows:
        for k, _t in times:
            for rho in panel_data(coeffs, harms, k):
                ymax = max(ymax, rho.max())
    ymax *= 1.15

    fig, axes = plt.subplots(2, 3, figsize=(11.5, 5.6), sharex=True, sharey=True)
    for r, (row_title, coeffs, harms) in enumerate(rows):
        for c, (k, t_label) in enumerate(times):
            ax = axes[r, c]
            rho_a, rho_b = panel_data(coeffs, harms, k)
            ax.plot(chi, rho_a, color=C_A, lw=2, label="A 波" if (r == 0 and c == 0) else None)
            ax.plot(chi, rho_b, color=C_B, lw=2, label="B 波" if (r == 0 and c == 0) else None)
            ax.set_ylim(0, ymax)
            ax.set_xlim(-np.pi, np.pi)
            if r == 0:
                ax.set_title(t_label, fontsize=11, pad=8)
            if c == 0:
                ax.text(-0.24, 0.5, row_title.split("：")[0].split("（")[0],
                        transform=ax.transAxes, rotation=90, va="center", ha="center", fontsize=11)
            ax.set_xticks([-np.pi, 0, np.pi])
            ax.set_xticklabels(["−π", "0", "π"])
            ax.tick_params(labelsize=9)
        axes[r, 0].set_ylabel("密度 ρ(χ)", fontsize=10)
    # 行の説明・注釈
    axes[0, 2].annotate("局在が A へ移乗", xy=(0.06, 0.80), xycoords="axes fraction",
                        fontsize=10, color=C_A)
    axes[0, 0].annotate("B が局在", xy=(0.60, 0.80), xycoords="axes fraction",
                        fontsize=10, color=C_B)
    axes[1, 1].annotate("両波とも一様なまま（見かけ上、何も起きない）\n※ A と B は同形のため線が完全に重なっている", xy=(0.05, 0.72),
                        xycoords="axes fraction", fontsize=10, color=MUTED)
    axes[2 - 2, 0].set_xlabel("")
    for ax in axes[1]:
        ax.set_xlabel("位相座標 χ", fontsize=10)
    fig.legend(loc="upper right", bbox_to_anchor=(0.995, 1.0), frameon=False, fontsize=10)
    fig.suptitle("同一の相互作用・同一の R での衝突前後比較（実データ, R = R₁₃₇）", fontsize=12, y=1.03)
    fig.text(0.01, -0.04,
             "上段: B=奇数63次フルカーネル（強局在）。衝突後、局在構造は完全に A 側へ移乗する（転送位相 v=1）。\n"
             "下段: B=基底波（局在なし）。同じ相互作用でも、乗り移る構造が無いため波形は変化しない。\n"
             "注: 奇数倍音カーネルは半周期反対称 f(χ+π)=−f(χ) を持つため、密度は π 周期となり主峰（χ=0）と対蹠峰（χ=±π）の二峰に見える。",
             fontsize=9, color=MUTED, va="top")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"fig2_waveforms_before_during_after_v1.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ---------- 図3: 局在の往復ヒートマップ（B63、A/B チャネル） ----------
    K_SHOW = 33
    cmap = LinearSegmentedColormap.from_list("blue_seq", BLUE_RAMP)
    win = np.abs(chi) <= 0.9
    maps = []
    for ch in (0, 1):
        M = np.stack([centered(density_chi(cf[k, ch], hf), shift)[win] for k in range(K_SHOW)])
        maps.append(M)
    vmax = max(m.max() for m in maps)
    peak_a = np.array([centered(density_chi(cf[k, 0], hf), shift)[chi_n // 2] for k in range(K_SHOW)])
    peak_b = np.array([centered(density_chi(cf[k, 1], hf), shift)[chi_n // 2] for k in range(K_SHOW)])
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 6.8), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1, 0.9]})
    for i, (ax, M, name) in enumerate(zip(axes[:2], maps, ("A 波", "B 波"))):
        im = ax.imshow(M.T, aspect="auto", origin="lower", cmap=cmap, vmin=0, vmax=vmax,
                       extent=[-0.5, K_SHOW - 0.5, -0.9, 0.9])
        ax.set_ylabel(f"{name}\nχ", fontsize=10, rotation=0, ha="right", va="center", labelpad=18)
        ax.set_yticks([-0.9, 0, 0.9])
        ax.tick_params(labelsize=9)
    ax3 = axes[2]
    ax3.plot(range(K_SHOW), peak_a, color=C_A, lw=2, marker="o", ms=4, label="A 波")
    ax3.plot(range(K_SHOW), peak_b, color=C_B, lw=2, marker="o", ms=4, label="B 波")
    ax3.set_ylabel("ピーク密度\nρ(χ=0)", fontsize=9, rotation=0, ha="right", va="center", labelpad=18)
    ax3.set_xlabel("衝突回数 k", fontsize=10)
    ax3.legend(frameon=False, fontsize=9, loc="upper right", ncol=2)
    ax3.tick_params(labelsize=9)
    cbar = fig.colorbar(im, ax=axes[:2], fraction=0.03, pad=0.02)
    cbar.set_label("密度 ρ(χ)", fontsize=9)
    fig.suptitle("フェルミオン型: 局在構造が A・B 間を周期的に往復する（B63, R = R₁₃₇, 実データ）",
                 fontsize=12)
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"fig3_transfer_map_v1.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ---------- 図1: 概念図（接近 → 相互作用 → 反跳） ----------
    rho_loc = centered(density_chi(cf[0, 1], hf), shift)   # 実データの局在形状を流用
    lobe = np.abs(chi) <= 0.35                              # 中心ローブのみ（対蹠峰・巻き付き端を除く）
    rho_loc = np.interp(np.linspace(-0.35, 0.35, chi_n), chi[lobe], rho_loc[lobe])
    rho_loc = rho_loc / rho_loc.max()
    rho_flat = np.full(chi_n, 0.10)
    x = np.linspace(0, 1, chi_n)

    def draw_wave(ax, x0, width, shape, color, label=None):
        xs = x0 + (x - 0.5) * width
        ax.plot(xs, shape, color=color, lw=2)
        ax.fill_between(xs, 0, shape, color=color, alpha=0.15, linewidth=0)
        if label:
            ax.text(x0, max(shape) + 0.12, label, ha="center", fontsize=10, color=color)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 5.0))
    scenes = [
        # (row, col, 左波形, 右波形, 左ラベル, 右ラベル, 矢印向き)
        (0, 0, rho_flat, rho_loc, "A（一様波）", "B（局在波束）", "in"),
        (0, 1, rho_loc, rho_flat, "A（局在を獲得）", "B（一様化）", "out"),
        (1, 0, rho_flat, rho_flat * 1.0, "A（一様波）", "B（一様波）", "in"),
        (1, 1, rho_flat, rho_flat * 1.0, "A（変化なし）", "B（変化なし）", "out"),
    ]
    for r, c, left, right, llab, rlab, mode in scenes:
        ax = axes[r, c]
        draw_wave(ax, 0.25, 0.42, left, C_A, llab)
        draw_wave(ax, 0.75, 0.42, right, C_B, rlab)
        y_ar = 0.62
        if mode == "in":
            ax.annotate("", xy=(0.44, y_ar), xytext=(0.30, y_ar),
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
            ax.annotate("", xy=(0.56, y_ar), xytext=(0.70, y_ar),
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
        else:
            ax.annotate("", xy=(0.30, y_ar), xytext=(0.44, y_ar),
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
            ax.annotate("", xy=(0.70, y_ar), xytext=(0.56, y_ar),
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.45)
        ax.axis("off")
    axes[0, 0].set_title("衝突前（接近）", fontsize=11)
    axes[0, 1].set_title("衝突後（反跳）", fontsize=11)
    fig.text(0.02, 0.72, "フェルミオン型\n（局在が乗り移る）", fontsize=11, va="center")
    fig.text(0.02, 0.28, "ボゾン型\n（乗り移る局在なし）", fontsize=11, va="center")
    fig.text(0.02, -0.02,
             "概念図。波束形状は実データ（B63 初期状態）を使用。本モデルの衝突は離散交換イベントであり、"
             "空間的接近・反跳はその抽象化である。",
             fontsize=8.5, color=MUTED)
    fig.tight_layout(rect=[0.12, 0.02, 1, 1])
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"fig1_schematic_v1.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps({"k0": k0, "k_mid": k_mid, "k_full": k_full,
                      "v_mid": float(v[k_mid]), "v_full": float(v[k_full]),
                      "ymax": float(ymax)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
