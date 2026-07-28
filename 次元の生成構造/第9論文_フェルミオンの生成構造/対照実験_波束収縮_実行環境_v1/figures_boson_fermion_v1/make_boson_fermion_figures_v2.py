#!/usr/bin/env python3
"""論文用図版 v2: フェルミオン型／ボゾン型衝突比較（5段階・最大振幅正規化・同一初期条件）

v1からの変更（発案者指示）:
  - 両行とも同一の初期条件（A=非局在の一様波、B=局在波束）。
    フェルミオン型 = 局在が乗り移る（衝突中は両方が局在を持つ状態を経る）
    ボゾン型      = そのまま透過（局在は B が持ったまま）
  - 各波を自身の最大振幅で正規化（形状比較。非局在波も見えるように）
  - 時間経過を5段階: 前 → 衝突前 → 衝突中 → 衝突後 → 後

データ源（取得済みダンプ、新走行なし）:
  フェルミオン型: 03_oddN_B63_keyR  R = R₁₃₇厳密値（完全移乗、k=0,2,4,6,8 で v: 0→1）
  ボゾン型:      同じ波形対の R = 1.0（交換なし条件、v ≡ 0、移乗ゼロ）

図1 fig1_schematic_v2   概念図5場面: 接近→衝突→分離。透過（局在はBのまま）vs 移乗（Aが局在を獲得）
図2 fig2_waveforms_v2   実データ5段階: 正規化波形 ρ(χ)/ρmax の推移
図3 は v1 のまま（局在往復ヒートマップ）
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROD = HERE.parent / "production_dump_v1"

C_A = "#2a78d6"      # A 波（青）
C_B = "#eb6834"      # B 波（橙）
INK = "#33322e"
MUTED = "#8a887d"

plt.rcParams.update({
    "font.family": "Hiragino Sans",
    "text.color": INK, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",
})

CHI_N = 512
# 転送位相 v が 0 → 1/4 → 1/2 → 3/4 → 1 に最も近づく衝突時点を全時系列から抽出
# （転送は周期≈5.4衝突で振動するため、単調な経過を示すには非等間隔の時点選択が必要）
STEP_KS = [0, 1, 4, 9, 19]
STEP_LABELS = ["前", "衝突前", "衝突中", "衝突後", "後"]


def load_npz(run: str, r_tag: str):
    fn = [p for p in glob.glob(str(PROD / run / "output/harmonic_dump_v1/*.npz"))
          if p.endswith(f"_{r_tag}_v1.npz")][0]
    z = np.load(fn)
    return z["coeffs"], z["harmonics"]


def density_chi(coeffs_k_ch: np.ndarray, harms: np.ndarray) -> np.ndarray:
    full = np.zeros((CHI_N, coeffs_k_ch.shape[1]), dtype=complex)
    for i, n in enumerate(harms):
        full[int(n) % CHI_N, :] = coeffs_k_ch[i, :]
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


def norm_shape(rho: np.ndarray) -> np.ndarray:
    m = rho.max()
    return rho / m if m > 0 else rho


def main() -> None:
    cf, hf = load_npz("03_oddN_B63_keyR", "R0p697177879128")   # フェルミオン型（移乗）
    cb, hb = load_npz("03_oddN_B63_keyR", "R1")                # ボゾン型（透過・移乗なし）
    chi = np.linspace(-np.pi, np.pi, CHI_N, endpoint=False)
    v = transfer_phase(cf)

    rho_ref = density_chi(cf[0, 1], hf)
    shift = CHI_N // 2 - int(np.argmax(rho_ref))

    def shapes(coeffs, harms, k):
        return (norm_shape(np.roll(density_chi(coeffs[k, 0], harms), shift)),
                norm_shape(np.roll(density_chi(coeffs[k, 1], harms), shift)))

    # ---------- 図2 v2: 実データ 2行×5列（正規化波形） ----------
    rows = [("フェルミオン型\n（局在が乗り移る）", cf, hf, True),
            ("ボゾン型\n（そのまま透過）", cb, hb, False)]
    fig, axes = plt.subplots(2, 5, figsize=(14.0, 5.4), sharex=True, sharey=True)
    for r, (row_label, coeffs, harms, is_f) in enumerate(rows):
        for c, (k, step) in enumerate(zip(STEP_KS, STEP_LABELS)):
            ax = axes[r, c]
            sa, sb = shapes(coeffs, harms, k)
            ax.plot(chi, sb, color=C_B, lw=2.4, ls=(0, (4, 2)), label="B 波（破線）" if (r, c) == (0, 0) else None)
            ax.plot(chi, sa, color=C_A, lw=2, label="A 波" if (r, c) == (0, 0) else None)
            ax.set_ylim(0, 1.18)
            ax.set_xlim(-np.pi, np.pi)
            ax.set_xticks([-np.pi, 0, np.pi])
            ax.set_xticklabels(["−π", "0", "π"])
            ax.tick_params(labelsize=9)
            if r == 0:
                ax.set_title(f"{step}\n(k={k},  v={v[k]:.2f})", fontsize=10.5, pad=6)
            else:
                ax.set_title(f"(k={k})", fontsize=9, pad=4, color=MUTED)
            if c == 0:
                ax.text(-0.42, 0.5, row_label, transform=ax.transAxes,
                        rotation=90, va="center", ha="center", fontsize=11)
                ax.set_ylabel("正規化密度 ρ/ρmax", fontsize=9)
        for ax in axes[1]:
            ax.set_xlabel("位相座標 χ", fontsize=9)
    axes[0, 0].annotate("B が局在", xy=(0.56, 0.86), xycoords="axes fraction", fontsize=9.5, color=C_B)
    axes[0, 2].annotate("両方が局在", xy=(0.54, 0.86), xycoords="axes fraction", fontsize=9.5, color=INK)
    axes[0, 4].annotate("A へ移乗完了", xy=(0.52, 0.86), xycoords="axes fraction", fontsize=9.5, color=C_A)
    axes[1, 4].annotate("B が局在のまま\n（何も乗り移らない）", xy=(0.50, 0.78), xycoords="axes fraction",
                        fontsize=9.5, color=C_B)
    fig.legend(loc="upper right", bbox_to_anchor=(0.995, 1.02), frameon=False, fontsize=10)
    fig.suptitle("同一の波形対（A=広がった基底波, B=局在波束）の衝突経過（実データ・各波を最大振幅で正規化）",
                 fontsize=12.5, y=1.05)
    fig.text(0.01, -0.05,
             "上段: 交換条件 R = R₁₃₇（フェルミオン型）。衝突中（v≈0.5）は両波が局在を分有し、後（v≈1）には局在が完全に A へ移乗する。\n"
             "下段: 交換なし条件 R = 1（ボゾン型）。同じ波形対でも局在は B が持ったまま透過し、全段階で波形が変化しない。\n"
             "各パネルは各波を自身の最大値で正規化した形状比較。奇数倍音カーネルは半周期反対称のため密度は π 周期（χ=0 と ±π の二峰）。\n"
             "転送は周期≈5.4衝突で振動するため、v が 0→1/4→1/2→3/4→1 に最も近い衝突時点（k=0,1,4,9,19）を抽出して表示。",
             fontsize=9, color=MUTED, va="top")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"fig2_waveforms_v2.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ---------- 図1 v2: 概念図 2行×5場面 ----------
    lobe = np.abs(chi) <= 0.35
    peak = np.interp(np.linspace(-0.35, 0.35, CHI_N), chi[lobe],
                     np.roll(rho_ref, shift)[lobe])
    peak = peak / peak.max()
    flat = np.full(CHI_N, 0.13)
    half = 0.5 * peak + 0.5 * flat          # 衝突中の「局在を分有」形状
    x = np.linspace(0, 1, CHI_N)

    def draw(ax, x0, width, shape, color, label=None, label_dy=0.10):
        xs = x0 + (x - 0.5) * width
        ax.plot(xs, shape, color=color, lw=2)
        ax.fill_between(xs, 0, shape, color=color, alpha=0.15, linewidth=0)
        if label:
            ax.text(x0, shape.max() + label_dy, label, ha="center", fontsize=9, color=color)

    def arrow(ax, x0, x1, y=1.28):
        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.5))

    # 位置プラン: A は左→右へ、B は右→左へ通過（5場面）
    posA = [0.16, 0.30, 0.50, 0.70, 0.84]
    posB = [0.84, 0.70, 0.50, 0.30, 0.16]
    fig, axes = plt.subplots(2, 5, figsize=(14.0, 4.6))
    for r, is_f in ((0, True), (1, False)):
        for c in range(5):
            ax = axes[r, c]
            xa, xb = posA[c], posB[c]
            if c == 2:  # 衝突中（重なり）
                if is_f:
                    draw(ax, 0.5, 0.40, half, C_A, "A", 0.16)
                    draw(ax, 0.5, 0.40, half * 0.96, C_B, "B", -0.34)
                else:
                    draw(ax, 0.5, 0.40, flat, C_A, "A", 0.13)
                    draw(ax, 0.5, 0.40, peak, C_B, "B", 0.10)
            else:
                a_shape = peak if (is_f and c >= 3) else flat
                b_shape = flat if (is_f and c >= 3) else peak
                draw(ax, xa, 0.34, a_shape, C_A, "A")
                draw(ax, xb, 0.34, b_shape, C_B, "B")
            if c < 2:
                arrow(ax, posA[c] - 0.06, posA[c] + 0.08)
                arrow(ax, posB[c] + 0.06, posB[c] - 0.08)
            elif c > 2:
                arrow(ax, posA[c] - 0.08, posA[c] + 0.06)
                arrow(ax, posB[c] + 0.08, posB[c] - 0.06)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1.55)
            ax.axis("off")
            if r == 0:
                ax.set_title(STEP_LABELS[c], fontsize=11)
    fig.text(0.015, 0.74, "フェルミオン型\n（局在が乗り移る）", fontsize=10.5, va="center")
    fig.text(0.015, 0.30, "ボゾン型\n（そのまま透過）", fontsize=10.5, va="center")
    fig.text(0.015, -0.03,
             "概念図。両行とも初期条件は同一（A=広がった基底波が右向き、B=局在波束が左向き）。上段は衝突中に局在を分有し、"
             "通過後は A が局在を獲得する（移乗）。下段は局在が B に留まったまま透過する。\n"
             "波束形状は実データ（B63初期状態の中心ローブ）。本モデルの衝突は離散交換イベントであり、空間的接近・通過はその抽象化である。",
             fontsize=8.5, color=MUTED, va="top")
    fig.tight_layout(rect=[0.10, 0.02, 1, 0.97])
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"fig1_schematic_v2.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps({"v_at_steps": {int(k): float(v[k]) for k in STEP_KS}}, ensure_ascii=False))


if __name__ == "__main__":
    main()
