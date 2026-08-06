#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""「空間軸3軸と固有時間の創生」図7枚

全データは保存済みJSON/npz（対照検定つき再走行系列）から読む。
fig_s1: 実体の2+1構造（特異値の対・接線固有値）
fig_s2: quadrature読出し（静止点・変位の保持）
fig_s3: 単一時計と引き込み
fig_s4: 次元の鋭さ 1/M と少数体領域
fig_s5: 計器が作った偽コム（分解能の教訓）
fig_s6: 質量=非コヒーレンス検証
fig_s7: 閉塞の錐と読出し構造（概念図）
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

plt.rcParams["font.family"] = ["Hiragino Sans", "Arial Unicode MS", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False
HERE = Path(__file__).resolve().parent
J = lambda n: json.loads((HERE / n).read_text())

v1 = J("pre_2plus1_structure_result_v1.json")
v2 = J("pre_chirality_ladder_resolution_result_v2.json")
v3b = J("pre_signed_comb_result_v3b.json")
v3c = J("pre_comb_highres_result_v3c.json")
mc = J("mass_noncoherence_check_result_v1.json")
ser = np.load(HERE / "pre_kinematics_series_v5.npz")

# ============ fig_s1 実体の2+1構造 ============
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
sv = v1["P1"]["sv_rel_top8"]
colors = ["tab:red", "tab:red", "tab:blue", "tab:blue", "tab:green", "tab:green", "tab:gray", "tab:gray"]
ax1.bar(range(1, 9), sv, color=colors)
for i in range(0, 8, 2):
    ax1.annotate("対", (i + 1.5, sv[i] + 0.02), ha="center", fontsize=9)
ax1.set_xlabel("特異値番号"); ax1.set_ylabel("σ/σ₁")
ax1.set_title("(a) 占有構造は厳密に対で並ぶ＝回転平面の階層（N=5）", fontsize=10)
eigs = v1["P2"]["top_eigs"]
th = np.linspace(0, 2 * np.pi, 200)
ax2.plot(np.cos(th), np.sin(th), "k:", lw=0.7)
for e in eigs:
    ms = 12 if e["overlap"] > 0.6 else 6
    c = "tab:red" if e["overlap"] > 0.6 else "tab:gray"
    ax2.plot(e["lam_re"], e["lam_im"], "o", color=c, ms=ms)
ax2.set_xlim(0.9800, 1.0100); ax2.set_ylim(-0.06, 0.06)
ax2.axhline(0, color="k", lw=0.5)
ax2.set_xlabel("Re λ"); ax2.set_ylabel("Im λ")
ax2.set_title("(b) 接線写像固有値：占有を担うのは複素共役対1つ（赤）\n＝回転平面1枚。第3方向は独立固有値を持たない", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "fig_s1_structure_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s1 saved")

# ============ fig_s2 quadrature読出し ============
r1b = ser["base_plane1"]; r1k = ser["kick1_plane1"]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
n_show = len(r1b)
sc = ax1.scatter(r1b.real, r1b.imag, c=np.arange(n_show), cmap="viridis", s=2)
plt.colorbar(sc, ax=ax1, label="サンプル時刻")
ax1.set_xlabel("x = Re(c·e^{−iφ})"); ax1.set_ylabel("y = Im(c·e^{−iφ})")
ax1.set_aspect("equal")
ax1.set_title("(a) quadrature読出し：時計同期成分は一定半径 |r|=0.983±0.005 のまま\n角度のみ緩慢に歳差する（レート0.12%/時計＝振幅不変の準静止）", fontsize=10)
dd = np.abs(r1k) - np.abs(r1b)
tt = np.arange(len(dd))
ax2.plot(tt, dd * 1e3, color="tab:red", lw=1)
ax2.axhline(0, color="k", lw=0.6, ls=":")
ax2.set_xlabel("キック後のサンプル時刻"); ax2.set_ylabel("変位 Δ|r| (×10⁻³)")
ax2.set_title("(b) 静止成分への変位は保持される（8000step・無減衰）", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "fig_s2_quadrature_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s2 saved")

# ============ fig_s3 単一時計と引き込み ============
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9))
for i, (n, d) in enumerate(sorted(v3c.items(), key=lambda kv: int(kv[0]))):
    ratios = [p_["abs_ratio"] for p_ in d["planes"]]
    ax1.plot([int(n)] * len(ratios), ratios, "o", ms=7,
             color=plt.cm.tab10(i), label=f"N={n}")
ax1.axhline(1.0, color="k", ls="--", lw=1)
ax1.set_ylim(0.85, 1.15)
ax1.set_xlabel("N"); ax1.set_ylabel("|ω_平面| / ω_clock")
ax1.set_title("(a) 全占有平面が単一時計で剛体回転\n（高分解能：ω≈π/72、比1±1%）", fontsize=10)
ax1.legend(fontsize=8)
r2b = ser["base_plane2"]; r2k = ser["kick2_plane2"]
om0 = float(ser["om_clock"][0])
ang_b = np.unwrap(np.angle(r2b)) / om0
ang_k = np.unwrap(np.angle(r2k)) / om0
tt = np.arange(len(ang_b))
ax2.plot(tt, ang_b - ang_b[0], color="tab:gray", lw=1.2,
         label="基線の自然揺らぎ（位相スリップを含む非同期）")
ax2.plot(tt, ang_k - ang_k[0], color="tab:red", lw=1.2,
         label="コヒーレントキック後（スリップなし・時計にロック −0.05%）")
ax2.set_xlabel("サンプル時刻"); ax2.set_ylabel("復調位相（時計単位）")
ax2.set_title("(b) 引き込み：摂動は凝縮体の固有時計に同期させられる", fontsize=10)
ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(HERE / "fig_s3_clock_entrainment_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s3 saved")

# ============ fig_s4 次元の鋭さ 1/M ============
sweep = v2["C3"]["sweep"]
fig, ax = plt.subplots(figsize=(7, 4.2))
Ns = [r["N"] for r in sweep]; mis = [r["misalign"] for r in sweep]
Ms = [n * (n - 1) // 2 for n in Ns]
ok = [(n, m_, ms_) for n, m_, ms_ in zip(Ns, Ms, mis) if n >= 5]
bad = [(n, m_, ms_) for n, m_, ms_ in zip(Ns, Ms, mis) if n < 5]
ax.loglog([m_ for _, m_, _ in ok], [ms_ for _, _, ms_ in ok], "o-",
          color="tab:blue", ms=9, label="空間が存在する側（N≥5）")
mm = np.array([8, 40])
ax.loglog(mm, 0.13 / mm, "k--", lw=1, label="~1/M（関係波数の逆数）")
for n, m_, ms_ in bad:
    ax.loglog([m_], [ms_], "x", color="tab:red", ms=14, mew=3,
              label=f"N={n}: 第3次元が結晶化しない（少数体領域）")
for n, m_, ms_ in ok + bad:
    ax.annotate(f"N={n}", (m_, ms_ * 1.3), ha="center", fontsize=9)
ax.set_xlabel("M = N(N−1)/2（関係波の本数）")
ax.set_ylabel("次元の曖昧さ 1 − |n̂·â|")
ax.set_title("次元の鋭さは分解能とともに増す——曖昧さ ~1/M\nN=4 では第3次元が結晶化しない（少数体領域は本稿の対象外）", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(HERE / "fig_s4_sharpness_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s4 saved")

# ============ fig_s5 計器が作った偽コム ============
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.9), sharey=True)
d5 = v3b["5"]
r5 = [p_["abs_ratio"] for p_ in d5["planes"]][:4]
ax1.plot(range(1, 5), r5, "s-", color="tab:red", ms=8)
for i in range(3):
    ax1.annotate(f"Δ={r5[i]-r5[i+1]:.4f}", (i + 1.5, (r5[i] + r5[i+1]) / 2 + 0.01),
                 ha="center", fontsize=9)
ax1.set_ylim(0.75, 1.15)
ax1.set_xlabel("平面番号"); ax1.set_ylabel("|ω| / ω_clock")
ax1.set_title("(a) 低分解能（窓2000step）：0.1%精度の「等間隔コム」\n——実はFFTビン幅そのもの（偽構造）", fontsize=10)
d5h = v3c["5"]
r5h = [p_["abs_ratio"] for p_ in d5h["planes"]]
ax2.plot(range(1, len(r5h) + 1), r5h, "o-", color="tab:blue", ms=8)
ax2.axhline(1.0, color="k", ls="--", lw=1)
ax2.set_xlabel("平面番号")
ax2.set_title("(b) 高分解能（窓10倍＋補間）：等間隔コムは消滅。\n上位平面は単一時計に一致（下位の残差はN=5の少数体性）", fontsize=10)
fig.suptitle("分解能の限界は構造をでっち上げる——本稿の主題が計器に現れた実例（反証と修正の記録）", fontsize=11)
fig.tight_layout()
fig.savefig(HERE / "fig_s5_false_comb_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s5 saved")

# ============ fig_s6 質量=非コヒーレンス検証 ============
fig, ax = plt.subplots(figsize=(7, 4.0))
labels = ["rank-1 合成\n（光的・機械零）", "非コヒーレント合成\n（独立2成分）",
          "親固有モード\n（凝縮前・光的）", "準安定凝縮体\n（質量あり）"]
vals = [max(mc["MC1"]["mass2_over_T2"], 1e-17), mc["MC2"]["mass2_over_T2"],
        mc["MC3"]["parent_lightlike"]["mass2_over_T2"], mc["MC3"]["metastable"]["mass2_over_T2"]]
cols = ["tab:gray", "tab:gray", "tab:blue", "tab:red"]
ax.bar(labels, vals, color=cols, width=0.55)
ax.set_yscale("log"); ax.set_ylim(1e-17, 10)
for i, v_ in enumerate(vals):
    ax.annotate(f"{v_:.1e}", (i, v_), ha="center", va="bottom", fontsize=9)
ax.set_ylabel("質量²/T² = detΓ/T²（非コヒーレンス）")
ax.set_title("質量読出し＝Gram行列の非コヒーレンス（本稿導出・機械検証）\n"
             "detΓ = T²−X²−Y²−Z² ≥ 0（光錐束縛は自動）", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "fig_s6_mass_check_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s6 saved")

# ============ fig_s7 閉塞の錐と読出し構造（概念図） ============
fig, ax = plt.subplots(figsize=(11, 5.2))
ax.set_xlim(0, 22); ax.set_ylim(0, 11); ax.axis("off")

def box(x, y, w, h, text, fc="#f0f0f0", ec="#555555", fs=10, weight="normal"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15", fc=fc, ec=ec, lw=1.4))
    ax.annotate(text, (x + w / 2, y + h / 2), ha="center", va="center", fontsize=fs, weight=weight)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=15,
                                  lw=1.5, color="#333333"))

ax.annotate("$x^2+y^2+z^2 \\; = \\; t^2+R^2+Q^2 \\; = \\; R'^2$", (11, 10.3),
            ha="center", fontsize=15)
box(1.2, 7.2, 5.6, 1.6, "凝縮体の実体\n回転平面1枚（複素固有対）＋回転軸1本", fc="#fff3c8", fs=9.5)
box(1.2, 4.9, 5.6, 1.4, "固有時計 φ(t)\n単一・引き込みで集団強制", fc="#fff3c8", fs=9.5)
box(1.2, 2.4, 5.6, 1.6, "波形の構造\n非コヒーレンス／巻き数・透過率", fc="#fff3c8", fs=9.5)

box(9.6, 8.0, 4.6, 1.4, "x, y ＝ quadrature\n（時計との内積の直交位相）", fc="#e2f0e2", fs=9)
box(9.6, 6.3, 4.6, 1.2, "z ＝ 回転軸方向の射影", fc="#e2f0e2", fs=9)
box(9.6, 4.6, 4.6, 1.2, "τ ＝ R′時計の刻み（可読）", fc="#e2f0e2", fs=9)
box(9.6, 2.9, 4.6, 1.2, "m ＝ 非コヒーレンス（可読）", fc="#e2f0e2", fs=9)
box(9.6, 1.2, 4.6, 1.2, "q ＝ 巻き数・透過率（可読）", fc="#e2f0e2", fs=9)
arrow(6.8, 8.0, 9.6, 8.6)
arrow(6.8, 7.4, 9.6, 6.9)
arrow(6.8, 5.6, 9.6, 5.2)
arrow(6.8, 3.2, 9.6, 3.5)
arrow(6.8, 2.8, 9.6, 1.8)

box(16.4, 4.4, 5.0, 2.6, "座標時間 t\n唯一の非可読\n残差 t²=R′²−m²−q²\nまたは大系の広域選択", fc="#f7dddd", ec="#a04040", fs=9.5, weight="bold")
ax.annotate("観測の全内容 ＝ 空間3軸 ＋ τ ＋ m ＋ q（＝ちょうど3+1＋物質属性）",
            (11, 0.35), ha="center", fontsize=11)
ax.annotate("固有時間（物理）と座標時間（規約）の区別が\n読出し構造から導出される",
            (18.9, 3.4), ha="center", fontsize=9, color="#a04040")
ax.set_title("閉塞の錐の読出し構造：何が読めて、何が読めないか", fontsize=13)
fig.tight_layout()
fig.savefig(HERE / "fig_s7_cone_readout_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s7 saved")
print("ALL FIGURES DONE")
