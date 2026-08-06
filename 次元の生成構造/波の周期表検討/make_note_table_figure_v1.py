#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""note用: 波の周期表の図版化（原始周期表＋標準模型62種割当を一枚に）"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for f in ("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
          "/System/Library/Fonts/Hiragino Sans GB.ttc"):
    if Path(f).exists():
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False

C_E = "#d9f0dd"; C_S = "#fdf3d0"; C_H = "#e8e8ec"; C_HD = "#3a4a6b"; C_HD2 = "#6b7a9b"

fig = plt.figure(figsize=(11, 12.5))
fig.suptitle("波の周期表（第0.4版）", fontsize=20, y=0.995, fontweight="bold")
fig.text(0.5, 0.962, "粒子種 ＝ 状態側の番地 × 観測時計 × 海との関係", ha="center", fontsize=13)

# ---- 上段: 原始周期表 ----
ax1 = fig.add_axes([0.03, 0.635, 0.94, 0.27]); ax1.axis("off")
ax1.set_title("表1  原始周期表（模型固有・実測）——標準模型の名前を使わない主表",
              fontsize=13, loc="left", pad=10)
head1 = ["生巻き類", "読める電荷\n(mod3)", "中性化\n歩数", "時計被覆\n(F帯/B帯)", "海中質量²\n(約数類)", "孤立\n安定性", "海中寿命"]
rows1 = [
 ["m=0（海・基底種）", "0（中性・自由）", "0", "2 / 1", "—", "安定", "基底種は安定（光的）"],
 ["奇数類 {±1,±3,±5,±7}", "±1 または 0", "4（最大）", "2 / 1", "0.787（5桁一致）", "厳密安定", "τ≈10⁴衝突（帯電）"],
 ["{±2, ±6}", "∓1 または 0", "3", "2 / 1", "0.764", "厳密安定", "保持率 0.226"],
 ["{±4}", "±1", "2", "2 / 1", "0.722", "厳密安定", "保持率 0.354（最長）"],
 ["{+8}", "−1", "1（最短）", "—", "未測", "—", "未測"],
]
t1 = ax1.table(cellText=rows1, colLabels=head1, loc="center", cellLoc="center")
t1.auto_set_font_size(False); t1.set_fontsize(10.5); t1.scale(1, 2.1)
for j in range(len(head1)):
    c = t1[0, j]; c.set_facecolor(C_HD); c.set_text_props(color="white", fontsize=10)
for i in range(1, len(rows1) + 1):
    for j in range(len(head1)):
        t1[i, j].set_facecolor("#f4f7fb" if i % 2 else "white")

# ---- 下段: 62種割当 ----
ax2 = fig.add_axes([0.03, 0.10, 0.94, 0.46]); ax2.axis("off")
ax2.set_title("表2  標準模型62種への割当案（仮説）——確度: E=実測錨／S=構造対応／H=仮説",
              fontsize=13, loc="left", pad=10)
head2 = ["粒子", "状態数", "巻き m（電荷 Q=m/3）", "統計（被覆）", "確度", "備考"]
rows2 = [
 ["u, c, t", "18", "+2（+2/3）", "フェルミオン（2）", "E", "閉じ込め＝mod3非可読・色＝隠れ残差", C_E],
 ["d, s, b", "18", "−1（−1/3）", "フェルミオン（2）", "E", "同上", C_E],
 ["e, μ, τ", "6", "−3（−1）", "フェルミオン（2）", "E", "自由・素電荷番地 sin²(23π/124)", C_E],
 ["ν×3", "6", "0（0）", "フェルミオン（2）", "E", "中性フェルミオンの成立を実測", C_E],
 ["光子 γ", "1", "0", "ボゾン（1）", "E", "基底種 ρ=1/1・無質量・安定", C_E],
 ["グルーオン g", "8", "色対（mod3中性）", "ボゾン（1）", "S", "色非一重項→単独非可読", C_S],
 ["W±, Z", "3", "±3, 0", "ボゾン（1）", "S", "t巻き（時計離調）→有質量", C_S],
 ["ヒッグス H", "1", "0", "ボゾン（1）", "S", "海（凝縮体）の集団モード", C_S],
 ["重力子 G", "1", "0", "ボゾン（1）", "H", "ℓ=1量子2個の複合のみ（予言）", C_H],
 ["（世代軸）", "—", "—", "—", "H", "3世代の分化軸は未同定", C_H],
]
t2 = ax2.table(cellText=[r[:6] for r in rows2], colLabels=head2, loc="center", cellLoc="center")
t2.auto_set_font_size(False); t2.set_fontsize(10.5); t2.scale(1, 2.0)
for j in range(len(head2)):
    c = t2[0, j]; c.set_facecolor(C_HD2); c.set_text_props(color="white", fontsize=10)
for i, r in enumerate(rows2, start=1):
    for j in range(len(head2)):
        t2[i, j].set_facecolor(r[6])
widths2 = [0.13, 0.07, 0.20, 0.16, 0.06, 0.38]
for j, w in enumerate(widths2):
    for i in range(len(rows2) + 1):
        t2[i, j].set_width(w)
widths1 = [0.20, 0.13, 0.09, 0.12, 0.15, 0.09, 0.22]
for j, w in enumerate(widths1):
    for i in range(len(rows1) + 1):
        t1[i, j].set_width(w)

fig.text(0.03, 0.035, "合計62種＝クォーク36＋レプトン12＋グルーオン8＋γ＋W2＋Z＋H＋G（反粒子は巻き m→−m）",
         fontsize=11)
fig.text(0.03, 0.015, "電荷=3で読む時計の読み値／統計=時計の二重被覆／質量・寿命=海との関係（詳細は論文 doi:10.5281/zenodo.21822358）",
         fontsize=9.5, color="#555555")
out = Path(__file__).resolve().parent / "fig_p11_periodic_table_visual_v1.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved:", out)
