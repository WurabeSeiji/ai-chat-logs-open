#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v2周期表の図版化（日英）: 原始周期表＋62種割当＋新・階層構造（バリオン階層・質量・寿命）"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

JA_FONT = None
for f in ("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
          "/System/Library/Fonts/Hiragino Sans GB.ttc"):
    if Path(f).exists():
        font_manager.fontManager.addfont(f)
        JA_FONT = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False

C_E = "#d9f0dd"; C_S = "#fdf3d0"; C_H = "#e8e8ec"
C_HD = "#3a4a6b"; C_HD2 = "#6b7a9b"; C_HD3 = "#7a5b3a"; C_BOX = "#eef3fa"

TXT = {
"ja": dict(
 title="波の周期表 v2",
 subtitle="粒子種 ＝ 状態側の番地 × 観測時計 × 海との関係　　（v2: 質量・寿命・階層を統合）",
 formula="一行定式（v2の幹）:  r = sin²θ = P奇 / (P奇 + P偶)　——　衝突率 ＝ 奇数倍音（フェルミオン帯）パワーの分率",
 t1="表1  原始周期表（模型固有・実測）——標準模型の名前を使わない主表",
 h1=["生巻き類", "読める電荷\n(mod3)", "中性化\n歩数", "時計被覆\n(F帯/B帯)", "海中質量²\n(約数類)", "孤立\n安定性", "海中寿命"],
 r1=[
  ["m=0（海・基底種）", "0（中性・自由）", "0", "2 / 1", "—", "安定", "基底種は安定（光的）"],
  ["奇数類 {±1,±3,±5,±7}", "±1 または 0", "4（最大）", "2 / 1", "0.787（5桁一致）", "厳密安定", "τ≈10⁴衝突（帯電）"],
  ["{±2, ±6}", "∓1 または 0", "3", "2 / 1", "0.764", "厳密安定", "保持率 0.226"],
  ["{±4}", "±1", "2", "2 / 1", "0.722", "厳密安定", "保持率 0.354（最長）"],
  ["{+8}", "−1", "1（最短）", "—", "未測", "—", "未測"]],
 t2="表2  標準模型62種への割当案（仮説）——確度: E=実測錨／S=構造対応／H=仮説",
 h2=["粒子", "状態数", "巻き m（電荷 Q=m/3）", "統計（被覆）", "確度", "備考"],
 r2=[
  ["u, c, t", "18", "+2（+2/3）", "フェルミオン（2）", "E", "閉じ込め＝mod3非可読・色＝隠れ残差", C_E],
  ["d, s, b", "18", "−1（−1/3）", "フェルミオン（2）", "E", "同上", C_E],
  ["e, μ, τ", "6", "−3（−1）", "フェルミオン（2）", "E", "自由・素電荷番地 sin²(23π/124)", C_E],
  ["ν×3", "6", "0（0）", "フェルミオン（2）", "E", "中性フェルミオンの成立を実測", C_E],
  ["光子 γ", "1", "0", "ボゾン（1）", "E", "基底種 ρ=1/1・無質量・安定＝真空固定点", C_E],
  ["グルーオン g", "8", "色対（mod3中性）", "ボゾン（1）", "S", "色非一重項→単独非可読", C_S],
  ["W±, Z", "3", "±3, 0", "ボゾン（1）", "S", "t巻き（時計離調）→有質量", C_S],
  ["ヒッグス H", "1", "0", "ボゾン（1）", "S", "海（凝縮体）の集団モード", C_S],
  ["重力子 G", "1", "0", "ボゾン（1）", "H", "ℓ=1量子2個の複合のみ（予言）", C_H],
  ["（世代軸）", "—", "—", "—", "H", "v2新候補: 世代＝位相ロック水準（基底=完全/第2=部分/第3=辺縁）", C_H]],
 law="質量・寿命の法則（v2・プロトタイプ実測）:  質量 = 〈ω〉（台上時計レートの平均）／ 寿命 = 1/σ_ω（線幅の逆数）／ τ·σ_ω ≈ 一定（CV30%）  ⇒  重いほど短寿命（世代現象論が自動で従う）",
 t3="表3  存在の階層構造（v2新設）——バリオン階層と質量・寿命の座席表",
 h3=["階層", "構成原理", "電荷", "統計", "質量", "寿命・崩壊", "確度"],
 r3=[
  ["海（真空）", "純ボゾン帯（P奇=0）", "0", "—", "時を刻まない\n（r=0の厳密固定点）", "∞", "E", C_E],
  ["素種（62種）", "純パリティ＝対蹠二点構造\n（一点局在は不可能）", "Q=m/3", "確定\n（被覆1/2）", "〈ω〉", "1/σ_ω", "S", C_S],
  ["ハドロン／バリオン", "Σm≡0 (mod3) のωロック類\np=uud:+3→+1／n=udd:0→0", "整数\n（自動）", "複合", "類の〈ω〉", "分裂＝類の分岐\nt_split=π/Δω", "S", C_S],
  ["古典物体・凝縮体", "偶奇等量混合（全倍音・一点局在）", "整数", "なし\n（混合）", "Σ〈ω〉\n（重力源）", "巨視的", "S", C_S]],
 foot1="合計62種＝クォーク36＋レプトン12＋グルーオン8＋γ＋W2＋Z＋H＋G（反粒子は巻き m→−m）",
 foot2="電荷=3で読む時計の読み値／統計=時計の二重被覆（空間では対蹠二点構造）／質量・寿命=時計場ω(x)の平均と幅／分裂=ωロック類の分岐（πの自然基準）",
 foot3="波の周期表 v2（doi:10.5281/zenodo.21822358）　確度: E=実測錨／S=構造対応／H=仮説",
 out="fig_v2_periodic_table_ja_v1.png"),
"en": dict(
 title="The Periodic Table of Waves, v2",
 subtitle="particle species = state-side address × observation clock × relation to the sea    (v2: mass, lifetime, hierarchy unified)",
 formula="The one-line law (the trunk of v2):  r = sin²θ = P_odd / (P_odd + P_even)  —  collision rate = fraction of odd-harmonic (fermionic-band) power",
 t1="Table 1  Native periodic table (model-intrinsic, measured) — the primary table without Standard-Model names",
 h1=["Raw winding class", "Readable charge\n(mod 3)", "Steps to\nneutralize", "Clock covering\n(F/B band)", "In-sea mass²\n(divisor class)", "Isolated\nstability", "In-sea lifetime"],
 r1=[
  ["m=0 (sea, ground species)", "0 (neutral, free)", "0", "2 / 1", "—", "stable", "ground species stable (lightlike)"],
  ["odd class {±1,±3,±5,±7}", "±1 or 0", "4 (max)", "2 / 1", "0.787 (5-digit match)", "exactly stable", "τ≈10⁴ collisions (charged)"],
  ["{±2, ±6}", "∓1 or 0", "3", "2 / 1", "0.764", "exactly stable", "retention 0.226"],
  ["{±4}", "±1", "2", "2 / 1", "0.722", "exactly stable", "retention 0.354 (longest)"],
  ["{+8}", "−1", "1 (shortest)", "—", "unmeasured", "—", "unmeasured"]],
 t2="Table 2  Assignment of the 62 Standard-Model species (hypothesis) — confidence: E=measured anchor / S=structural / H=hypothesis",
 h2=["Particle", "States", "Winding m (charge Q=m/3)", "Statistics (cover)", "Conf.", "Notes"],
 r2=[
  ["u, c, t", "18", "+2 (+2/3)", "fermion (2)", "E", "confinement = mod-3 unreadability; color = hidden residue", C_E],
  ["d, s, b", "18", "−1 (−1/3)", "fermion (2)", "E", "same as above", C_E],
  ["e, μ, τ", "6", "−3 (−1)", "fermion (2)", "E", "free; elementary-charge address sin²(23π/124)", C_E],
  ["ν×3", "6", "0 (0)", "fermion (2)", "E", "neutral fermion established by measurement", C_E],
  ["photon γ", "1", "0", "boson (1)", "E", "ground species ρ=1/1, massless, stable = vacuum fixed point", C_E],
  ["gluons g", "8", "color pair (mod-3 neutral)", "boson (1)", "S", "non-singlet color → unreadable alone", C_S],
  ["W±, Z", "3", "±3, 0", "boson (1)", "S", "t-winding (clock detuning) → massive", C_S],
  ["Higgs H", "1", "0", "boson (1)", "S", "collective mode of the sea (condensate)", C_S],
  ["graviton G", "1", "0", "boson (1)", "H", "only as a composite of two ℓ=1 quanta (prediction)", C_H],
  ["(generation axis)", "—", "—", "—", "H", "v2 candidate: generations = phase-lock levels (full/partial/marginal)", C_H]],
 law="Mass-lifetime law (v2, prototype):  mass = ⟨ω⟩ (mean clock rate) / lifetime = 1/σ_ω (inverse linewidth) / τ·σ_ω ≈ const (CV 30%)  ⇒  heavier ⇒ shorter-lived (generation phenomenology automatic)",
 t3="Table 3  The hierarchy of existence (new in v2) — seating chart of the baryon hierarchy, mass, and lifetime",
 h3=["Level", "Organizing principle", "Charge", "Statistics", "Mass", "Lifetime / decay", "Conf."],
 r3=[
  ["sea (vacuum)", "pure bosonic band (P_odd = 0)", "0", "—", "does not tick\n(exact fixed point, r=0)", "∞", "E", C_E],
  ["elementary (62)", "pure parity = antipodal two-point\nstructure (no point localization)", "Q=m/3", "definite\n(cover 1/2)", "⟨ω⟩", "1/σ_ω", "S", C_S],
  ["hadrons / baryons", "ω-locked class with Σm≡0 (mod 3)\np=uud:+3→+1 / n=udd:0→0", "integer\n(automatic)", "composite", "⟨ω⟩ of the class", "splitting = class bifurcation\nt_split = π/Δω", "S", C_S],
  ["classical bodies", "equal parity mixture\n(all harmonics, point-localized)", "integer", "none\n(mixed)", "Σ⟨ω⟩\n(gravity source)", "macroscopic", "S", C_S]],
 foot1="62 in total = 36 quarks + 12 leptons + 8 gluons + γ + 2 W + Z + H + G (antiparticles: winding m→−m)",
 foot2="charge = reading of the divide-by-3 clock / statistics = double cover of the clock (spatially: antipodal two-point structure) / mass & lifetime = mean and width of the clock field ω(x) / splitting = bifurcation of ω-locked classes (natural π criterion)",
 foot3="The Periodic Table of Waves v2 (doi:10.5281/zenodo.21822358)   confidence: E=measured anchor / S=structural / H=hypothesis",
 out="fig_v2_periodic_table_en_v1.png"),
}

def make(lang):
    d = TXT[lang]
    plt.rcParams["font.family"] = JA_FONT if lang == "ja" else "DejaVu Sans"
    fig = plt.figure(figsize=(12.5, 14.2))
    fig.suptitle(d["title"], fontsize=19, y=0.998, fontweight="bold")
    fig.text(0.5, 0.972, d["subtitle"], ha="center", fontsize=11.5)
    fig.text(0.5, 0.948, d["formula"], ha="center", fontsize=11,
             bbox=dict(boxstyle="round,pad=0.5", fc=C_BOX, ec="#8899bb"))

    # 表1
    ax1 = fig.add_axes([0.03, 0.775, 0.94, 0.145]); ax1.axis("off")
    ax1.set_title(d["t1"], fontsize=12.5, loc="left", pad=8)
    t1 = ax1.table(cellText=d["r1"], colLabels=d["h1"], loc="center", cellLoc="center")
    t1.auto_set_font_size(False); t1.set_fontsize(9.8); t1.scale(1, 1.75)
    for j in range(len(d["h1"])):
        c = t1[0, j]; c.set_facecolor(C_HD); c.set_text_props(color="white", fontsize=9.3)
    for i in range(1, len(d["r1"]) + 1):
        for j in range(len(d["h1"])):
            t1[i, j].set_facecolor("#f4f7fb" if i % 2 else "white")
    for j, w in enumerate([0.20, 0.13, 0.09, 0.12, 0.15, 0.09, 0.22]):
        for i in range(len(d["r1"]) + 1):
            t1[i, j].set_width(w)

    # 表2
    ax2 = fig.add_axes([0.03, 0.475, 0.94, 0.275]); ax2.axis("off")
    ax2.set_title(d["t2"], fontsize=12.5, loc="left", pad=8)
    t2 = ax2.table(cellText=[r[:6] for r in d["r2"]], colLabels=d["h2"], loc="center", cellLoc="center")
    t2.auto_set_font_size(False); t2.set_fontsize(9.8); t2.scale(1, 1.9)
    for j in range(len(d["h2"])):
        c = t2[0, j]; c.set_facecolor(C_HD2); c.set_text_props(color="white", fontsize=9.3)
    for i, r in enumerate(d["r2"], start=1):
        for j in range(len(d["h2"])):
            t2[i, j].set_facecolor(r[6])
    for j, w in enumerate([0.12, 0.06, 0.19, 0.14, 0.05, 0.44]):
        for i in range(len(d["r2"]) + 1):
            t2[i, j].set_width(w)

    # 法則ボックス
    fig.text(0.5, 0.452, d["law"], ha="center", fontsize=10,
             bbox=dict(boxstyle="round,pad=0.5", fc="#fdf3d0", ec="#bb9944"), wrap=True)

    # 表3
    ax3 = fig.add_axes([0.03, 0.23, 0.94, 0.175]); ax3.axis("off")
    ax3.set_title(d["t3"], fontsize=12.5, loc="left", pad=8)
    t3 = ax3.table(cellText=[r[:7] for r in d["r3"]], colLabels=d["h3"], loc="center", cellLoc="center")
    t3.auto_set_font_size(False); t3.set_fontsize(9.8); t3.scale(1, 2.6)
    for j in range(len(d["h3"])):
        c = t3[0, j]; c.set_facecolor(C_HD3); c.set_text_props(color="white", fontsize=9.3)
    for i, r in enumerate(d["r3"], start=1):
        for j in range(len(d["h3"])):
            t3[i, j].set_facecolor(r[7])
    for j, w in enumerate([0.14, 0.30, 0.08, 0.09, 0.14, 0.17, 0.05]):
        for i in range(len(d["r3"]) + 1):
            t3[i, j].set_width(w)

    fig.text(0.03, 0.19, d["foot1"], fontsize=10.5)
    fig.text(0.03, 0.17, d["foot2"], fontsize=9, color="#555555")
    fig.text(0.03, 0.152, d["foot3"], fontsize=9, color="#555555")
    out = Path(__file__).resolve().parent / d["out"]
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out)

for lang in ("ja", "en"):
    make(lang)
