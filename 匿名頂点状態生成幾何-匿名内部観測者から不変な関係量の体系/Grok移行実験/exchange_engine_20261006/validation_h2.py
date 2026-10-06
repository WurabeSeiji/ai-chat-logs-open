#!/usr/bin/env python3
"""H₂ 部品の検証：文献値との対照。出力 validation_h2.md
  1. BO フィット vs Pachucki 2010 の表（82 点）
  2. 準位 vs H2SPECTRE 7.4 の E2(FULL)（非相対論・非断熱まで、93 準位）と Roueff+2019（CDS、全補正、302 準位）
  3. 遷移 vs Roueff+2019 の σ・A_qua・A_ma（4712 行）
  4. 格子の収束、オルト／パラの分離
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import system_h2 as S  # noqa: E402

CM = S.HARTREE_CM
L = ["# H₂ 部品の検証：文献値との対照", "",
     "準位：BO（Pachucki 2010）＋断熱補正（Pachucki–Komasa 2014）、μ_n = m_p/2、sinc-DVR dr = 0.025 bohr、R ≤ 50 bohr。",
     "率：Θ(R) = Q_WSD(R)/2（WSD 1998）、g(R)（PK 2011）。比較対象の Roueff+2019 は V = BO＋断熱（同じ）、Θ・g は Pachucki–Komasa の値、σ は PK 2018 の全補正（非断熱・相対論・QED）。", ""]

# 1. BO フィット
R, E = np.loadtxt(HERE / "h2_data/h2_V_BO_pachucki2010.dat", unpack=True)
d = S.v_el(R) - 1.0 - E
L += ["## 1. BO ポテンシャルのフィット vs Pachucki 2010 表 III（82 点）", "",
      "最大 |差| = %.2e hartree（R = %.2f）、rms = %.2e hartree。表の精度は ~1e-15、フィットの精度が ~1e-12。" % (np.max(np.abs(d)), R[np.argmax(np.abs(d))], np.sqrt(np.mean(d * d))),
      "断熱補正 E_a(1.4) = %.4f cm⁻¹（2H 原点）、E_a(R→∞) → 0（R = 30 で %.1e cm⁻¹）。" % (S.e_ad(1.4) * CM, S.e_ad(30.0) * CM), ""]

# 準位・遷移の計算
states, channels, init, aux = S.build()
levels = aux["levels"]; idx = aux["idx"]
Emine = {s["label"]: s["E_eV"] / S.HARTREE_EV * CM for s in states}       # cm⁻¹、2H 原点（負）

# 2a. H2SPECTRE E2(FULL)
rows = np.loadtxt(HERE / "h2_data/validation/h2spectre74_H2_E2full.dat")
L += ["## 2. 準位", "", "### 2a. vs H2SPECTRE 7.4 data/H2.dat の E2(FULL)（非相対論、非断熱込み）", "",
      "H2.dat は全エネルギー（hartree）。解離エネルギー D = −E − m_p/(m_p+1)（有限質量の 2H）。本計算との差 = 非断熱補正（本計算に無い）。", "",
      "| v | 準位数 | 本計算 − E2(FULL) [cm⁻¹] 最小 | 最大 | (v,0) の値 |", "|---|---|---|---|---|"]
byv = {}
for v, J, Eh in rows:
    key = (int(v), int(J))
    if key not in Emine:
        continue
    Dfull = -(Eh + S.MPROT / (S.MPROT + 1)) * CM
    byv.setdefault(int(v), []).append((key, Emine[key] + Dfull))      # Emine = −D_mine
for v in sorted(byv):
    diffs = [x[1] for x in byv[v]]
    j0 = [x[1] for x in byv[v] if x[0][1] == 0]
    L.append("| %d | %d | %+.4f | %+.4f | %s |" % (v, len(diffs), min(diffs), max(diffs), ("%+.4f" % j0[0]) if j0 else ""))
D0_mine = -Emine[(0, 0)]
out24 = [x[1] for x in byv.get(2, []) if x[0] == (2, 24)]
L += ["", "D₀：本計算 %.4f、E2(FULL) 36118.7977、Roueff+2019（全補正）36118.0695 cm⁻¹。(0,1)−(0,0)：本計算 %.4f、H2SPECTRE 全補正 118.4868 cm⁻¹。" % (D0_mine, Emine[(0, 1)] - Emine[(0, 0)]),
      "v = 2 の最大 %+.1f cm⁻¹ は (2,24) 一点だけ（H2.dat の −1.02598733666042 hartree は隣の J から外れており、配布データ側の誤記と見る。2b の CDS との比較では (2,24) も +1.3〜+3.4 の範囲内）。" % (out24[0] if out24 else float("nan")), ""]

# 2b. CDS
cds = []
for line in (HERE / "h2_data/validation/roueff2019_table2.dat").read_text().splitlines():
    f = line.split()
    if len(f) < 16:
        continue
    cds.append(dict(vu=int(f[0]), Ju=int(f[1]), vl=int(f[2]), Jl=int(f[3]), sigma=float(f[4]), Aq=float(f[8]), Am=float(f[9]), Eup=float(f[12])))
Ecds = {}
for r in cds:
    Ecds[(r["vu"], r["Ju"])] = r["Eup"]
    Ecds.setdefault((r["vl"], r["Jl"]), r["Eup"] - r["sigma"])
missing = sorted(set(Ecds) - set(Emine)); extra = sorted(set(Emine) - set(Ecds))
L += ["### 2b. vs Roueff+2019（CDS J/A+A/630/A58、準位エネルギーは PK 2018 の全補正）", "",
      "CDS の準位 %d、本計算の束縛準位 %d。CDS にあって本計算に無い：%s。本計算にあって CDS に無い：%s。" % (len(Ecds), len(Emine), missing or "なし", extra or "なし"), "",
      "| v | 準位数 | 本計算 − CDS [cm⁻¹] 最小 | 最大 | (v,0) |", "|---|---|---|---|---|"]
byv = {}
for key, e in Ecds.items():
    if key in Emine:
        byv.setdefault(key[0], []).append((key, Emine[key] - e))
for v in sorted(byv):
    diffs = [x[1] for x in byv[v]]
    j0 = [x[1] for x in byv[v] if x[0][1] == 0]
    L.append("| %d | %d | %+.4f | %+.4f | %s |" % (v, len(diffs), min(diffs), max(diffs), ("%+.4f" % j0[0]) if j0 else ""))
L.append("")

# 3. 遷移
chan = {}
for a, b, kind, k, A in channels:
    chan[(states[a]["label"], states[b]["label"], kind)] = A
rel_s, rel_q, rel_q5, rel_m = [], [], [], []
strong_q, strong_m = [], []
worst_s, worst_q = [], []
n_missing = 0
for r in cds:
    u, l = (r["vu"], r["Ju"]), (r["vl"], r["Jl"])
    if u not in Emine or l not in Emine:
        n_missing += 1
        continue
    s_mine = Emine[u] - Emine[l]
    rel_s.append(s_mine / r["sigma"] - 1.0)
    worst_s.append((abs(rel_s[-1]), u, l, r["sigma"], s_mine))
    if r["Aq"] > 0:
        Aq = chan.get((u, l, "e2"), 0.0)
        rel_q.append(Aq / r["Aq"] - 1.0)
        rel_q5.append(Aq * (r["sigma"] / s_mine) ** 5 / r["Aq"] - 1.0)
        if r["Aq"] > 1e-9:
            strong_q.append(rel_q5[-1])
            worst_q.append((abs(rel_q5[-1]), u, l, r["Aq"], Aq * (r["sigma"] / s_mine) ** 5))
    if r["Am"] > 0:
        Am = chan.get((u, l, "m1"), 0.0)
        rel_m.append(Am * (r["sigma"] / s_mine) ** 3 / r["Am"] - 1.0)
        if r["Am"] > 1e-11:
            strong_m.append(rel_m[-1])
rel_s, rel_q, rel_q5, rel_m = map(np.array, (rel_s, rel_q, rel_q5, rel_m))
worst_s.sort(reverse=True); worst_q.sort(reverse=True)
def stat(x):
    x = np.abs(np.asarray(x))
    return "中央値 %.1e、最大 %.1e（%d 本）" % (np.median(x), np.max(x), len(x))
L += ["## 3. 遷移 vs Roueff+2019（CDS 4712 行、本計算に両準位がある %d 行）" % (len(cds) - n_missing), "",
      "- 波数 σ の相対差 |σ_mine/σ_CDS − 1|：%s。差の中身は非断熱・相対論・QED（本計算に無い）。" % stat(rel_s),
      "- A_qua の相対差（そのまま）：%s。" % stat(rel_q),
      "- A_qua の相対差（σ⁵ を CDS の σ に揃えて、Θ の行列要素だけの差）：%s。A_qua > 1e-9 s⁻¹ の強い線に限ると %s。" % (stat(rel_q5), stat(strong_q)),
      "- A_ma の相対差（σ³ を揃えて、g の行列要素だけの差）：%s。A_ma > 1e-11 s⁻¹ では %s。" % (stat(rel_m), stat(strong_m)), "",
      "CDS の A は有効 4 桁。", "",
      "σ の相対差が最大の 5 本（近接準位間の小さな σ に、準位の数 cm⁻¹ の差が効く）：", ""]
for x, u, l, sc, sm in worst_s[:5]:
    L.append("- (%d,%d)→(%d,%d)：σ_CDS %.4f、σ_mine %.4f cm⁻¹、相対差 %.1e" % (u[0], u[1], l[0], l[1], sc, sm, x))
L += ["", "強い線（A_qua > 1e-9）で Θ 行列要素の差が最大の 5 本：", ""]
for x, u, l, ac, am in worst_q[:5]:
    L.append("- (%d,%d)→(%d,%d)：A_CDS %.3e、A_mine（σ 揃え）%.3e s⁻¹、相対差 %.1e" % (u[0], u[1], l[0], l[1], ac, am, x))
L += ["", "### 3a. 代表線", "", "| 線 | σ_CDS [cm⁻¹] | σ_mine | A_qua CDS [s⁻¹] | A_E2 mine | 比 | A_ma CDS | A_M1 mine | 比 |", "|---|---|---|---|---|---|---|---|---|"]
picks = [("0-0 S(0)", (0, 2), (0, 0)), ("0-0 S(1)", (0, 3), (0, 1)), ("1-0 S(1)", (1, 3), (0, 1)), ("1-0 S(0)", (1, 2), (0, 0)),
         ("1-0 Q(1)", (1, 1), (0, 1)), ("1-0 Q(3)", (1, 3), (0, 3)), ("1-0 O(3)", (1, 1), (0, 3)), ("1-0 O(5)", (1, 3), (0, 5)),
         ("2-0 S(1)", (2, 3), (0, 1)), ("2-1 S(1)", (2, 3), (1, 1)), ("3-0 S(3)", (3, 5), (0, 3)), ("1-0 Q(15)", (1, 15), (0, 15))]
cdsmap = {((r["vu"], r["Ju"]), (r["vl"], r["Jl"])): r for r in cds}
for name, u, l in picks:
    r = cdsmap.get((u, l))
    if r is None or u not in Emine or l not in Emine:
        continue
    Aq, Am = chan.get((u, l, "e2"), 0.0), chan.get((u, l, "m1"), 0.0)
    L.append("| %s | %.4f | %.4f | %.3e | %.3e | %.4f | %s | %s | %s |" % (
        name, r["sigma"], Emine[u] - Emine[l], r["Aq"], Aq, Aq / r["Aq"] if r["Aq"] > 0 else float("nan"),
        ("%.3e" % r["Am"]) if r["Am"] > 0 else "—", ("%.3e" % Am) if Am > 0 else "—", ("%.4f" % (Am / r["Am"])) if r["Am"] > 0 else "—"))
L.append("")

# 4. 収束・オルトパラ
E0, _, _ = S.dvr_levels(0); E0b, _, _ = S.dvr_levels(0, dr=0.02, rmax=60.0)
E4, _, _ = S.dvr_levels(4); E4b, _, _ = S.dvr_levels(4, dr=0.02, rmax=60.0)
cross = sum(1 for a, b, kind, k, A in channels if (states[a]["J"] - states[b]["J"]) % 2)
L += ["## 4. 格子の収束とオルト／パラ", "",
      "J = 0：dr 0.025/R 50 → 0.02/60 で準位の最大変化 %.1e cm⁻¹（%d 準位）。J = 4：%.1e cm⁻¹（%d → %d 準位）。" % (
          np.max(np.abs(E0b[:len(E0)] - E0)) * CM, len(E0), np.max(np.abs(E4b[:len(E4)] - E4)) * CM, len(E4), len(E4b)),
      "(14,4)：CDS では束縛（H2SPECTRE 全補正で解離の 0.027 cm⁻¹ 下）。本計算（BO＋断熱、非断熱なし）では E = %s。" % (
          ("%.4f cm⁻¹（束縛）" % Emine[(14, 4)]) if (14, 4) in Emine else "束縛準位として現れない（E > 0）"),
      "J の偶奇を変える辺の本数：%d（E2・M1・GW はすべて偶奇を保つ）。準位 %d、辺 %d（e2 %d、m1 %d、gw %d）。" % (
          cross, len(states), len(channels), sum(c[2] == "e2" for c in channels), sum(c[2] == "m1" for c in channels), sum(c[2] == "gw" for c in channels)), ""]
out = "\n".join(L) + "\n"
(HERE / "validation_h2.md").write_text(out, encoding="utf-8")
print(out)
