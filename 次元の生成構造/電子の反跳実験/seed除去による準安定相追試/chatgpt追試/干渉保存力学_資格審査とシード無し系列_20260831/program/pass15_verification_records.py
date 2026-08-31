#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス15（再現性の回復・木原指摘）：チャット内 heredoc / 手計算で行い記録を残していなかった 2 つの検証を
プログラム化し、全数値をファイルに保存する。
(A) W（合成波）と辺読みの照合（harmonic_ladder.md 訂正節の根拠）：mp_N3 / ne_N4 の終窓について、
    候補周波数ごとの辺別 DTFT 振幅と、線ごとの同相度 |Σ_e c_e|/Σ_e|c_e| を全列挙。
(B) 側帯理論の照合（wavelength_table.md 追記節の根拠）：λ_sat/λ_star の理論
    (N−2)/((N−2)−√((N−2)²−4)) と実測（wavelength_table.csv の波長グループ）、
    および mp_N5 の速い枝 (Ω+√(Ω²−P̃²))a²Δ と実測線（harmonic_ladder.csv）の対照。"""
import os, csv, math, json
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
LW = 16384

md = ["# 検証記録（パス15）：W×辺読みの照合と側帯理論の照合（全数値）", ""]

# ---- (A) 照合 ----
CAND = {"mp_N3": [0.0033781, 0.007661, 0.019363, 0.011283],
        "ne_N4": [-0.0208989, 0.015382, 0.010566]}
rowsA = []
for tag, oms in CAND.items():
    Z = np.load(os.path.join(DATA, tag, "states_treatment.npz"))["Z"][-LW:]
    M = Z.shape[1]; t = np.arange(LW); h = np.hanning(LW)
    md.append(f"## (A) {tag}：辺別 DTFT 振幅（規格化=その辺の最大値）と線の同相度")
    md.append("")
    md.append("| 辺 | " + " | ".join(f"ω={om:+.6f}" for om in oms) + " |")
    md.append("|---|" + "---|"*len(oms))
    C = np.zeros((M, len(oms)), complex)
    for e in range(M):
        for a, om in enumerate(oms):
            C[e, a] = np.sum(Z[:, e]*h*np.exp(-1j*om*t))
        mx = max(abs(C[e, a]) for a in range(len(oms)))
        md.append(f"| {e} | " + " | ".join(f"{abs(C[e,a])/mx:.4f}" for a in range(len(oms))) + " |")
        for a, om in enumerate(oms):
            rowsA.append(dict(tag=tag, edge=e, omega=om, dtft_abs=abs(C[e, a])))
    md.append("")
    md.append("| ω | 同相度 \\|Σc\\|/Σ\\|c\\| | \\|Σc\\| |")
    md.append("|---|---|---|")
    for a, om in enumerate(oms):
        coh = abs(C[:, a].sum())/max(np.abs(C[:, a]).sum(), 1e-300)
        md.append(f"| {om:+.6f} | {coh:.6f} | {abs(C[:,a].sum()):.3e} |")
        rowsA.append(dict(tag=tag, edge=-1, omega=om, dtft_abs=abs(C[:, a].sum())))
    md.append("")
with open(os.path.join(ROOT, "results", "crosscheck_W_vs_edges.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["tag", "edge", "omega", "dtft_abs"]); w.writeheader(); w.writerows(rowsA)

# ---- (B) 側帯理論 ----
md.append("## (B) 側帯理論の照合")
md.append("")
md.append("| 対象 | 理論 | 実測 | 相対差 |")
md.append("|---|---|---|---|")
rowsB = []
def add(name, theo, meas):
    dev = abs(meas-theo)/abs(theo)
    md.append(f"| {name} | {theo:.6f} | {meas:.6f} | {dev:.2e} |")
    rowsB.append(dict(item=name, theory=theo, measured=meas, rel_dev=dev))
# λ 比（wavelength_table.csv から波長グループの代表値を再計算）
wt = list(csv.DictReader(open(os.path.join(ROOT, "results", "wavelength_table.csv"))))
def lam_ratio(tag):
    lams = sorted(float(r["lam_norm"]) for r in wt if r["tag"] == tag)
    lo = [x for x in lams if x < 2]; hi = [x for x in lams if x >= 2]
    return (sum(hi)/len(hi))/(sum(lo)/len(lo)) if hi else None
th5 = 3.0/(3.0-math.sqrt(5.0)); th12 = 10.0/(10.0-math.sqrt(96.0))
for tag in ("mp_N5", "ne_N5", "rb_N5"):
    add(f"λ比 {tag}", th5, lam_ratio(tag))
for tag in ("ne_N12", "rb_N12", "mp_N12", "hm_N12"):
    add(f"λ比 {tag}", th12, lam_ratio(tag))
# mp_N5 の速い枝：Ω=(N−2)a²、P̃=2a²、a²Δ=|ν_star|/(N−2)
nu_star = 0.0326647  # wavelength_table の mp_N5 星辺（λ=1 群の |ν|）
a2D = nu_star/3.0
fast_theory = (3.0+math.sqrt(5.0))*a2D
add("mp_N5 速い枝 |ω|（理論 (3+√5)a²Δ）", fast_theory, 0.057009)
slow_theory = (3.0-math.sqrt(5.0))*a2D
add("mp_N5 遅い枝 |ω|（理論 (3−√5)a²Δ）", slow_theory, 0.008313)
md.append("")
md.append(f"理論式：λ_sat/λ_star = (N−2)/((N−2)−√((N−2)²−4))。N=5: {th5:.6f}（= (3/2)φ²）、N=12: {th12:.6f}。")
with open(os.path.join(ROOT, "results", "sideband_theory_check.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["item", "theory", "measured", "rel_dev"]); w.writeheader(); w.writerows(rowsB)
with open(os.path.join(ROOT, "results", "verification_records.md"), "w") as f:
    f.write("\n".join(md)+"\n")
for r in rowsB:
    print(f"{r['item']}: 理論 {r['theory']:.6f} 実測 {r['measured']:.6f} 差 {r['rel_dev']:.2e}")
print("PASS15 OK")
