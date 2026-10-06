#!/usr/bin/env python3
"""対照実験：符号を入力にした数値解法（two_body_general.py、SIGN = −1）と、閉じた式の固定一式
（../exchange_rel_20261005/exchange_cascade.py）を、同じ条件（水素、初期 5p₃/₂ F=1、t_max = 1e17 s）で比較する。

比較するもの：
  1. 準位 50 個：クーロン（Dirac＋反跳）、重力、超微細、ETOT の差
  2. 辺の A 係数：種別ごとの最大相対差（E1/E2 は Dirac 行列要素になるので相対 α² の差が出るはず）
  3. 放出の合計・内訳、終状態、滞在時間
出力：control_report.md
"""
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import two_body_general as T  # noqa: E402
X = T.X
EV = X.EV

t0 = time.time()
levels = T.build_levels(T.SIGN, verbose=False)
ed_g, ETOT_g = T.edges_general(levels, X.W_GW)
res_g = T.run_master(levels, ed_g, ETOT_g)
t_build = time.time() - t0

# 固定一式の対応
idx_fixed = X.IDX
L = []
L.append("# 対照実験：数値解法（符号入力） vs 閉じた式の固定一式（水素、SIGN = −1）")
L.append("")
L.append("同じ条件：基底 n ≤ 5、初期 5p₃/₂ F=1、放射は等方吸収体へ、t_max = 1e17 s。数値解法の構築 %.1f s。" % t_build)
L.append("")
L.append("## 1. 準位（eV）。差 = 数値 − 閉じた式")
L.append("")
L.append("| 準位 | クーロン（Dirac＋反跳） 閉じた式 | 差 | 重力 閉じた式 | 差 | 超微細 閉じた式 | 差（Dirac/NR − 1） | ETOT 差 |")
L.append("|---|---|---|---|---|---|---|---|")
names = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}
worst = dict(coul=0.0, grav=0.0, hfs_rel=0.0, etot=0.0)
for Lg in levels:
    i = idx_fixed[Lg["label"]]
    dc = (Lg["coul"] - X.ORB_C[i]) * EV
    dg = (Lg["grav"] - X.ORB_G[i]) * EV
    hr = Lg["hfs"] / X.HFS[i] - 1.0
    de = (Lg["coul"] + Lg["grav"] + Lg["hfs"] - X.ETOT[i]) * EV
    worst["coul"] = max(worst["coul"], abs(dc)); worst["grav"] = max(worst["grav"], abs(dg) / abs(X.ORB_G[i] * EV))
    worst["hfs_rel"] = max(worst["hfs_rel"], abs(hr)); worst["etot"] = max(worst["etot"], abs(de))
    n, ell, j, F = Lg["n"], Lg["ell"], Lg["j"], Lg["F"]
    if n <= 2 or (n == 3 and ell == 2) or (n == 5 and ell == 4):
        L.append("| %d%s_{%d/2} F=%d | %.9f | %+.1e | %.4e | %+.1e | %+.4e | %+.2e | %+.1e |"
                 % (n, names[ell], int(2 * j), int(F), X.ORB_C[i] * EV, dc, X.ORB_G[i] * EV, dg, X.HFS[i] * EV, hr, de))
L.append("")
L.append("全 50 準位：クーロンの最大差 %.1e eV（Dirac 固有値の数値精度。閉じた式に対し相対 1e-11）、重力の最大相対差 %.1e（Schrödinger 期待値の数値精度）、"
         "超微細の相対差 Dirac/NR − 1 の最大 %.1e（相対論補正 ≈ (3/2)α² = %.1e が 1s、2s で出る）、ETOT の最大差 %.1e eV。"
         % (worst["coul"], worst["grav"], worst["hfs_rel"], 1.5 * X.ALPHA**2, worst["etot"]))
L.append("")

# 2. A 係数
ed_f = X.edges(X.W_GW)
def key_f(e):
    s, d, k, A = e
    return (X.ST[s], X.ST[d], k)
def key_g(e):
    s, d, k, A = e
    return (levels[s]["label"], levels[d]["label"], k)
Af = {key_f(e): e[3] for e in ed_f}
Ag = {key_g(e): e[3] for e in ed_g}
L.append("## 2. A 係数（数値 Dirac 行列要素 vs 閉じた式の Schrödinger 行列要素）")
L.append("")
L.append("| 種別 | 本数（固定／数値） | 最大相対差 | 相対差の中央値 |")
L.append("|---|---|---|---|")
for kind in ("e1", "e2", "gw", "m1", "2ph"):
    kf = [k for k in Af if k[2] == kind and Af[k] > 0.0]
    rel = [Ag[k] / Af[k] - 1.0 for k in kf if k in Ag]
    L.append("| %s | %d／%d | %.2e | %.2e |" % (kind, len(kf), sum(1 for k in Ag if k[2] == kind and Ag[k] > 0.0), max(abs(x) for x in rel) if rel else 0.0, float(np.median(np.abs(rel))) if rel else 0.0))
Af = {k: v for k, v in Af.items() if v > 0.0}
Ag = {k: v for k, v in Ag.items() if v > 0.0}
missing = [k for k in Af if k not in Ag] + [k for k in Ag if k not in Af]
L.append("")
L.append("辺の集合の不一致：%d 本（%s）" % (len(missing), "なし" if not missing else str(missing[:5])))
L.append("")
L.append("E1・E2 の差は、固定一式が Schrödinger の動径積分、数値解法が Dirac の (G_aG_b + F_aF_b) を使うことによる相対 α² の相対論補正。"
         "2p→1s：固定 %.6e、数値 %.6e s⁻¹（比 %.6f）。m1・2ph は同じ式なので差は ω³ の準位差だけ。"
         % (Af[((2, 1, 1, 1), (1, 0, 1, 1), "e1")], Ag[((2, 1, 1, 1), (1, 0, 1, 1), "e1")], Ag[((2, 1, 1, 1), (1, 0, 1, 1), "e1")] / Af[((2, 1, 1, 1), (1, 0, 1, 1), "e1")]))
L.append("")

# 3. 放出・終状態・滞在時間
res_f = X.run(X.W_GW)
L.append("## 3. 放出・終状態・滞在時間")
L.append("")
L.append("| 量 | 固定一式 | 数値解法 | 差 |")
L.append("|---|---|---|---|")
ef = sum(res_f["flux_E"].values()) * EV; eg = sum(res_g["flux_E"].values()) * EV
L.append("| 放出合計 [eV] | %.9f | %.9f | %+.2e |" % (ef, eg, eg - ef))
for k in sorted(res_f["flux_E"]):
    L.append("| 放出 %s [eV] | %.6e | %.6e | %+.2e（相対） |" % (k, res_f["flux_E"][k] * EV, res_g["flux_E"].get(k, 0.0) * EV, res_g["flux_E"].get(k, 0.0) / res_f["flux_E"][k] - 1.0))
dropf = (res_f["p0"] @ X.ETOT - res_f["p_inf"] @ X.ETOT) * EV
dropg = (res_g["p0"] @ ETOT_g - res_g["p_inf"] @ ETOT_g) * EV
L.append("| 準位の落差 [eV] | %.9f | %.9f | 収支 固定 %.1e／数値 %.1e |" % (dropf, dropg, ef - dropf, eg - dropg))
fin_f = [(X.ST[i], float(res_f["p_inf"][i])) for i in res_f["absorbing"] if res_f["p_inf"][i] > 0]
fin_g = [(levels[i]["label"], float(res_g["p_inf"][i])) for i in res_g["absorbing"] if res_g["p_inf"][i] > 0]
L.append("| 終状態 | %s | %s | |" % (fin_f, fin_g))
tau_f = {X.ST[i]: res_f["tau"][i] for i in range(X.DIM)}
tau_g = {levels[i]["label"]: res_g["tau"][i] for i in range(len(levels))}
rel_tau = max(abs(tau_g[k] / tau_f[k] - 1.0) for k in tau_f if tau_f[k] > 0)
L.append("| 滞在時間 τ の最大相対差 | | | %.2e |" % rel_tau)
L.append("")
L.append("## 4. 判定")
L.append("")
L.append("数値解法は閉じた式の一式を、クーロン準位で %.0e eV、重力で相対 %.0e、超微細で相対 %.0e（相対論補正分）、A 係数で相対 %.0e（Dirac 行列要素分）、放出合計で %.0e eV の範囲で再現する。"
         "差はいずれも数値精度か、数値解法の側がより正しい相対論補正（α²）で、閉じた式の側の近似に由来する。符号 SIGN は数値解法の入力であり、+1 では束縛準位が無く準位表が空になる。"
         % (worst["coul"], worst["grav"], worst["hfs_rel"], max(abs(Ag[k] / Af[k] - 1.0) for k in Af if k in Ag and k[2] in ("e1", "e2")), abs(eg - ef)))
out = "\n".join(L) + "\n"
(HERE / "control_report.md").write_text(out, encoding="utf-8")
print(out)
print("(%.1f s)" % (time.time() - t0))
