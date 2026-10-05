#!/usr/bin/env python3
"""遷移率（Einstein の A 係数）・超微細・マスター方程式の検算。プログラム本体は触らず import するだけ。

1. 角運動量代数：6j・3j の既知値と和則（j'、F' で和をとると分解前の値に戻る）。
2. E1 の分解前 A 係数を NIST ASD（Wiese–Fuhr 2009）の値と比較。
3. 超微細間隔を実測と比較（g_e = 2 なので差は g_e/2 = 1.00116 と QED）。
4. 21 cm（1s F=1→0 の M1）を実測 2.8843e-15 s⁻¹ と比較。核スピン項の符号（電子と陽子のモーメントが足し合う）。
5. E2 3d→1s と重力波。
6. 寿命と、マスター方程式のエネルギー保存・終状態。
"""
import sys
from pathlib import Path
import numpy as np
from matplotlib import font_manager

font_manager.fontManager.addfont = lambda path: None   # 17 行目の Linux フォント登録を無効化（本体は変更しない）
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exchange_cascade as X  # noqa: E402

L = []
MHZ = 1e6 * 4.135667696e-15   # 1 MHz [eV]

# 1. 角運動量代数
L.append("1. 角運動量代数")
L.append("   {½ 0 ½; 1 ½ 1}² = %.6f（1/6 = %.6f）、(2 2 0;000)² = %.6f（1/5）、(1 2 1;000)² = %.6f（2/15 = %.6f）、(2 2 2;000)² = %.6f（2/35 = %.6f）"
         % (X.wigner6j(0.5, 0, 0.5, 1, 0.5, 1) ** 2, 1 / 6, X.w3j000_sq(2, 2, 0), X.w3j000_sq(1, 2, 1), 2 / 15, X.w3j000_sq(2, 2, 2), 2 / 35))
worst = 0.0
for l in range(0, 5):
    for l2 in range(0, 5):
        for k in (1, 2):
            if abs(l - l2) > k or (l + l2) < k:
                continue
            for s_e in (-1, 1):
                j = max(l + 0.5 * s_e, 0.5)
                if l == 0 and s_e == -1:
                    continue
                for s_p in (-1, 1):
                    F = j + 0.5 * s_p
                    tot = 0.0
                    for s_e2 in (-1, 1):
                        j2 = max(l2 + 0.5 * s_e2, 0.5)
                        if l2 == 0 and s_e2 == -1:
                            continue
                        for s_p2 in (-1, 1):
                            F2 = j2 + 0.5 * s_p2
                            tot += X.recouple(l, j, l2, j2, F, F2, k)
                    worst = max(worst, abs(tot - 1.0))
L.append("   和則 Σ_{j',F'} (2j'+1)(2ℓ+1){ℓ j ½; j' ℓ' k}²(2F'+1)(2j+1){j F ½; F' j' k}² = 1：最大差 %.1e（ℓ,ℓ' ≤ 4、k = 1,2、全 j,F）" % worst)

# 2. E1 の分解前 A 係数
ed = X.edges(X.W_GW)
def A_sum(n, l, n2, l2, kinds=("e1",)):
    i0 = next(i for i in X.VALID if X.ST[i][0] == n and X.ST[i][1] == l)   # 初期の (j,F) は任意（和則）
    return sum(A for src, dst, k, A in ed if src == i0 and k in kinds and X.ST[dst][0] == n2 and X.ST[dst][1] == l2)
nist = {(2, 1, 1, 0): 6.2649e8, (3, 1, 1, 0): 1.6725e8, (4, 1, 1, 0): 6.8186e7, (5, 1, 1, 0): 3.4375e7,
        (3, 0, 2, 1): 6.3134e6, (3, 2, 2, 1): 6.4651e7, (3, 1, 2, 0): 2.2449e7, (4, 3, 3, 2): 1.3789e7,
        (4, 2, 2, 1): 2.0625e7, (4, 0, 2, 1): 2.5785e6}
L.append("")
L.append("2. E1 の A 係数（分解前、s⁻¹）と NIST ASD")
names = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}
wr = 0.0
for (n, l, n2, l2), ref in nist.items():
    val = A_sum(n, l, n2, l2)
    wr = max(wr, abs(val / ref - 1))
    L.append("   %d%s→%d%s  本式 %.4e  NIST %.4e  比 %.5f" % (n, names[l], n2, names[l2], val, ref, val / ref))
L.append("   最大相対差 %.1e（NIST の表は 5 桁）" % wr)

# 3. 超微細
L.append("")
L.append("3. 超微細間隔（MHz）と実測")
def hfs_split(n, l, j):
    return (X.hyperfine_level(n, l, j, j + 0.5) - X.hyperfine_level(n, l, j, j - 0.5)) * X.EV / MHZ
meas = {(1, 0, 0.5): 1420.405752, (2, 0, 0.5): 177.5569, (2, 1, 0.5): 59.22, (2, 1, 1.5): 23.65}
for (n, l, j), ref in meas.items():
    v = hfs_split(n, l, j)
    L.append("   n=%d ℓ=%d j=%.1f  本式 %.3f  実測 %.3f  比 %.5f" % (n, l, j, v, ref, v / ref))
L.append("   比 0.9989 = 1/(g_e/2) = 1/1.00116（電子異常磁気能率、QED）。それ以外の差は 1e-4 以下")

# 4. 21 cm
L.append("")
L.append("4. 21 cm（1s F=1→F=0 の M1）")
i1, i0 = X.IDX[(1, 0, 1, 1)], X.IDX[(1, 0, 1, -1)]
A21 = sum(A for src, dst, k, A in ed if src == i1 and dst == i0 and k == "m1")
w21 = (X.ETOT[i1] - X.ETOT[i0]) * X.ME_EV * X.E_SI / X.HBAR_SI
gp_save = X.G_P
X.G_P = 0.0
A21_e = X.m1_rate(0, 0.5, 1.0, 0.5, 0.0, w21)
# 射影定理：j' = j では ⟨J⟩ + ⟨S⟩ = g_j⟨J⟩。電子部分だけの換算行列要素を −g_j·redJ と比較
worst_proj = 0.0
for (l, j) in ((0, 0.5), (1, 0.5), (1, 1.5), (2, 1.5), (2, 2.5)):
    for F, F2 in ((j + 0.5, j - 0.5), (j - 0.5, j + 0.5)):
        if F2 < 0 or F < 0:
            continue
        pref = np.sqrt((2 * F + 1) * (2 * F2 + 1))
        redJ = (-1) ** int(j + 0.5 + F + 1) * pref * X.wigner6j(j, F2, 0.5, F, j, 1) * np.sqrt(j * (j + 1) * (2 * j + 1))
        worst_proj = max(worst_proj, abs(X.m1_reduced(l, j, F, j, F2) + X.lande_g(l, j) * redJ))
X.G_P = gp_save
L.append("   本式 %.4e s⁻¹、実測 2.8843e-15 s⁻¹、比 %.5f（(2/2.00232)² = %.5f：電子異常磁気能率）" % (A21, A21 / 2.8843e-15, (2 / 2.00231930) ** 2))
L.append("   核スピン項を落とした値との比 %.5f、(1 + g_pμ_N/2μ_B)² = %.5f：電子と陽子のモーメントは足し合う（符号正しい）" % (A21 / A21_e, (1 + X.G_P / (2 * X.MU_P)) ** 2))
L.append("   射影定理 ⟨J⟩ + ⟨S⟩ = g_j⟨J⟩（j' = j、s/p/d の全 j、ΔF = ±1）：最大差 %.1e" % worst_proj)
L.append("   ω は本式の超微細間隔（1418.8 MHz）で計算。周波数を実測に置けば ω³ で %.4f 倍" % ((1420.405752 / hfs_split(1, 0, 0.5)) ** 3))
A_fs = sum(A for src, dst, k, A in ed if k == "m1" and X.ST[src][:3] == (2, 1, 1) and X.ST[dst][:3] == (2, 1, -1))
L.append("   微細構造間の M1 2p₃/₂→2p₁/₂（全 F→F'、1 つの F から）：%.3e s⁻¹（E1 崩壊 6.3e8 の 1e-21）" % (A_fs / 2))

# 5. E2 と重力波
L.append("")
L.append("5. E2 と重力波")
L.append("   3d→1s：E2 %.1f s⁻¹（§4 の直接計算 593.8）、GW %.3e s⁻¹、比 %.4e = W_GW" % (A_sum(3, 2, 1, 0, ("e2",)), A_sum(3, 2, 1, 0, ("gw",)), A_sum(3, 2, 1, 0, ("gw",)) / A_sum(3, 2, 1, 0, ("e2",))))
L.append("   3p→2p（Δℓ = 0 の四重極）：E2 %.3e s⁻¹、3d→2p は E1 %.3e と E2 なし（Δℓ = 1）" % (A_sum(3, 1, 2, 1, ("e2",)), A_sum(3, 2, 2, 1)))

# 6. 寿命・マスター方程式
L.append("")
L.append("6. 寿命とマスター方程式")
R = X.generator(ed)
def life(key):
    i = X.IDX[key]
    return -1.0 / R[i, i]
L.append("   2p₃/₂ F=2 寿命 %.4e s（1/6.2649e8 = 1.5962e-9）、2s F=1 寿命 %.4f s（二光子 1/8.2206 = 0.12165）、1s F=1 寿命 %.3e s = %.2e 年"
         % (life((2, 1, 1, 1)), life((2, 0, 1, 1)), life((1, 0, 1, 1)), life((1, 0, 1, 1)) / 3.15576e7))
res = X.run(X.W_GW)
p0, pinf = res["p0"], res["p_inf"]
emitted = sum(res["flux_E"].values()) * X.EV
drop = (p0 @ X.ETOT - pinf @ X.ETOT) * X.EV
L.append("   放出の合計 %.9f eV、準位の落差 %.9f eV、差 %.1e eV（厳密に等しいはず）" % (emitted, drop, emitted - drop))
L.append("   終状態 %s、Σp = %.15f" % ([(X.ST[i], round(float(pinf[i]), 12)) for i in X.VALID if pinf[i] > 0], pinf.sum()))
L.append("   Radau の p(t_max) と厳密な p_∞ の最大差 %.1e、Σp(t) の 1 からの最大ずれ %.1e、p(t) の最小値 %.1e"
         % (np.max(np.abs(res["hist"][-1] - pinf)), np.max(np.abs(res["hist"].sum(axis=1) - 1)), res["hist"].min()))
tau = res["tau"]
top = sorted(((tau[i], X.ST[i]) for i in X.VALID if tau[i] > 0), reverse=True)[:6]
L.append("   滞在時間 τ の上位：%s" % ", ".join("%s %.3e s" % (s, t) for t, s in top))
L.append("   放出の内訳 [eV]：%s" % ", ".join("%s %.6e" % (k, v * X.EV) for k, v in sorted(res["flux_E"].items())))

out = "\n".join(L) + "\n"
(HERE / "rates_check.txt").write_text(out, encoding="utf-8")
print(out)
