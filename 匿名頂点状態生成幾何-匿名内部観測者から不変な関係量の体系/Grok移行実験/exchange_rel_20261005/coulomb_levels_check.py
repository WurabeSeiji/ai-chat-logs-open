#!/usr/bin/env python3
"""クーロン準位（parts() の orb_c）の検算。プログラム本体は触らず、import して値を読むだけ。

比較対象：
  Bohr      元の 49 行目の式 −μ_red·81/(2ħ²n²)
  Dirac     電子だけ相対論 μ_red(f−1)
  +反跳     −μ_red²(f−1)²/(2(m+M))
  +Breit    α⁴μ_red³/(2n³M²)[1/(j+½) − 1/(ℓ+½)]（ℓ ≥ 1）
  合計      プログラムの ORB_C × EV
外部参照値（CODATA 2022 / 実測）：
  1S 電離エネルギー 13.598434599702 eV（実測）、1S Lamb シフト 8172.9 MHz（QED、本式の外）
  2P3/2 − 2P1/2 = 10969.04 MHz（実測、Hagley–Pipkin 1994 + Lamb）
  2S1/2 − 2P1/2 = 1057.845 MHz（Lamb、QED、本式の外）
"""
import sys
from pathlib import Path
import numpy as np
from matplotlib import font_manager

font_manager.fontManager.addfont = lambda path: None   # 17 行目の Linux フォント登録を無効化（本体は変更しない）
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exchange_cascade as X  # noqa: E402

H_EV_S = 4.135667696e-15      # Planck 定数 [eV·s]
MHZ = 1e6 * H_EV_S            # 1 MHz を eV に
A = X.ALPHA
mu = X.MU_P / (1.0 + X.MU_P)
m_tot = 1.0 + X.MU_P


def terms(n, ell, s_e):
    j = max(ell + 0.5 * s_e, 0.5)
    kappa = j + 0.5
    delta = kappa - np.sqrt(kappa**2 - A**2)
    f = (1.0 + (A / (n - delta))**2) ** -0.5
    bohr = -mu * 81.0 / (2.0 * X.HBAR**2 * n**2)
    dirac = mu * (f - 1.0)
    recoil = -mu**2 * (f - 1.0)**2 / (2.0 * m_tot)
    breit = (A**4 * mu**3 / (2.0 * n**3 * X.MU_P**2) * (1.0 / (j + 0.5) - 1.0 / (ell + 0.5))) if ell > 0 else 0.0
    sommerfeld = (-mu * A**2 / (2 * n**2) - mu * A**4 / (2 * n**3) * (1 / (j + 0.5) - 3 / (4 * n) + mu / (4 * n * m_tot))
                  + breit)
    return j, bohr, dirac, recoil, breit, sommerfeld


lines = []
lines.append("単位 eV。EV = m_e c² = %.5f eV（釘付けなし）" % X.EV)
lines.append("")
lines.append("準位 (n,ℓ,s_e→j)        Bohr(元)        Dirac        反跳(2項)      Breit(3項)     合計=ORB_C×EV    展開式との差")
worst_prog = 0.0
worst_exp = 0.0
for n in range(1, 4):
    for ell in range(n):
        for s_e in ((1,) if ell == 0 else (-1, 1)):
            j, bohr, dirac, recoil, breit, somm = terms(n, ell, s_e)
            total = dirac + recoil + breit
            prog = X.ORB_C[X.IDX[(n, ell, s_e, -1)]]
            worst_prog = max(worst_prog, abs(prog - total))
            worst_exp = max(worst_exp, abs(total - somm) / abs(total))
            lines.append("(%d,%d,%+d→%.1f)  %14.9f %14.9f %+14.3e %+14.3e %16.9f   %+.1e"
                         % (n, ell, s_e, j, bohr * X.EV, dirac * X.EV, recoil * X.EV, breit * X.EV, prog * X.EV,
                            (total - somm) * X.EV))
lines.append("")
lines.append("プログラムの ORB_C と本検算の合計の最大差 = %.1e（m_e 単位）" % worst_prog)
lines.append("閉じた式と α⁴ 展開式の最大相対差 = %.1e（次の項は準位に対して相対 O(α⁴)。1s では α⁴/8 = %.1e）" % (worst_exp, A**4 / 8))
lines.append("")

e1s = -X.ORB_C[X.IDX[(1, 0, 1, -1)]] * X.EV
e1s_meas = 13.598434599702
lamb_1s = 8172.9 * MHZ
lines.append("1S 電離エネルギー：本式 %.7f eV、実測 %.7f eV、差 %.3e eV = %.1f MHz" % (e1s, e1s_meas, e1s - e1s_meas, (e1s - e1s_meas) / MHZ))
lines.append("   1S Lamb シフト（QED、本式の外）= %.3e eV = 8172.9 MHz。残差との差 %.1f MHz" % (lamb_1s, (e1s - e1s_meas - lamb_1s) / MHZ))
lines.append("   元の Bohr 釘付け 13.605693 eV は実測と %.3e eV（%.0f MHz）ずれていた" % (13.605693 - e1s_meas, (13.605693 - e1s_meas) / MHZ))
lines.append("")
fs = (X.ORB_C[X.IDX[(2, 1, 1, -1)]] - X.ORB_C[X.IDX[(2, 1, -1, -1)]]) * X.EV
_, _, d32, _, _, _ = terms(2, 1, 1)
_, _, d12, _, _, _ = terms(2, 1, -1)
lines.append("2P3/2 − 2P1/2：本式 %.2f MHz（Dirac 単体 %.2f MHz）、実測 10969.04 MHz、差 %.1f MHz（電子異常磁気能率などの QED）"
             % (fs / MHZ, (d32 - d12) * X.EV / MHZ, fs / MHZ - 10969.04))
ls = (X.ORB_C[X.IDX[(2, 0, 1, -1)]] - X.ORB_C[X.IDX[(2, 1, -1, -1)]]) * X.EV
lines.append("2S1/2 − 2P1/2：本式 %.2f MHz（Dirac では縮退）、実測 1057.845 MHz（Lamb、QED）" % (ls / MHZ))
lines.append("")
lines.append("1s の内訳：Dirac シフト（相対 α²/4）= %.3e eV、反跳 = %.3e eV、Breit = 0（ℓ = 0）"
             % ((terms(1, 0, -1)[2] - (-mu * A**2 / 2)) * X.EV, terms(1, 0, -1)[3] * X.EV))
lines.append("3項（Breit）の最大値 = %.2e eV（倍精度での ORB_C の分解能 ≈ %.0e eV）"
             % (max(abs(terms(n, l, s)[4]) for n in range(1, 6) for l in range(n) for s in (-1, 1)) * X.EV,
                2.7e-5 * 2.2e-16 * X.EV))

out = "\n".join(lines) + "\n"
(HERE / "coulomb_levels_check.txt").write_text(out, encoding="utf-8")
print(out)
