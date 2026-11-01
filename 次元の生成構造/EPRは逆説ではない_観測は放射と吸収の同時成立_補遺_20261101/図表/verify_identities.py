# -*- coding: utf-8 -*-
"""恒等式の数値確認（固定格子、乱数なし）。結果を results/checks.txt に書く。
C1：q_B = −q_A から、読出し平面への像 x と その時間変化率 v の和が全時刻でゼロ
C2：片側の符号（軸 n_A に対する cos の符号）のゲージ平均が 1/2
C3：交差項 Sc[q_A q̄_B] = −|q_A|² が全体の位相 φ と時刻 t によらない
C4：A だけに作用する操作（A の軸への射影・符号の確定）の前後で、B の符号分布が変わらない
記法：q = A e^{nθ}、読出し平面は span(1, n)。位置 x = (Re q, Im_n q)、速度 v = dx/dt。"""
import numpy as np, pathlib
out = pathlib.Path(__file__).resolve().parent / "results"; out.mkdir(exist_ok=True)
lines = []
def rec(s): lines.append(s); print(s)
A, w = 1.3, 0.7                      # 振幅、角速度（任意の固定値）
t = np.linspace(0.0, 20.0, 2001)
phi = np.linspace(0.0, 2 * np.pi, 3601)[:-1]  # ゲージ（全体の位相）の格子
# C1
theta0 = 0.4
qA = A * np.exp(1j * (w * t + theta0)); qB = -qA
xA = np.c_[qA.real, qA.imag]; xB = np.c_[qB.real, qB.imag]
vA = np.gradient(xA, t, axis=0); vB = np.gradient(xB, t, axis=0)
rec("C1 max|x_A+x_B| = %.3e, max|v_A+v_B| = %.3e  (expected 0)" % (np.abs(xA + xB).max(), np.abs(vA + vB).max()))
rec("C1 also: x_A−x_B は一定の大きさ %.6f = 2A, p_A+p_B=0 と同時に成立" % np.hypot(*(xA - xB).T).mean())
# C2
nA = 0.0                                      # 軸 n_A の向き（読出し平面内の角）
sA = np.sign(np.cos(phi - nA)); pA_plus = (sA > 0).mean()
rec("C2 P_A(+) over gauge = %.6f  (expected 0.5)" % pA_plus)
# C3
cross = np.array([[(qa * np.conj(-qa)).real for qa in (A * np.exp(1j * (w * tt + ph)) for tt in t[::200])] for ph in phi[::600]])
rec("C3 Sc[q_A conj(q_B)] : min %.6f max %.6f  (expected −A² = %.6f, constant)" % (cross.min(), cross.max(), -A * A))
# C4：A を軸 n_A に落とす（破壊的読出し）。B の状態は q_B のまま。B の軸 n_B に対する符号分布の前後比較
nB = 0.9
sB_before = (np.sign(np.cos(phi + np.pi - nB)) > 0).mean()
qA_read = np.sign(np.cos(phi - nA)) * A * np.exp(1j * nA)   # A は ± の一方に落ちる
sB_after = (np.sign(np.cos(phi + np.pi - nB)) > 0).mean()    # q_B は変わらない
corr = np.mean(np.sign(np.cos(phi - nA)) * np.sign(np.cos(phi + np.pi - nB)))
rec("C4 P_B(+) before = %.6f, after = %.6f  (expected equal, 0.5); sign correlation = %.4f (= −(1−2|n_A−n_B|/π) = %.4f for signs)" % (sB_before, sB_after, corr, -(1 - 2 * abs(nA - nB) / np.pi)))
(out / "checks.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
