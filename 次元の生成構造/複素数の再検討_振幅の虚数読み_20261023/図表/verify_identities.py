# -*- coding: utf-8 -*-
"""本稿で使う恒等式を数値で確かめ、results/checks.txt と results/checks.json に残す。

記法：w ＝ log u ＝ φ ＋ iθ（φ ＝ 対数振幅、θ ＝ 位相）。複素角は ω ＝ θ − iφ で u ＝ e^{iω}、w ＝ iω。
     cosh(w/2) ＝ cos(ω/2)、sinh(w/2) ＝ i sin(ω/2) なので、w で書けば cosh・sinh、ω で書けば cos・sin。
C1  |cosh(w/2)|² ＋ |sinh(w/2)|² ＝ cosh φ            （＝ |cos(ω/2)|² ＋ |sin(ω/2)|²）
C2  |cosh(w/2)|² − |sinh(w/2)|² ＝ cos θ、φ ＝ 0 で cos²(θ/2) ＋ sin²(θ/2) ＝ 1
C3  |1 ＋ e^w|² ＝ 2e^φ(cosh φ ＋ cos θ) ＝ 4e^φ |cosh(w/2)|²
C4  (1 − u)/(1 ＋ u) ＝ −tanh(w/2)、u ＝ e^w（複素）
C5  |k|² ＝ (cosh φ − cos θ)/(cosh φ ＋ cos θ)。|k| ＝ 1 ⟺ cos θ ＝ 0（φ によらない）
C6  単位四元数：u ＝ e^{nα} で (1＋u)⁻¹(1−u) ＝ −n tan(α/2) ＝ −Im u/(1 ＋ Re u)（修正ロドリゲス・パラメータ）
C7  交換子 e^{ma}e^{nθ} − e^{nθ}e^{ma} ＝ 2 sin a sin θ (m × n)（m ⟂ n、単位純虚）
C8  e^{ma ＋ nθ} ≠ e^{ma}e^{nθ}（BCH の補正。順序の宣言が要る）
C9  ℂ では虚数読みした振幅は位相と合流する：e^{i(a＋θ)} の大きさは 1、a は θ のずれと区別できない
C10 距離 D ＝ cot(α/2) と結合 |k| ＝ tan(α/2) は D ＝ 1/|k|、α ＝ π/2 で両方 1（自己双対点）
乱数は使わない。格子は固定。
"""
from pathlib import Path
import json
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"; RES.mkdir(exist_ok=True)
out = {}; lines = []


def rec(key, ok, detail):
    out[key] = {"ok": bool(ok), "detail": detail}
    lines.append(f"{key}: {'OK' if ok else 'NG'}  {detail}")


th = np.linspace(-np.pi, np.pi, 721); ph = np.linspace(-3, 3, 301)
TH, PH = np.meshgrid(th, ph)
w = PH + 1j*TH                      # log u
om = TH - 1j*PH                     # 複素角 ω（u ＝ e^{iω}）
c2 = np.abs(np.cosh(w/2))**2; s2 = np.abs(np.sinh(w/2))**2
c2o = np.abs(np.cos(om/2))**2; s2o = np.abs(np.sin(om/2))**2
rec("C0_w_vs_omega", np.max(np.abs(c2 - c2o)) < 1e-12 and np.max(np.abs(s2 - s2o)) < 1e-12, f"|cosh(w/2)|²=|cos(ω/2)|² max dev {np.max(np.abs(c2 - c2o)):.2e}; |sinh(w/2)|²=|sin(ω/2)|² max dev {np.max(np.abs(s2 - s2o)):.2e}")
rec("C1_sum_is_cosh", np.max(np.abs(c2 + s2 - np.cosh(PH))) < 1e-10, f"max|c2+s2-cosh φ| = {np.max(np.abs(c2 + s2 - np.cosh(PH))):.2e}")
rec("C2_diff_is_cos", np.max(np.abs(c2 - s2 - np.cos(TH))) < 1e-12, f"max|c2-s2-cos θ| = {np.max(np.abs(c2 - s2 - np.cos(TH))):.2e}; at φ=0 sum-1 max = {np.max(np.abs((c2+s2)[PH==0]-1)):.2e}")
u = np.exp(w)
den = np.abs(1 + u)**2
rec("C3_denominator", (np.max(np.abs(den - 2*np.exp(PH)*(np.cosh(PH) + np.cos(TH)))) < 1e-10) and (np.max(np.abs(den - 4*np.exp(PH)*c2)) < 1e-10),
    f"max|den-2e^φ(cosh φ+cos θ)| = {np.max(np.abs(den - 2*np.exp(PH)*(np.cosh(PH)+np.cos(TH)))):.2e}; max|den-4e^φ|cosh(w/2)|²| = {np.max(np.abs(den - 4*np.exp(PH)*c2)):.2e}")
mask = np.abs(1 + u) > 1e-3            # 極（u ＝ −1、接触）の近くは除く
k = np.where(mask, (1 - u)/np.where(mask, 1 + u, 1), np.nan)
rec("C4_k_is_minus_tanh", np.nanmax(np.abs(k + np.tanh(w/2))) < 1e-9, f"max|k+tanh(w/2)| = {np.nanmax(np.abs(k + np.tanh(w/2))):.2e} (|1+u|>1e-3)")
ak2 = np.abs(k)**2; pred = (np.cosh(PH) - np.cos(TH))/(np.cosh(PH) + np.cos(TH))
rel = np.nanmax(np.abs(ak2 - pred)/(1 + np.abs(pred)))
rec("C5_abs_k_squared", rel < 1e-10, f"max rel|k²-pred| = {rel:.2e}; |k|=1 at θ=±π/2 for φ∈[-3,3]: max dev = {np.nanmax(np.abs(np.abs((1-np.exp(ph+1j*np.pi/2))/(1+np.exp(ph+1j*np.pi/2)))-1)):.2e}")


def qmul(p, q):
    w1, x1, y1, z1 = p; w2, x2, y2, z2 = q
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2, w1*x2 + x1*w2 + y1*z2 - z1*y2, w1*y2 - x1*z2 + y1*w2 + z1*x2, w1*z2 + x1*y2 - y1*x2 + z1*w2])


def qinv(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])/np.dot(q, q)


def qexp(v):
    """純虚ベクトル v（3 成分）の指数 e^{v} ＝ cos|v| ＋ (v/|v|) sin|v|"""
    t = np.linalg.norm(v)
    if t < 1e-15: return np.array([1.0, 0, 0, 0])
    return np.concatenate([[np.cos(t)], np.sin(t)*v/t])


# C6: MRP
n = np.array([0.36, 0.48, 0.80]); n /= np.linalg.norm(n)
errs = []
for al in np.linspace(0.05, 3.0, 60):
    u_q = qexp(n*al)
    kq = qmul(qinv(np.array([1, 0, 0, 0]) + u_q), np.array([1, 0, 0, 0]) - u_q)
    e1 = np.linalg.norm(kq - np.concatenate([[0], -n*np.tan(al/2)]))
    e2 = np.linalg.norm(kq[1:] + u_q[1:]/(1 + u_q[0]))
    errs.append(max(e1, e2))
rec("C6_MRP_identity", max(errs) < 1e-12, f"max err over α∈[0.05,3.0] = {max(errs):.2e}")
# C7: commutator
m = np.array([1.0, 0, 0]); nn = np.array([0, 1.0, 0]); mxn = np.cross(m, nn)
errs = []
for a in np.linspace(0, np.pi, 25):
    for t in np.linspace(0, np.pi, 25):
        comm = qmul(qexp(m*a), qexp(nn*t)) - qmul(qexp(nn*t), qexp(m*a))
        errs.append(np.linalg.norm(comm - np.concatenate([[0], 2*np.sin(a)*np.sin(t)*mxn])))
rec("C7_commutator", max(errs) < 1e-12, f"max err = {max(errs):.2e}")
# C8: BCH
a, t = np.pi/4, np.pi/4
d = np.linalg.norm(qexp(m*a + nn*t) - qmul(qexp(m*a), qexp(nn*t)))
rec("C8_BCH_not_equal", d > 1e-3, f"|e^(ma+nθ) − e^(ma)e^(nθ)| = {d:.4f} at a=θ=π/4 (non-zero: ordering must be declared)")
# C9: in ℂ, imaginary amplitude merges with phase
aa, tt = 0.7, 0.4
z = np.exp(1j*(aa + tt))
rec("C9_C_merge", abs(abs(z) - 1) < 1e-15, f"|e^{{i(a+θ)}}| = {abs(z):.15f}; a and θ enter only as a+θ")
# C10: D = 1/|k|, both 1 at α = π/2
al = np.pi/2
rec("C10_selfdual", abs(np.tan(al/2) - 1) < 1e-15 and abs(1/np.tan(al/2) - 1) < 1e-15, f"tan(π/4) = {np.tan(al/2):.15f}, cot(π/4) = {1/np.tan(al/2):.15f}")

allok = all(v["ok"] for v in out.values())
lines.append(f"ALL: {'OK' if allok else 'NG'} ({sum(v['ok'] for v in out.values())}/{len(out)})")
(RES / "checks.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
(RES / "checks.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
print("\n".join(lines))
