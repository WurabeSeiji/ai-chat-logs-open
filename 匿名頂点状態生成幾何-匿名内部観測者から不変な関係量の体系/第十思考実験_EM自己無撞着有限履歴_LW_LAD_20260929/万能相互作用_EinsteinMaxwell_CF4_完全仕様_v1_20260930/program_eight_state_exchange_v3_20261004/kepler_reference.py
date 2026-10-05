#!/usr/bin/env python3
"""対照：--em kepler --f 0 の走行を Kepler の閉じた式と比べる（v2）。
節点 n の離心近点角は E_K = 2nθ₀（近点から）。ξ_K = a(1 − e cos E_K)。
走行の ξ、状態から読んだ周期内の位置 E_cycle と時刻 t_cycle を、閉じた式と比べる。不変量 a、L、ε_K の揺れも出す。"""
import csv, gzip, math, sys

ALPHA = 7.2973525643e-3; M_A = 7.0e-4; RATIO = 1836.15267343; M_B = M_A * RATIO
GM = (M_A + M_B) / M_B
KC = 2 * ALPHA * GM; KG = KC / (9.0 / (M_A * M_B)); KAPPA = KC + KG
TH0 = 23 * math.pi / 124
TWO_PI = 2 * math.pi

p = sys.argv[1]
op = gzip.open if p.endswith(".gz") else open
with op(p, "rt") as fh:
    rows = list(csv.DictReader(fh))
r0 = rows[0]
a = float(r0["a_orb"]); e = float(r0["ecc"]); om = float(r0["omega"])
n_mean = math.sqrt(KAPPA / a**3)
print(f"a = {a:.9f}, e = {e:.9e}, omega = {om:.9e}, 2omega/a = {2*om/a:.9e}, sqrt(kappa/a^3) = {n_mean:.9e}")
print(f"theta0 = {TH0:.9f} (R = sin^2 = {math.sin(TH0)**2:.6f}), nodes per orbit = {math.pi/TH0:.6f}")
mx_xi = mx_E = mx_t = 0.0
for r in rows:
    n = int(r["node"]); EK = (2 * n * TH0) % TWO_PI
    mx_xi = max(mx_xi, abs(float(r["xi"]) - a * (1 - e * math.cos(EK))))
    if e > 1e-4:
        dE = (float(r["E_cycle"]) - EK + math.pi) % TWO_PI - math.pi
        mx_E = max(mx_E, abs(dE))
        tK = (EK - e * math.sin(EK)) / n_mean
        per = TWO_PI / n_mean
        dt = abs(float(r["t_cycle"]) - tK); dt = min(dt, per - dt)
        mx_t = max(mx_t, dt)
def spread(key):
    v = [float(r[key]) for r in rows]; return max(v) - min(v), v[0]
sN, N0 = spread("a_orb"); sL, L0 = spread("L"); sE, E0 = spread("eps_K")
print(f"max|xi - xi_K| = {mx_xi:.3e}  (xi range {min(float(r['xi']) for r in rows):.6f}..{max(float(r['xi']) for r in rows):.6f})")
if e > 1e-4:
    print(f"max|E_cycle - E_K| = {mx_E:.3e} rad;  max|t_cycle - t_K| = {mx_t:.3e}  (period {2*math.pi/n_mean:.6e})")
else:
    print("E_cycle / t_cycle: e < 1e-4 なので周期内の位置は定義できない（円）")
print(f"spread a = {sN:.3e} (of {N0:.6f}); spread L = {sL:.3e} (of {L0:.9f}); spread eps_K = {sE:.3e} (of {E0:.9e})")
# 124 節点ごとに ξ が戻るか（関係の再帰。状態から直接）
rec = [abs(float(rows[k]["xi"]) - float(rows[0]["xi"])) for k in range(124, len(rows), 124)]
if rec:
    print(f"124-node recurrence of xi: max |xi_n - xi_0| over {len(rec)} checks = {max(rec):.3e}")
