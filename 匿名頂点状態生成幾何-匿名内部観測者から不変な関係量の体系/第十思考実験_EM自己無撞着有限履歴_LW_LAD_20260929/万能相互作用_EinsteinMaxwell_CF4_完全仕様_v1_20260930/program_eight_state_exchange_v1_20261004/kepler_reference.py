#!/usr/bin/env python3
"""対照：--em kepler --f 0 の走行を Kepler の閉じた式と比べる。
節点 n で離心近点角 E = 2nθ₀（θ₀ = 23π/124）。ξ_K = a(1 − e cos E)、t_K = (E − e sin E)/n_mean、n_mean = √(κ/a³)。
走行の ξ、T との最大差、不変量 a・L・ε_K の揺れ、124 節点ごとの関係 z の戻りを出す。"""
import csv, gzip, math, sys

ALPHA = 7.2973525643e-3; M_A = 7.0e-4; RATIO = 1836.15267343; M_B = M_A * RATIO
GM = (M_A + M_B) / M_B
KC = 2 * ALPHA * GM; KG = KC / (9.0 / (M_A * M_B)); KAPPA = KC + KG
TH0 = 23 * math.pi / 124

p = sys.argv[1]
op = gzip.open if p.endswith(".gz") else open
with op(p, "rt") as fh:
    rows = list(csv.DictReader(fh))
r0 = rows[0]
a = float(r0["a_orb"]); e = float(r0["ecc"]); om = float(r0["omega"]); dtau = float(r0["dtau"])
n_mean = math.sqrt(KAPPA / a**3)
print(f"a = {a:.9f}, e = {e:.9e}, omega = {om:.9e}, dtau = {dtau:.6f}, 2omega/a = {2*om/a:.9e}, sqrt(kappa/a^3) = {n_mean:.9e}")
print(f"theta0 = {TH0:.9f} (R = sin^2 = {math.sin(TH0)**2:.6f}), nodes per orbit = {math.pi/TH0:.6f}, 124 nodes = {124*TH0/math.pi:.1f} orbits")
mx_xi = mx_t = 0.0
for r in rows:
    n = int(r["node"]); E = 2 * n * TH0
    mx_xi = max(mx_xi, abs(float(r["xi"]) - a * (1 - e * math.cos(E))))
    mx_t = max(mx_t, abs(float(r["T"]) - (E - e * math.sin(E)) / n_mean))
def spread(key):
    v = [float(r[key]) for r in rows]; return max(v) - min(v), v[0]
sN, N0 = spread("a_orb"); sL, L0 = spread("L"); sE, E0 = spread("eps_K")
print(f"max|xi - xi_K| = {mx_xi:.3e}  (xi range {min(float(r['xi']) for r in rows):.6f}..{max(float(r['xi']) for r in rows):.6f})")
print(f"max|T - t_K|   = {mx_t:.3e}  (T final {float(rows[-1]['T']):.6e}; period 2pi/n = {2*math.pi/n_mean:.6e})")
print(f"spread a = {sN:.3e} (of {N0:.6f}); spread L = {sL:.3e} (of {L0:.9f}); spread eps_K = {sE:.3e} (of {E0:.9e})")
rec = [(int(r["node"]), float(r["rec124_dz"])) for r in rows if r["rec124_dz"] not in ("nan", "")]
if rec:
    print(f"124-node recurrence of z: max relative |z_n - z_0| over {len(rec)} checks = {max(d for _, d in rec):.3e}")
