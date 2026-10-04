#!/usr/bin/env python3
"""対照：--em kepler --f 0 の走行を Kepler の閉じた式と比べる。
初期行の a、e から、E = 2ωs（離心近点角）、ξ_K = a(1 − e cos E)、t_K = (E − e sin E)/n、n = √(κ/a³) を作り、
走行の ξ、t との最大差を出す。不変量 N_p + N_q、L、ε_K の揺れも出す。"""
import csv, math, sys

KC = 2 * 7.2973525643e-3 * (1 + 7.0e-4 / (7.0e-4 * 1836.15267343))
M_A = 7.0e-4; M_B = M_A * 1836.15267343
KG = KC / (9.0 / (M_A * M_B))
KAPPA = KC + KG

p = sys.argv[1]
with open(p) as fh:
    rows = list(csv.DictReader(fh))
r0 = rows[0]
a = float(r0["a_orb"]); e = float(r0["ecc"]); om = float(r0["omega"])
n_mean = math.sqrt(KAPPA / a**3)
print(f"a = {a:.9f}, e = {e:.9e}, omega = {om:.9e}, 2omega/a = {2*om/a:.9e}, sqrt(kappa/a^3) = {n_mean:.9e}")
mx_xi = 0.0; mx_t = 0.0
for r in rows:
    s = float(r["s"]); E = 2 * om * s
    xi_K = a * (1 - e * math.cos(E))
    t_K = (E - e * math.sin(E)) / n_mean
    mx_xi = max(mx_xi, abs(float(r["xi"]) - xi_K))
    mx_t = max(mx_t, abs(float(r["t"]) - t_K))
def spread(key):
    v = [float(r[key]) for r in rows]
    return max(v) - min(v), v[0]
sN, N0 = spread("a_orb"); sL, L0 = spread("L"); sE, E0 = spread("eps_K")
print(f"max|xi - xi_K| = {mx_xi:.3e}  (xi range {min(float(r['xi']) for r in rows):.6f}..{max(float(r['xi']) for r in rows):.6f})")
print(f"max|t - t_K|   = {mx_t:.3e}  (t final {float(rows[-1]['t']):.6e}; period 2pi/n = {2*math.pi/n_mean:.6e})")
print(f"spread N_p+N_q = {sN:.3e} (of {N0:.6f}); spread L = {sL:.3e} (of {L0:.9f}); spread eps_K = {sE:.3e} (of {E0:.9e})")
th = float(r0["theta"]); print(f"theta = {th:.9f}, 2theta/2pi = {th/math.pi:.9f}, nodes per orbit = {math.pi/th:.4f}, dt(first node) = {float(r0['dt_node']):.4f}")
