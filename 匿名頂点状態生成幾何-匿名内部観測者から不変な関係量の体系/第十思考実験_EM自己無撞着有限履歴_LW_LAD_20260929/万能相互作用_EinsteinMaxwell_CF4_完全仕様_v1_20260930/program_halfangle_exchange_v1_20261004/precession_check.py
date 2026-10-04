#!/usr/bin/env python3
"""対照：近点移動の一次摂動（∂⟨δε⟩/∂L）を、v7 §11.1 の表（円軌道の近点移動 (γ_cl − 1)/α²）と比べる。
ε の項の部分集合ごとに、円軌道 a = ξ₁、L = √(κa) で一周あたりの近点移動を 2πα² で割る。純粋な算術。"""
import math
import halfangle_exchange_node_map as H

ALPHA2 = H.ALPHA**2
a = H.XI1
L = math.sqrt(H.KAPPA * a)
T = 2 * math.pi * a**1.5 / math.sqrt(H.KAPPA)

def gamma_minus_1_over_alpha2(mono_subset, sa=1, sb=-1):
    saved = H.MONO
    H.MONO = mono_subset
    try:
        rate = H.precession_rate(a, 0.0, L, sa, sb)
    finally:
        H.MONO = saved
    return rate * T / (2 * math.pi) / ALPHA2

full = H.MONO
kep = [t for t in full if (t[0], t[1], t[2]) in {(0, 2, 0), (1, 0, 0), (2, 0, 2)}]          # T1〜T3（T1 の m² 部は p=2,b=2）
sr = kep + [t for t in full if (t[0], t[1], t[2]) == (0, 4, 0)] + [t for t in full if (t[0], t[1], t[2]) in {(2, 2, 2), (4, 0, 4)}]
# T4 = (3ν−1)𝐏⁴/8 は単項式では (0,4,0),(2,2,2),(4,0,4) の三つ
kn = sr + [t for t in full if t[0] == 3 and t[1] == 0 and t[2] == 0 and t[3] == 0 and t[4] == 0]   # T19〜T22 と T12 は p=3 の m,s なし項
def show(name, subset):
    print(f"{name:28s} (γ_cl − 1)/α² = {gamma_minus_1_over_alpha2(subset):.4f}")
show("Newton+Coulomb (expect 0)", kep)
show("+SR p^4 (v7: 0.4993)", sr)
show("+KN quadrupole (v7: 0.8741)", kn)
show("all 22, up/down (v7: 0.3748)", full)
print(f"all 22, spin 0 (v7: 0.8750): {gamma_minus_1_over_alpha2(full, 0, 0):.4f}")
print(f"T = {T:.6e}, L = {L:.9f}, a = {a}")
