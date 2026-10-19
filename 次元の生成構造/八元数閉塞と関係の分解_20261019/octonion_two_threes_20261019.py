#!/usr/bin/env python3
"""
二つの「三」― 閉塞の三相と内部の Z₃（2026-10-19）

`octonion_closure_check_20261019.py` と同じフォルダに置いて実行する。

問い：e/3 の三相は、N = 3 の閉塞配置 (1, ω, ω²) のことか、それとも内部半分の
Q₁Q₂Q₃ のことか。この二つは同じものか。

計算で出た答え：別のもので、thirds を与えるのは後者。正確には、観測者が内部の時計軸
ℓ を固定したときに残る SU(3) = stab(ℓ) の**中心 Z₃**が thirds を与える。

  1. stab(ℓ) ⊂ G₂ は 8 次元（su(3)）。その Cartan は 2 次元。
  2. 7 本の虚方向の重みは、ゼロ 1 本（= ℓ）と、正六角形の 6 本。
     一つおきの三つ組が 120° 間隔で和ゼロ（3 と 3̄）。
  3. 三つの重みの位相がすべて 2π/3 になる torus の点 Z を解くと、Z は
     G₂ の元、stab(ℓ) の中心、Z³ = I、固有値は ±120° が 3 ペアと 0° が 2 本。
  4. Z は 1 と ℓ を固定し、三つの平面 (v₁,Q₁)、(v₂,Q₂)、(v₃,Q₃) を 120° 回す。
     **空間と内部を混ぜる。** 内部だけの回転ではない。
  5. Z は閉塞解を閉塞解に写す。固定される配置は span(1, ℓ) の中のものだけ。
  6. どの配置でも S + ZS + Z²S = 3 ×（span(1, ℓ) への射影）。
     一つの像だけでは 3 の成分が残り、三つそろうと消える（色の中和）。
  7. 三相 (1, ω, ω²) を内部軸 span(1, ℓ) に置くと Z で不変（トライアリティ 0）。
     空間軸 span(1, e₁) に置くと Z で動き、三つの像を持つ。

再現性：固定シード。出力は標準出力のみ。
"""
import os
import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "octonion_closure_check_20261019.py")).read().split("results = {")[0])

def nullbasis(A, tol=1e-9):
    u, s, vt = np.linalg.svd(A)
    return vt[int(np.sum(s > tol * s[0])):]

# --- g₂ と stab(ℓ) = su(3) -------------------------------------------------
rows = []
for i in range(8):
    for j in range(8):
        prod = omul(E[i], E[j])
        for c in range(8):
            row = np.zeros(64)
            for m in range(8): row[c*8 + m] += prod[m]
            for m in range(8): row[m*8 + i] -= omul(E[m], E[j])[c]
            for m in range(8): row[m*8 + j] -= omul(E[i], E[m])[c]
            rows.append(row)
M = np.array(rows)
g2 = [v.reshape(8, 8) for v in nullbasis(M)]
fix_l = np.array([np.eye(64)[c*8 + 4] for c in range(8)])          # D(e4) = 0
su3 = [v.reshape(8, 8) for v in nullbasis(np.vstack([M, fix_l]))]
print(f"dim g₂ = {len(g2)}, dim stab(ℓ) = {len(su3)}  （ℓ = e4）")

# --- Cartan と重み ---------------------------------------------------------
rng = np.random.default_rng(11)
X = sum(c * D for c, D in zip(rng.normal(size=len(su3)), su3))
ns = nullbasis(np.array([(X @ D - D @ X).reshape(-1) for D in su3]).T)
H = [sum(c * D for c, D in zip(v, su3)) for v in ns]
print(f"Cartan の次元 = {len(H)}, [H₁,H₂] = {np.linalg.norm(H[0]@H[1]-H[1]@H[0]):.1e}")

w1, V = np.linalg.eig(H[0]); H2d = np.linalg.inv(V) @ H[1] @ V
wts = np.array([[w1[k].imag, H2d[k, k].imag] for k in range(8)])
idx = [k for k in range(8) if np.hypot(*wts[k]) > 1e-9]
ang = np.degrees(np.arctan2(wts[idx, 1], wts[idx, 0])) % 360
order = np.argsort(ang)
print("\n非ゼロ重みの角度（度）:", np.round(np.sort(ang), 2))
print(f"ゼロ重みの本数: {8 - len(idx)}  （実軸 1 と固定軸 ℓ）")
W = wts[[idx[order[i]] for i in [0, 2, 4]]]
t3 = np.degrees(np.arctan2(W[:, 1], W[:, 0])) % 360
v3 = np.array([[np.cos(np.radians(a)), np.sin(np.radians(a))] for a in t3])
print(f"  3 の三つ組: 角度 {np.round(t3,2)}, 差 {np.round(np.diff(np.sort(t3)),2)}, 和 {np.linalg.norm(v3.sum(0)):.1e}")

# --- 中心 Z₃ ---------------------------------------------------------------
# 三つの重みの位相がすべて 2π/3 (mod 2π)。和がゼロなので (2π/3, 2π/3, −4π/3)。
t, _, _, _ = np.linalg.lstsq(W, np.array([2*np.pi/3, 2*np.pi/3, -4*np.pi/3]), rcond=None)
Z = expm(t[0] * H[0] + t[1] * H[1])
print("\n=== 中心 Z = exp(t·H) ===")
print(f"  残差 {np.linalg.norm(W@t - np.array([2*np.pi/3,2*np.pi/3,-4*np.pi/3])):.1e}")
print(f"  直交            : {np.linalg.norm(Z.T@Z-np.eye(8)):.1e}")
print(f"  Z(1) = 1        : {np.linalg.norm(Z@E[0]-E[0]):.1e}")
print(f"  Z(ℓ) = ℓ        : {np.linalg.norm(Z@E[4]-E[4]):.1e}")
print(f"  Z³ = I          : {np.linalg.norm(np.linalg.matrix_power(Z,3)-np.eye(8)):.1e}")
print(f"  Z ≠ I           : ‖Z−I‖ = {np.linalg.norm(Z-np.eye(8)):.3f}")
mx = max(np.linalg.norm(Z@omul(x, y) - omul(Z@x, Z@y)) for x, y in rng.normal(size=(400, 2, 8)))
print(f"  積の自己同型     : {mx:.1e}  → G₂ の元")
print(f"  su(3) と可換     : {max(np.linalg.norm(Z@D-D@Z) for D in su3):.1e}  → stab(ℓ) の中心")
print(f"  固有値の偏角（度）: {np.sort(np.round(np.degrees(np.angle(np.linalg.eigvals(Z))),2))}")
print("\n  Z の行列（観測者が時空 span(1,e1,e2,e3)、内部 span(e4,e5,e6,e7) を取った場合）")
np.set_printoptions(precision=3, suppress=True, linewidth=150)
print(np.round(Z, 3))
print(f"  時空ブロックの変化 |Z_qq − I| = {np.linalg.norm(Z[:4,:4]-np.eye(4)):.3f}, "
      f"時空↔内部の混合 = {np.linalg.norm(Z[:4,4:]):.3f}")
print("  → 三つの平面 (v₁,Q₁), (v₂,Q₂), (v₃,Q₃) を 120° 回す。空間と内部を混ぜる。")

# --- 閉塞解への作用 --------------------------------------------------------
def cres(S):
    return float(np.linalg.norm(S.sum(0))), float(np.linalg.norm(sum(ojordan(v, v) for v in S)))

def gen_k(N, k, r, basis):
    one = np.ones(N) / np.sqrt(N)
    q0 = r.normal(size=N); q0 -= one * (one @ q0)
    P = np.eye(N) - np.outer(one, one) - np.outer(q0, q0) / (q0 @ q0)
    C = np.array([P @ r.normal(size=N) for _ in range(k)]); C *= np.sqrt((q0 @ q0) / np.sum(C * C))
    return np.outer(q0, E[0]) + C.T @ basis

def rand_basis(k, r):
    bs = []
    for _ in range(k):
        v = r.normal(size=8); v[0] = 0
        for b in bs: v -= (b @ v) * b
        bs.append(v / np.linalg.norm(v))
    return np.array(bs)

print("\n=== Z₃ は閉塞解を閉塞解に写すか ===")
for k in [1, 2, 3]:
    r = np.random.default_rng(100 + k)
    S = gen_k(6, k, r, rand_basis(k, r)); SZ = S @ Z.T
    print(f"  k={k}: 元 {cres(S)[0]:.0e},{cres(S)[1]:.0e} → Z 後 {cres(SZ)[0]:.0e},{cres(SZ)[1]:.0e}")

ev, evec = np.linalg.eig(Z)
fix = np.real(evec[:, np.abs(ev - 1) < 1e-9]); Qf, _ = np.linalg.qr(fix)
print(f"\n=== Z の固定部分空間 ===\n  次元 {fix.shape[1]}、基底:")
for c in range(fix.shape[1]): print("   ", np.round(Qf[:, c], 3))
print("  → span(1, ℓ)。実軸と内部の時計軸で、ℂ の一つのコピー（ℓ² = −1）。")

P1 = Qf @ Qf.T
print("\n=== 三つの像の和は固定部分空間に落ちるか（色の中和）===")
for k in [1, 2, 3]:
    r = np.random.default_rng(200 + k)
    S = gen_k(6, k, r, rand_basis(k, r))
    tot = S + S @ Z.T + S @ np.linalg.matrix_power(Z, 2).T
    print(f"  k={k}: ‖S + ZS + Z²S − 3·P₁S‖ = {np.linalg.norm(tot - 3*(S@P1.T)):.1e}  "
          f"(‖S‖={np.linalg.norm(S):.3f}, ‖P₁S‖={np.linalg.norm(S@P1.T):.3f})")
print("  → 三つの像の和 = span(1, ℓ) への射影のちょうど 3 倍。")
print("     一つの像だけでは 3 の成分が残り、三つそろうと消える。")

print("\n=== 三相 (1, ω, ω²) をどこに置くか ===")
w = np.exp(2j * np.pi / 3)
for name, col in [("内部軸 span(1, ℓ) ", 4), ("空間軸 span(1, e₁)", 1)]:
    T = np.zeros((3, 8))
    for n in range(3):
        z = w ** n; T[n, 0] = z.real; T[n, col] = z.imag
    print(f"  {name}: 閉塞 {cres(T)[0]:.0e},{cres(T)[1]:.0e}, ‖ZT − T‖ = {np.linalg.norm(T@Z.T - T):.3f}")
print("  → 内部軸に置けば Z₃ 不変（トライアリティ 0、中性）。")
print("     空間軸に置けば Z で動き、三つの像を持つ（トライアリティ ≠ 0）。")
