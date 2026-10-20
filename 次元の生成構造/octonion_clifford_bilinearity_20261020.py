#!/usr/bin/env python3
"""
双線形性が八元数の積から出ること（2026-10-20）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか隣のフォルダに置いて実行する（代数の定義だけを読み込む）。

前段：ツィレルソン限界と光円錐_v2_20261020.md の §4-1（設定が振幅に線形に入る理由、
最優先の未導出項目）と §3.4（連鎖ベルでは判定に向かない）。

確かめること
  1. クリフォード関係 L_a L_b + L_b L_a = -2<a,b>（a, b 純虚単位）。L_a は反対称。
     |L_a x| = |a||x|（ノルム乗法性）。したがって a -> L_a x は線形等長。
  2. 反対称部分 Re(conj(x) [a(bx) - b(ax)]) は、特別な状態を要さず恒等的に消える
     （[L_a, L_b] が反対称行列で <x, Ax> = 0 だから）。ただし状態依存性も消えるので
     この一体の形ではベル相関にならない。
  3. 二者の形 E(a,b) = <a x, b y> / (|x||y|) はツィレルソンの形そのもの。
     設定について双線形、|E| <= 1、状態 (x,y) に依存、CHSH は 2 sqrt(2) に届く。
  4. 相関行列は 7x7 で x = y = 1 のとき単位行列。読み出しは平面に閉じていない。
     非平面のベル不等式（I3322、プラトン型）も扱える見込み。

再現性：乱数シード 20261020、最適化は固定シードの多出発点。出力は標準出力のみ。
"""
import os, math
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
CAND = [os.path.join(HERE, "octonion_closure_check_20261019.py"),
        os.path.join(HERE, "..", "八元数閉塞と関係の分解_20261019", "octonion_closure_check_20261019.py"),
        os.path.join(HERE, "八元数閉塞と関係の分解_20261019", "octonion_closure_check_20261019.py")]
SRC = next((p for p in CAND if os.path.exists(p)), None)
if SRC is None:
    raise SystemExit("octonion_closure_check_20261019.py が見つかりません")
exec(open(SRC).read().split("results = {")[0])          # omul, oconj, E

rng = np.random.default_rng(20261020)
def Lmat(a):
    """左乗算 L_a の 8x8 行列"""
    return np.column_stack([omul(a, E[k]) for k in range(8)])
def imunit(r):
    v = r.normal(size=8); v[0] = 0.0; return v/np.linalg.norm(v)

# ---------------------------------------------------------------------------
print("=== 1. クリフォード関係（純虚の左乗算は Cl(7) を作る）===")
mx = max(np.max(np.abs(Lmat(a) @ Lmat(b) + Lmat(b) @ Lmat(a) + 2*(a @ b)*np.eye(8)))
         for a, b in [(imunit(rng), imunit(rng)) for _ in range(500)])
print(f"  L_a L_b + L_b L_a = -2<a,b> の誤差: {mx:.1e}")
mx = max(np.max(np.abs(Lmat(a).T + Lmat(a))) for a in [imunit(rng) for _ in range(300)])
print(f"  L_a の反対称性 |L_a^T + L_a|: {mx:.1e}")
mx = 0.0
for _ in range(300):
    a = imunit(rng); x = rng.normal(size=8)
    mx = max(mx, abs(np.linalg.norm(omul(a, x)) - np.linalg.norm(x)))
print(f"  |L_a x| = |a||x| の誤差: {mx:.1e}  -> a |-> L_a x は線形等長")

# ---------------------------------------------------------------------------
print("\n=== 2. 反対称部分は恒等的に消える（特別な状態は要らない）===")
m1 = m2 = 0.0
for _ in range(500):
    a = imunit(rng); b = imunit(rng); x = rng.normal(size=8)
    m1 = max(m1, abs(float(x @ omul(a, omul(b, x))) - float(x @ omul(b, omul(a, x)))))
    m2 = max(m2, abs(float(x @ omul(a, omul(b, x))) + (a @ b)*(x @ x)))
print(f"  Re(conj(x)[a(bx) - b(ax)]) の最大: {m1:.1e}")
print(f"  Re(conj(x) a(bx)) + <a,b>|x|^2 の最大: {m2:.1e}  -> 全ての x で成立")
print("  理由：[L_a, L_b] は反対称行列で、任意の反対称 A と任意の x について <x, Ax> = 0。")
print("  ただし状態依存性も無いので、この一体の形ではベル相関にならない。")

# ---------------------------------------------------------------------------
print("\n=== 3. 二者の形はツィレルソンの形そのもの ===")
mx = 0.0
for _ in range(300):
    a1, a2 = imunit(rng), imunit(rng); c1, c2 = rng.normal(size=2)
    b = imunit(rng); x, y = rng.normal(size=(2, 8))
    lhs = float(omul(c1*a1 + c2*a2, x) @ omul(b, y))
    rhs = c1*float(omul(a1, x) @ omul(b, y)) + c2*float(omul(a2, x) @ omul(b, y))
    mx = max(mx, abs(lhs - rhs))
print(f"  E(a,b) = <a x, b y> の設定 a についての線形性: 誤差 {mx:.1e}")
print("  u_a = a x、v_b = b y は単位ベクトル（§1 のノルム乗法性）なので E = <u_a, v_b>。")

def chsh_oct(p):
    v = p.reshape(4, 8).copy(); v[:, 0] = 0.0
    v = v/np.linalg.norm(v, axis=1, keepdims=True)
    a, ap, b, bp = v; x = E[0]; y = E[0]
    Ec = lambda u, w: float(omul(u, x) @ omul(w, y))
    return -(Ec(a, b) + Ec(a, bp) + Ec(ap, b) - Ec(ap, bp))
best = None
for t in range(60):
    r = minimize(chsh_oct, np.random.default_rng(t).normal(size=32), method="BFGS")
    if best is None or r.fun < best.fun: best = r
print(f"  CHSH 最大 = {-best.fun:.10f}   2 sqrt(2) = {2*math.sqrt(2):.10f}   "
      f"差 {abs(-best.fun - 2*math.sqrt(2)):.1e}")

# ---------------------------------------------------------------------------
print("\n=== 4. 相関行列は 7 次元。平面に閉じていない ===")
for nm, (x, y) in [("x = y = 1", (E[0], E[0])),
                   ("ランダムな x, y", (rng.normal(size=8), rng.normal(size=8)))]:
    x = x/np.linalg.norm(x); y = y/np.linalg.norm(y)
    M = np.array([[float(omul(E[i], x) @ omul(E[j], y)) for j in range(1, 8)] for i in range(1, 8)])
    print(f"  {nm:<16}: 7x7 相関行列の特異値 = {np.round(np.linalg.svd(M, compute_uv=False), 6)}")
v = best.x.reshape(4, 8).copy(); v[:, 0] = 0.0
v = v/np.linalg.norm(v, axis=1, keepdims=True)
print(f"  CHSH 最適解の設定 4 本が張る次元 = {np.linalg.matrix_rank(v[:, 1:], tol=1e-8)}（平面なら 2）")
M3 = np.array([[float(omul(E[i], E[0]) @ omul(E[j], E[0])) for j in (1, 2, 3)] for i in (1, 2, 3)])
print(f"  三次元の直交設定での相関行列の特異値 = {np.round(np.linalg.svd(M3, compute_uv=False), 6)}")
print("  -> 設定空間は 7 次元で、キュービット相関（3 次元）を完全に含む。")
print("     I3322 やプラトン型の非平面不等式も量子値に届く見込み。")
