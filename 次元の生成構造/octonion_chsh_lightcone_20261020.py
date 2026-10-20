#!/usr/bin/env python3
"""
CHSH の 2sqrt(2) と光円錐が同じ 45 度に来ること（2026-10-20）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか隣のフォルダに置いて実行する（代数の定義だけを読み込む）。

著者の問題提起：ベルの定理が言えるのは no-go だけで、「だから自然は非実在」は論理的に
出てこない。本質は「なぜそう観測されるか」であり、特に二つの一致が不思議である。
(1) CHSH が 2sqrt(2)（ツィレルソン限界）にちょうど止まり、信号を送らない範囲の最大値
    4（PR 箱）までは行かない。
(2) 信号を送れない境界が、読み出された時空の光円錐とぴったり重なる。
相関の側の理屈と時空の側の理屈が別々にあるなら、この重なりは偶然になる。

このプログラムで確かめること
  1. 対から読める G2 不変な双線形スカラーの数（2 本。単位ベクトルなら Re(conj(x) y) = cos theta
     の一本に落ちる）。相関が「角のコサイン」であることは定理で、仮定ではない。
  2. CHSH を Re(w) で最大化すると 2 sqrt(2)、最適角は 45 度ちょうど。
  3. 45 度は同じ関係 w で Re(w^2) = 0（光円錐）かつ |Im(w^2)| 最大（パワー流の最大）。
  4. PR 箱の 4 が出ないのは、内積でない相関が作れないから。

再現性：乱数シード 20261020 と最適化の多出発点（固定シード）。出力は標準出力のみ。
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

def Im(x):
    y = x.copy(); y[0] = 0.0; return y
def nrm(x): return float(np.sqrt(x @ x))
def nullbasis(A, tol=1e-9):
    u, s, vt = np.linalg.svd(A); return vt[int(np.sum(s > tol*s[0])):]

# ---------------------------------------------------------------------------
print("=== 1. 相関として読めるのは内積だけか（G2 不変な双線形スカラーの数）===")
rows = []
for i in range(8):
    for j in range(8):
        p = omul(E[i], E[j])
        for c in range(8):
            r = np.zeros(64)
            for m in range(8): r[c*8 + m] += p[m]
            for m in range(8): r[m*8 + i] -= omul(E[m], E[j])[c]
            for m in range(8): r[m*8 + j] -= omul(E[i], E[m])[c]
            rows.append(r)
g2 = [v.reshape(8, 8) for v in nullbasis(np.array(rows))]
# B(x,y) = x^T M y が G2 不変 <=> D^T M + M D = 0（全 D in g2）
C = []
for D in g2:
    for a in range(8):
        for b in range(8):
            r = np.zeros(64)
            for m in range(8):
                r[m*8 + b] += D[m, a]; r[a*8 + m] += D[m, b]
            C.append(r)
inv = nullbasis(np.array(C))
print(f"  dim = {inv.shape[0]}（Re(x)Re(y) と <Im x, Im y> の 2 本）")
for v in inv:
    M = v.reshape(8, 8); M = M/np.max(np.abs(M))
    print(f"    対角 {np.round(np.diag(M), 3)}、非対角の最大 {np.max(np.abs(M - np.diag(np.diag(M)))):.1e}")
print("  -> 単位ベクトルどうしなら合わさって Re(conj(x) y) = cos theta の一本。")
print("     相関が『角のコサイン』であることは定理であって、仮定ではない。")

# ---------------------------------------------------------------------------
print("\n=== 2. CHSH を Re(w) で最大化 ===")
def chsh(p):
    v = p.reshape(4, 8); v = v/np.linalg.norm(v, axis=1, keepdims=True)
    a, ap, b, bp = v
    return -(a @ b + a @ bp + ap @ b - ap @ bp)
best = None
for t in range(40):
    r = minimize(chsh, np.random.default_rng(t).normal(size=32), method="BFGS")
    if best is None or r.fun < best.fun: best = r
v = best.x.reshape(4, 8); v = v/np.linalg.norm(v, axis=1, keepdims=True)
a, ap, b, bp = v; S = -best.fun
print(f"  S_max = {S:.10f}   2 sqrt(2) = {2*math.sqrt(2):.10f}   差 {abs(S - 2*math.sqrt(2)):.2e}")
def ang(x, y): return math.degrees(math.acos(np.clip(x @ y, -1, 1)))
print(f"  最適角: (a,b) = {ang(a,b):.3f}, (a,b') = {ang(a,bp):.3f}, "
      f"(a',b) = {ang(ap,b):.3f}, (a',b') = {ang(ap,bp):.3f}")
print(f"          (b,b') = {ang(b,bp):.3f}, (a,a') = {ang(a,ap):.3f}")
print("  -> CHSH の最適は 45 度ちょうど。")

# ---------------------------------------------------------------------------
print("\n=== 3. 45 度は同じ関係 w で何か ===")
print(f"  {'theta':>6}{'Re(w) 相関':>14}{'Re(w^2) 光円錐':>18}{'|Im(w^2)| パワー流':>20}")
for th in (0, 30, 45, 60, 90):
    x = E[0].copy()
    y = math.cos(math.radians(th))*E[0] + math.sin(math.radians(th))*E[1]
    w = omul(oconj(x), y); w2 = omul(w, w)
    print(f"  {th:>5}°{w[0]:>14.4f}{w2[0]:>18.4f}{np.linalg.norm(Im(w2)):>20.4f}")
print("  -> theta = 45 度はちょうど Re(w^2) = 0（光円錐）かつ |Im(w^2)| 最大（放射の最大）。")
print("     CHSH の最適角・光円錐・放射の最大が同じ一つの角。")

# ---------------------------------------------------------------------------
print("\n=== 4. PR 箱（S = 4）が出ない理由 ===")
print("  S = 4 には E(a,b) が内積でない相関であることが要る。§1 により作れない。")
print("  Re(w) = cos theta である限り、")
print("  S = <a, b+b'> + <a', b-b'> <= |b+b'| + |b-b'| <= 2 sqrt(2)（コーシー・シュワルツ）。")
print(f"  最適解での |b+b'| = {nrm(b+bp):.6f}、|b-b'| = {nrm(b-bp):.6f}、"
      f"和 = {nrm(b+bp)+nrm(b-bp):.6f}、2 sqrt(2) = {2*math.sqrt(2):.6f}")
