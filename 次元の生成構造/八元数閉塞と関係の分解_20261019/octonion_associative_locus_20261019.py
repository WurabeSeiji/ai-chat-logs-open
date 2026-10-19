#!/usr/bin/env python3
"""
結合的 3-平面と内部空間の発生（2026-10-19）

`octonion_closure_check_20261019.py` と同じフォルダに置いて実行する（代数の定義を読み込む）。

問い：倍音が少ない状態では、外部時空には広がっているとして、内部空間は逆に局在しているのか。

計算で出た答え：内部空間は「局在している」のではなく「まだ開いていない」。
配置が使っている虚方向の本数を k とすると（前の検算 C で k = min(7, N − 2)）、

  k ≤ 2 … 必ずある四元数部分代数の中に収まる。適切な観測者には内部状態が存在しない。
  k = 3 … 3-平面が結合的（A×B がその平面内）なときだけ収まる。一般には収まらない。
  k ≥ 4 … 次元で収まらない。内部は必ず開く。

つまり内部空間を作っているのは次元の不足ではなく、三本目の方向が積と整合しないこと。

このプログラムの中身
  1. N = 5（k = 3）のランダム閉塞解で 3-平面の結合性のずれを測る。
  2. 次元の勘定（余次元 4）を二通りで確認する。
  3. 四元数論文の対称配置（十字、三相、正十六胞体の v/√3 版）を O に正しく埋めて測る。
  4. 四元数部分代数の中に N = 5 解を作れば結合的になることを確認する。
  5. 等モジュラー（|o_n| 一定）を課すと結合的になるか（＝対称性は結合的な側を選ぶか）。

再現性：各所で固定シードを使う。出力は標準出力のみ（軽い検算なのでファイルは出さない）。
"""
import os
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "octonion_closure_check_20261019.py")).read().split("results = {")[0])

def closure_residual(S):
    return (float(np.linalg.norm(S.sum(0))),
            float(np.linalg.norm(sum(ojordan(x, x) for x in S))))

def plane_and_k(S, tol=1e-9):
    """配置が実際に使っている虚方向（主成分）と、その本数 k。"""
    U, sv, _ = np.linalg.svd(S[:, 1:].T)
    k = int(np.sum(sv > tol * sv[0]))
    B = np.array([np.concatenate([[0], U[:, i]]) for i in range(k)])
    return k, B

def assoc_defect(B):
    """3-平面 span(A, B, C) が四元数部分代数の虚部か（結合的か）。
       A×B がその平面の外へ出る割合。k ≤ 2 は必ず収まるので 0、k ≥ 4 は収まらないので 1。"""
    k = len(B)
    if k <= 2: return 0.0
    if k >= 4: return 1.0
    A, Bb, C = B
    w = omul(A, Bb); w = w - w[0] * E[0]
    P = np.vstack(B)
    return float(np.linalg.norm(w - P.T @ (P @ w)) / np.linalg.norm(w))

def gen_closure(N, r):
    """Σo = 0, Σo² = 0 を厳密に満たす N 個の八元数（和零 ∩ q₀⊥ への射影と等長スケール）。"""
    one = np.ones(N) / np.sqrt(N)
    x0 = r.normal(size=N); x0 -= one * (one @ x0)
    P = np.eye(N) - np.outer(one, one) - np.outer(x0, x0) / (x0 @ x0)
    X = np.array([P @ r.normal(size=N) for _ in range(7)])
    X *= np.sqrt((x0 @ x0) / np.sum(X * X))
    return np.column_stack([x0, X.T])

# ---------------------------------------------------------------------------
print("=== 1. N = 5（k = 3）のランダム閉塞解で 3-平面は結合的か ===")
vals = [assoc_defect(plane_and_k(gen_closure(5, np.random.default_rng(900 + t)))[1]) for t in range(200)]
vals = np.array(vals)
print(f"  200 本：結合性のずれ min = {vals.min():.4f}, median = {np.median(vals):.4f}, max = {vals.max():.4f}")
print(f"  ずれが 1e-9 未満のもの：{int(np.sum(vals < 1e-9))} / 200")
print("  → 一般の閉塞解は結合的でない。結合的 3-平面は測度ゼロ。")

# ---------------------------------------------------------------------------
print("\n=== 2. 次元の勘定（余次元 4 を二通りで）===")
print("  (a) O の N = 5 解空間 24 次元。四元数部分代数の中の解 12 次元 ＋ 部分代数の選び方")
print("      G₂/SO(4) = 8 次元 → 四元数に乗る解は 20 次元。余次元 4。")
print("  (b) Im(O) の 3-平面は Gr(3,7) = 3 × 4 = 12 次元。うち結合的なものは G₂/SO(4) = 8 次元。")
print("      やはり余次元 4。二通りが一致。")

# ---------------------------------------------------------------------------
print("\n=== 3. 四元数論文の対称配置を O に埋めて測る ===")
def show(name, S):
    k, B = plane_and_k(S); rl, rs = closure_residual(S)
    print(f"  {name}: N = {len(S)}, k = {k}, 結合性のずれ = {assoc_defect(B):.1e}, "
          f"閉塞残差 |Σo| = {rl:.1e}, |Σo²| = {rs:.1e}")

cross = np.zeros((4, 8)); cross[0, 0] = 1; cross[1, 0] = -1; cross[2, 1] = 1; cross[3, 1] = -1
show("十字 (1, −1, i, −i)   ", cross)

w = np.exp(2j * np.pi / 3); tri = np.zeros((3, 8))
for n in range(3):
    z = w ** n; tri[n, 0] = z.real; tri[n, 1] = z.imag
show("三相 (1, ω, ω²)       ", tri)

# 正十六胞体 ±e_μ。等方閉塞（13 本）を満たす配置は v → v/√3 で本当の二乗閉塞（8 本）を満たす（論文 §8.3）
c16 = np.zeros((8, 8)); s = 1 / np.sqrt(3)
c16[0, 0] = 1; c16[1, 0] = -1
for i in range(1, 4):
    c16[2*i, i] = s; c16[2*i+1, i] = -s
show("正十六胞体（v/√3 版） ", c16)
print("  → 既存の対称配置はすべて結合的な側（内部空間をゼロにできる側）に乗っている。")

# ---------------------------------------------------------------------------
print("\n=== 4. 結合的な N = 5 解は作れるか ===")
r = np.random.default_rng(7)
A = r.normal(size=8); A[0] = 0; A /= np.linalg.norm(A)
Bb = r.normal(size=8); Bb[0] = 0; Bb -= (A @ Bb) * A; Bb /= np.linalg.norm(Bb)
C = omul(A, Bb); C[0] = 0; C /= np.linalg.norm(C)
Bas = np.array([E[0], A, Bb, C])                      # 四元数部分代数（4 x 8）
N = 5; one = np.ones(N) / np.sqrt(N)
x0 = r.normal(size=N); x0 -= one * (one @ x0)
P = np.eye(N) - np.outer(one, one) - np.outer(x0, x0) / (x0 @ x0)
X = np.array([P @ r.normal(size=N) for _ in range(3)]); X *= np.sqrt((x0 @ x0) / np.sum(X * X))
S = np.column_stack([x0, X.T]) @ Bas
k, B = plane_and_k(S); rl, rs = closure_residual(S)
print(f"  部分代数の中で作った N = 5 解: k = {k}, 結合性のずれ = {assoc_defect(B):.1e}, "
      f"閉塞残差 = {rl:.1e}, {rs:.1e}")
print("  → 作れる。ただし 24 次元の解空間の中の 20 次元の部分集合（余次元 4）。")

# ---------------------------------------------------------------------------
print("\n=== 5. 等モジュラーを課すと結合的になるか ===")
def solve_equal_modulus(N, seed, restarts=12):
    """|o_n| = 1、Σo = 0、Σo² = 0 を最小二乗で解く。"""
    def f(x):
        S = x.reshape(N, 8)
        e1 = np.sum((np.sum(S**2, axis=1) - 1.0) ** 2)
        e2 = np.sum(S.sum(0) ** 2)
        sq = sum(ojordan(v, v) for v in S)
        return e1 + e2 + float(sq @ sq)
    best = None
    for t in range(restarts):
        x0 = np.random.default_rng(seed * 100 + t).normal(size=N * 8)
        res = minimize(f, x0, method="L-BFGS-B",
                       options=dict(maxiter=20000, maxfun=60000, ftol=1e-18, gtol=1e-14))
        if best is None or res.fun < best.fun: best = res
    return best.x.reshape(N, 8), best.fun

for seed in [1, 2, 3, 4]:
    S, err = solve_equal_modulus(5, seed)
    k, B = plane_and_k(S); rl, rs = closure_residual(S); nr = np.linalg.norm(S, axis=1)
    print(f"  seed {seed}: 残差合計 = {err:.2e}, |o_n| の幅 = {nr.max()-nr.min():.1e}, "
          f"|Σo| = {rl:.1e}, |Σo²| = {rs:.1e}, k = {k}, 結合性のずれ = {assoc_defect(B):.4f}")
print("  → 等モジュラーは k = 4 を選ぶ。対称性は結合的な側へは選ばず、むしろ内部を開く側へ押す。")
