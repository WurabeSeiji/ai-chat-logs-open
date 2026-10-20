#!/usr/bin/env python3
"""
八元数の干渉と二重スリット（2026-10-20）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか隣のフォルダに置いて実行する（代数の定義だけを読み込む）。

著者の引き戻し：八元数は所詮ただの実数 8 個、複素数 4 個、四元数 2 個。複素数 2 個以上
である以上、波（進行しない定常波）でしかない。すべての作用は干渉結果が読み出しで、
それが位相または振幅を変えるという干渉縞の記録にすぎない。ポイントは常に交差項。
八元数の干渉問題、二重スリット問題はどうなるのか。

計算で出た答え
  1. 重ね合わせ（並列のスリット）は八元数でも複素数と同型。和は常に結合的なので
     |Sum x_i|^2 = Sum |x_i|^2 + 2 Sum_{i<j} Re(conj(x_i) x_j)。本数によらない。
  2. 二波は必ず可干渉（アルティンの定理：どの二つも C か H を張る）。
     二重スリットは八元数でも厳密に複素数と同じ。
  3. 八元数が効くのは並列ではなく直列。振幅が積で合成される段で、三つ目から
     (xy)z != x(yz) となり「どう括るか」で結果が変わる。
  4. 大きさ（干渉の強さ）は経路に依らず、向き（位相）だけが経路に依る。
     S^7 は群でない＝八元数には干渉があるが位相群がない。
  5. 同じ関係 w = conj(x) y の実部が干渉縞の交差項、虚部が作用・生成子が使う量。
     干渉縞と力学は相補な半分を読んでいる。
  6. 可干渉（全三つ組が結合的）<=> k <= 2 または k = 3 結合的 <=> 四元数部分代数に収まる
     <=> 内部空間がゼロ <=> 質量ゼロ・電荷ゼロ。別々に出た三つの条件が同じだった。

再現性：乱数シード 20261020。出力は標準出力のみ。
"""
import os, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CAND = [os.path.join(HERE, "octonion_closure_check_20261019.py"),
        os.path.join(HERE, "..", "八元数閉塞と関係の分解_20261019", "octonion_closure_check_20261019.py"),
        os.path.join(HERE, "八元数閉塞と関係の分解_20261019", "octonion_closure_check_20261019.py")]
SRC = next((p for p in CAND if os.path.exists(p)), None)
if SRC is None:
    raise SystemExit("octonion_closure_check_20261019.py が見つかりません")
exec(open(SRC).read().split("results = {")[0])          # omul, oconj, E

rng = np.random.default_rng(20261020)
def Re(x): return float(x[0])
def nrm(x): return float(np.sqrt(x @ x))

# ---------------------------------------------------------------------------
print("=== 1. 重ね合わせ（並列のスリット）は八元数でも複素数と同じ ===")
mx = 0.0
for n in (2, 3, 5, 8):
    for _ in range(300):
        X = rng.normal(size=(n, 8))
        lhs = nrm(X.sum(0))**2
        rhs = (sum(nrm(x)**2 for x in X)
               + 2*sum(Re(omul(oconj(X[i]), X[j])) for i in range(n) for j in range(i+1, n)))
        mx = max(mx, abs(lhs - rhs))
print(f"  |Sum x_i|^2 = Sum |x_i|^2 + 2 Sum_{{i<j}} Re(conj(x_i) x_j) の誤差（n = 2,3,5,8）: {mx:.1e}")
print("  -> 和は常に結合的。干渉縞は対ごとの実スカラー Re(conj(x) y) の和で、本数によらず複素数と同型。")

# ---------------------------------------------------------------------------
print("\n=== 2. 二波は必ず可干渉（アルティンの定理）===")
mx = 0.0
for _ in range(500):
    x, y = rng.normal(size=(2, 8))
    for a in (x, y, omul(x, y)):
        for b in (x, y, omul(x, y)):
            for c in (x, y):
                mx = max(mx, np.linalg.norm(omul(omul(a, b), c) - omul(a, omul(b, c))))
print(f"  二つの元が生成する部分代数での結合子の最大: {mx:.1e}")
print("  -> どの二つも C か H を張る。二重スリットは八元数でも厳密に複素数と同じ。")

# ---------------------------------------------------------------------------
print("\n=== 3. 三波：重ね合わせは平気だが、順次の作用は括弧で変わる ===")
tri = []
for _ in range(500):
    x, y, z = rng.normal(size=(3, 8))
    x, y, z = x/nrm(x), y/nrm(y), z/nrm(z)
    tri.append(np.linalg.norm(omul(omul(x, y), z) - omul(x, omul(y, z))))
tri = np.array(tri)
print(f"  単位三元の結合子 |(xy)z - x(yz)|: 中央値 {np.median(tri):.3f}、最大 {tri.max():.3f}")
print("  -> 三つ目から「どう括るか」で結果が変わる。振幅が積で合成される段で初めて効く。")

# ---------------------------------------------------------------------------
print("\n=== 4. 位相の合成が経路に依る（S^7 は群でない）===")
us = []
for _ in range(400):
    x, y, z = rng.normal(size=(3, 8))
    x, y, z = x/nrm(x), y/nrm(y), z/nrm(z)
    p = omul(omul(x, y), z); q = omul(x, omul(y, z))
    us.append((abs(nrm(p) - 1), abs(nrm(q) - 1), nrm(p - q)))
us = np.array(us)
print(f"  |(xy)z| = 1 も |x(yz)| = 1 も保たれる（誤差 {us[:, :2].max():.1e}）が、")
print(f"  二つの向きは違う（差の中央値 {np.median(us[:, 2]):.3f}）。")
print("  -> 大きさ（干渉の強さ）は経路に依らず、向き（位相）だけが経路に依る。")

# ---------------------------------------------------------------------------
print("\n=== 5. 干渉縞が読むのは Re(w)、作用が読むのは Im(w) ===")
m1 = m2 = 0.0
for _ in range(300):
    x, y = rng.normal(size=(2, 8))
    w = omul(oconj(x), y)
    th = np.arccos(np.clip((x @ y)/(nrm(x)*nrm(y)), -1, 1))
    m1 = max(m1, abs(Re(w) - nrm(x)*nrm(y)*np.cos(th)))
    im = w.copy(); im[0] = 0.0
    m2 = max(m2, abs(np.linalg.norm(im) - nrm(x)*nrm(y)*abs(np.sin(th))))
print(f"  Re(w) = |x||y| cos theta の誤差: {m1:.1e}（干渉縞の交差項そのもの）")
print(f"  |Im(w)| = |x||y| |sin theta| の誤差: {m2:.1e}（生成子・作用が使う量）")
print("  -> 同じ関係 w の実部と虚部。干渉縞と力学は相補な半分を読んでいる。")

# ---------------------------------------------------------------------------
print("\n=== 6. 可干渉の条件 = 結合性 = k の条件 ===")
def k_of(S):
    s = np.linalg.svd(S[:, 1:].T, compute_uv=False)
    return 0 if s[0] < 1e-300 else int(np.sum(s > 1e-9*s[0]))

def assoc_max(S):
    m = 0.0
    for i in range(len(S)):
        for j in range(len(S)):
            for l in range(len(S)):
                m = max(m, np.linalg.norm(omul(omul(S[i], S[j]), S[l]) - omul(S[i], omul(S[j], S[l]))))
    return m

def quat_basis(r, k):
    A = r.normal(size=8); A[0] = 0; A /= np.linalg.norm(A)
    B = r.normal(size=8); B[0] = 0; B -= (A @ B)*A; B /= np.linalg.norm(B)
    C = omul(A, B); C[0] = 0; C /= np.linalg.norm(C)
    return np.array([A, B, C])[:k]

def gen_k(N, k, r, basis=None):
    one = np.ones(N)/np.sqrt(N); q0 = r.normal(size=N); q0 -= one*(one @ q0)
    P = np.eye(N) - np.outer(one, one) - np.outer(q0, q0)/(q0 @ q0)
    C = np.array([P @ r.normal(size=N) for _ in range(k)]); C *= np.sqrt((q0 @ q0)/np.sum(C*C))
    if basis is None:
        bs = []
        for i in range(k):
            v = r.normal(size=8); v[0] = 0
            for b in bs: v -= (b @ v)*b
            bs.append(v/np.linalg.norm(v))
        basis = np.array(bs)
    return np.outer(q0, E[0]) + C.T @ basis

print(f"  {'配置':<28}{'k':>3}{'全三つ組の結合子の最大':>24}")
for k in (1, 2, 3):
    S = gen_k(6, k, np.random.default_rng(50 + k))
    print(f"  一般の k = {k} の閉塞解            {k_of(S):>3}{assoc_max(S):>24.1e}")
r = np.random.default_rng(99)
S = gen_k(6, 3, r, basis=quat_basis(r, 3))
print(f"  k = 3 結合的（H の中）            {k_of(S):>3}{assoc_max(S):>24.1e}")
for k in (4, 5):
    S = gen_k(8, k, np.random.default_rng(70 + k))
    print(f"  一般の k = {k} の閉塞解            {k_of(S):>3}{assoc_max(S):>24.1e}")
print("  -> 可干渉（全三つ組が結合的）<=> k <= 2 または k = 3 結合的")
print("     <=> 四元数部分代数に収まる <=> 内部空間ゼロ <=> 質量ゼロ・電荷ゼロ")
