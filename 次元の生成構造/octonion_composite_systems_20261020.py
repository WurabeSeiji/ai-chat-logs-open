#!/usr/bin/env python3
"""
複合系とは何か ― テンソル積を許すのは複素数だけ（2026-10-20）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか隣のフォルダに置いて実行する（代数の定義だけを読み込む）。

著者の問い（2026-10-20）：J_3(O) がテンソル積を許さないというのは面白いが、
逆に、想定している「複合系」とは何なのか。では四元数や複素数ではできるのか。

計算で出た答え
  1. 「複合系」の操作的な中身は局所トモグラフィ。これを課すと
     dim V_AB = dim V_A x dim V_B が必然になる。
  2. n x n エルミートの実次元 n + d n(n-1)/2 で、(2+d)^2 = 4+6d を解くと d = 2 だけ。
     R は 9 対 10 で足りず、H は 36 対 28 で多すぎ、O は J_4(O) が存在しない。
     テンソル積を許すのは複素数だけ。
  3. したがって「八元数にテンソル積が無い」のは欠陥ではなく、C 以外に共通の事実。
  4. この枠組みの複合は加法的（要素を足す。次元は 8N）。閉塞は大域の拘束なので
     局所トモグラフィを構成で満たさない。再構成定理の仮定の外にいる。
  5. 閉塞は和でそのまま合成する（N=5 と N=6 を合わせて N=11 でも残差 1e-15）。

再現性：固定シード。出力は標準出力のみ。
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CAND = [os.path.join(HERE, "octonion_closure_check_20261019.py"),
        os.path.join(HERE, "..", "八元数閉塞と関係の分解_20261019", "octonion_closure_check_20261019.py"),
        os.path.join(HERE, "八元数閉塞と関係の分解_20261019", "octonion_closure_check_20261019.py")]
SRC = next((p for p in CAND if os.path.exists(p)), None)
if SRC is None:
    raise SystemExit("octonion_closure_check_20261019.py が見つかりません")
exec(open(SRC).read().split("results = {")[0])          # omul, oconj, ojordan, E

# ---------------------------------------------------------------------------
print("=== 1. 「複合系」の操作的な中身：局所トモグラフィ ===")
print("  標準の要請：AB の状態は、A と B への局所測定の統計だけで決まる。")
print("  これを課すと必然的に dim V_AB = dim V_A x dim V_B になる")
print("  （局所測定の組 {A の効果} x {B の効果} が AB の状態空間を張るから）。")

def dimJ(n, d):
    """n x n エルミート行列の実次元。d = 1(R), 2(C), 4(H), 8(O)"""
    return n + d*n*(n-1)//2

print("\n=== 2. 各代数で dim(AB) = dim(A) x dim(B) が成り立つか ===")
print(f"  {'代数':<6}{'1 粒子 n=2':>12}{'2 粒子 n=4':>12}{'積 dim(A)^2':>14}{'一致':>8}")
for nm, d in [("R", 1), ("C", 2), ("H", 4), ("O", 8)]:
    a = dimJ(2, d); ab = dimJ(4, d)
    note = "（J_4(O) はジョルダン代数にならない）" if nm == "O" else ""
    print(f"  {nm:<6}{a:>12}{ab:>12}{a*a:>14}{('はい' if ab == a*a else 'いいえ'):>8}  {note}")
print("\n  R : 3 x 3 = 9 だが実 QM の二粒子は 10。局所トモグラフィが破れる（荒木、Wootters）")
print("  C : 4 x 4 = 16 で一致。【ここだけ成り立つ】")
print("  H : 6 x 6 = 36 だが四元数エルミート 4x4 は 28。テンソル積が圏として定義できない")
print("  O : アルバートの定理により J_3(O) までしかジョルダン代数にならない。複合系を作れない")

print("\n=== 3. なぜ C だけか（一行の方程式）===")
print("  n x n エルミートの実次元は n + d n(n-1)/2。dim(n=2) = 2+d、dim(n=4) = 4+6d。")
print("  (2+d)^2 = 4+6d  ->  d^2 - 2d = 0  ->  d = 0 または 2")
for d in (1, 2, 4, 8):
    print(f"    d = {d}: (2+d)^2 = {(2+d)**2:>3}、4+6d = {4+6*d:>3}"
          f"{'   一致' if (2+d)**2 == 4+6*d else ''}")
print("  -> d = 2、すなわち複素数だけ。『量子力学はなぜ複素数か』の再構成の一本。")

print("\n=== 4. この枠組みの『複合』は加法的 ===")
print("  状態は N 個の八元数の配置。二つを合わせるとは要素を足すこと。次元は 8N。")
for N1, N2 in [(3, 3), (4, 5), (5, 8)]:
    print(f"    N = {N1} と N = {N2}: 8*{N1} + 8*{N2} = {8*(N1+N2):>4}"
          f"（テンソル積なら {8*N1} x {8*N2} = {8*N1*8*N2}）")
print("  閉塞 Sum o = 0, Sum o^2 = 0 は大域の拘束で、局所の測定からは決まらない。")
print("  -> この枠組みは局所トモグラフィを構成で満たさない。再構成定理の仮定の外。")

print("\n=== 5. 閉塞は和でそのまま合成する ===")
def gen(N, r, k=7):
    one = np.ones(N)/np.sqrt(N); q0 = r.normal(size=N); q0 -= one*(one @ q0)
    P = np.eye(N) - np.outer(one, one) - np.outer(q0, q0)/(q0 @ q0)
    C = np.array([P @ r.normal(size=N) for _ in range(k)]); C *= np.sqrt((q0 @ q0)/np.sum(C*C))
    bs = []
    for i in range(k):
        v = r.normal(size=8); v[0] = 0
        for b in bs: v -= (b @ v)*b
        bs.append(v/np.linalg.norm(v))
    return np.outer(q0, E[0]) + C.T @ np.array(bs)
def cres(S):
    return np.linalg.norm(S.sum(0)), np.linalg.norm(sum(ojordan(v, v) for v in S))
A = gen(5, np.random.default_rng(1)); B = gen(6, np.random.default_rng(2))
U = np.vstack([A, B])
print(f"  A（N=5）の閉塞      {cres(A)[0]:.0e}, {cres(A)[1]:.0e}")
print(f"  B（N=6）の閉塞      {cres(B)[0]:.0e}, {cres(B)[1]:.0e}")
print(f"  A ∪ B（N=11）の閉塞 {cres(U)[0]:.0e}, {cres(U)[1]:.0e}  -> 和でそのまま閉じる")
print("  テンソル積なしに合成できている。ただし『もつれ』に相当する自由度は入らない。")
