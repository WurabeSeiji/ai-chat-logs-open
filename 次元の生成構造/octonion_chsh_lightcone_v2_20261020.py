#!/usr/bin/env python3
"""
ツィレルソン限界と光円錐 ― 定理・計算・仮定を分ける（2026-10-20 改訂版）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか隣のフォルダに置いて実行する（代数の定義だけを読み込む）。

初版（octonion_chsh_lightcone_20261020.py）は「G2 不変量が内積しかないから 2sqrt(2)」と
書いたが、これは誤り。不変性だけでは 2sqrt(2) は出ない。選んでいるのは双線形性である。
本版はその訂正と、著者の指摘（1〜4）を反映したもの。

なお本検討は Joy Christian の「局所実在論で量子相関を再現した」という主張とは
まったく別である。あれは局所実在的な再現の主張で、強く反駁されている。ここでは
局所実在論を擁護せず、ツィレルソン限界の起源を読み出しの構造から説明することだけを扱う。

確かめること
  1. 不変性だけでは 2sqrt(2) は出ない（角度だけの関数の反例：三角波で 2、符号型で 4）。
  2. 双線形なら CHSH の最大化条件は sin 3theta = sin theta、すなわち theta = pi/4、
     w^4 = -1（w は 1 の原始 8 乗根）。delta = pi/4 を別に仮定しなくてよい。
  3. 角の定義：w の位相はブロッホ角。スピン 1/2 で 45 度、光子の偏光角では 22.5 度。
  4. 連鎖ベル不等式：N 設定での最適間隔は pi/(2N)、最大は 2N cos(pi/2N)。
     読み出しの刻みが N に依らないなら、一致は CHSH の 2x2 構造に固有。
  5. 振幅の重ね合わせから出るのは「振幅についての双線形性」。
     設定についての双線形性は別に要る（この枠組みではまだ出ていない）。

再現性：固定シード。出力は標準出力のみ。
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

# ---------------------------------------------------------------------------
print("=== 1. 不変性だけでは 2sqrt(2) は出ない（反例）===")
print("  角度だけの関数はどれも G2 不変。CHSH の最大値を三つの相関で比べる。")
def S_of(f):
    def neg(p):
        a, ap, b, bp = p
        return -(f(a-b) + f(a-bp) + f(ap-b) - f(ap-bp))
    best = 1e9
    for t in range(200):
        x0 = np.random.default_rng(t).uniform(0, 2*np.pi, 4)
        r = minimize(neg, x0, method="Nelder-Mead",
                     options=dict(maxiter=4000, fatol=1e-12, xatol=1e-12))
        best = min(best, r.fun)
    return -best
def f_cos(t): return math.cos(t)
def f_tri(t):
    t = abs(((t + math.pi) % (2*math.pi)) - math.pi)
    return 1 - 2*t/math.pi
def f_sgn(t):
    t = ((t + math.pi) % (2*math.pi)) - math.pi
    return 1.0 if abs(t) < math.pi/2 else -1.0
for nm, f, note in [("cos theta（双線形）", f_cos, "量子"),
                    ("1 - 2theta/pi（三角波）", f_tri, "古典・局所"),
                    ("sgn(cos theta)（符号型）", f_sgn, "PR 的")]:
    print(f"  {nm:<26}（{note}）: CHSH 最大 = {S_of(f):.6f}")
print(f"  参考：2 sqrt(2) = {2*math.sqrt(2):.6f}、古典上限 = 2、無信号上限 = 4")
print("  -> どれも角度だけの関数＝G2 不変。2 sqrt(2) を選ぶのは不変性ではなく双線形性。")

# ---------------------------------------------------------------------------
print("\n=== 2. 双線形なら CHSH の最大化条件は sin 3theta = sin theta ===")
print("  S(theta) = 3 cos theta - cos 3theta、dS/dtheta = -3 sin theta + 3 sin 3theta = 0")
th = np.linspace(1e-6, np.pi/2, 2000001)
S = 3*np.cos(th) - np.cos(3*th)
i = int(np.argmax(S))
print(f"  数値最大：theta = {math.degrees(th[i]):.6f} 度、S = {S[i]:.10f}")
print(f"  解析解：theta = pi/4 = 45 度、S = 2 sqrt(2) = {2*math.sqrt(2):.10f}")
w = np.exp(1j*np.pi/4)
print(f"  w = exp(i pi/4) について w^4 = {w**4:.6f}  -> w^4 = -1、w は 1 の原始 8 乗根")
print("  -> delta = pi/4 を別に仮定しなくても、CHSH の最大化が位数 8 を選ぶ。")

# ---------------------------------------------------------------------------
print("\n=== 3. 角の定義 ===")
print("  w の位相 theta はブロッホ球上の角。スピン 1/2 の物理角では 45 度、")
print("  光子の偏光角では半分の 22.5 度（偏光角 phi に対しブロッホ角 2 phi）。")

# ---------------------------------------------------------------------------
print("\n=== 4. 連鎖ベル不等式：最適間隔 pi/(2N)、最大 2N cos(pi/2N) ===")
print(f"  {'N':>3}{'最適間隔（度）':>18}{'量子最大':>14}{'2N cos(pi/2N)':>18}{'古典上限':>10}")
for N in (2, 3, 4, 6, 8):
    d = np.linspace(1e-6, np.pi/2, 400001)
    Sv = (2*N - 1)*np.cos(d) - np.cos((2*N - 1)*d)
    j = int(np.argmax(Sv))
    print(f"  {N:>3}{math.degrees(d[j]):>18.4f}{Sv[j]:>14.6f}{2*N*math.cos(math.pi/(2*N)):>18.6f}{2*N-2:>10}")
print("  -> 最適間隔は pi/(2N)。位数は 4N。読み出しの刻み delta が N に依らず pi/4 のままなら、")
print("     光円錐との一致は CHSH の 2x2 構造に固有ということになる（反証可能なテスト）。")

# ---------------------------------------------------------------------------
print("\n=== 5. 重ね合わせから出るのは『振幅についての』双線形性 ===")
rng = np.random.default_rng(1); mx = 0.0
for _ in range(300):
    x, y, z = rng.normal(size=(3, 8)); a, b = rng.normal(size=2)
    lhs = float(omul(oconj(a*x + b*y), z)[0])
    rhs = a*float(omul(oconj(x), z)[0]) + b*float(omul(oconj(y), z)[0])
    mx = max(mx, abs(lhs - rhs))
print(f"  Re(conj(a x + b y) z) = a Re(conj(x) z) + b Re(conj(y) z) の誤差: {mx:.1e}")
print("  -> 強度 |x + y|^2 の交差項 2 Re(conj(x) y) は、振幅について自動的に双線形。")
print("     ただし CHSH が要るのは『設定 a, b について双線形』であって、")
print("     設定が振幅に線形に入ることは別に要る（この枠組みではまだ出ていない）。")
