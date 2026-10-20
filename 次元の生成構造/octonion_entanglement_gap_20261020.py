#!/usr/bin/env python3
"""
残った穴の正体 ― テンソル積が無いこと（2026-10-20）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか隣のフォルダに置いて実行する（代数の定義だけを読み込む）。

前段：双線形性が八元数の積から出た_20261020.md の §4「読み」と、
そこで残した「二者性はどこから来るか」。

この検算で確定したこと
  1. 相関行列は 6 個の特異値が厳密に 1 に張り付き、7 番目だけが |<x,y>|。
     x, y をどう選んでも変わらない＝最大もつれを構成で入れている。
  2. 局所期待値 <x, a x>、<x, x a> は恒等的にゼロ（L_a も R_a も反対称）。
  3. 左右を入れ替えた四つの形はすべて同じスペクトル。
  4. 結合子を使う形は常に階数 4 で、非ゼロの 4 個は完全に縮退。
     混ぜると |E| <= 1 を破る。結合子だけでは部分もつれを作れない。
  5. (x, y) in S^7 x S^7（14 パラメータ）から到達する相関行列の族は 22 次元。
     7x7 の相関行列の全体 49 次元に対して薄い。階数 1（直積状態）には届かない。

結論：残った穴は「二者性」ではなく「テンソル積が無いこと」。二つの八元数の対は
テンソル積ではないので、シュミット係数を入れる場所がない。

再現性：乱数シード 20261020。出力は標準出力のみ。
"""
import os
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
def Im(x):
    y = x.copy(); y[0] = 0.0; return y
def unit(x): return x/np.linalg.norm(x)
def assoc(a, b, c): return omul(omul(a, b), c) - omul(a, omul(b, c))
def Mmat(x, y):
    return np.array([[float(omul(E[i], x) @ omul(E[j], y)) for j in range(1, 8)] for i in range(1, 8)])

# ---------------------------------------------------------------------------
print("=== 1. 相関行列の特異値は 6 個が 1、7 番目が |<x,y>| ===")
for t in range(5):
    x = unit(rng.normal(size=8)); y = unit(rng.normal(size=8))
    s = np.linalg.svd(Mmat(x, y), compute_uv=False)
    print(f"  試行 {t}: 特異値 = {np.round(s, 6)}   |<x,y>| = {abs(float(x @ y)):.6f}")
print("  理由：{L_e_i x} と {L_e_j y} は R^8 内の正規直交 7 本で、直交補空間はそれぞれ x と y。")
print("        7 次元部分空間 x_perp と y_perp の交わりは 6 次元なので主角が 6 個ゼロ。")

# ---------------------------------------------------------------------------
print("\n=== 2. 局所期待値は恒等的にゼロ ===")
m1 = m2 = 0.0
for _ in range(400):
    a = unit(Im(rng.normal(size=8))); x = unit(rng.normal(size=8))
    m1 = max(m1, abs(float(x @ omul(a, x))))
    m2 = max(m2, abs(float(x @ omul(x, a))))
print(f"  <x, a x> の最大 = {m1:.1e}、<x, x a> の最大 = {m2:.1e}")
print("  -> L_a も R_a も反対称。一回の乗算で作る観測量は局所期待値を持てない。")

# ---------------------------------------------------------------------------
print("\n=== 3. 左右を入れ替えても同じ ===")
forms = {
    "<a x, b y>": lambda a, x, b, y: omul(a, x) @ omul(b, y),
    "<a x, y b>": lambda a, x, b, y: omul(a, x) @ omul(y, b),
    "<x a, y b>": lambda a, x, b, y: omul(x, a) @ omul(y, b),
    "<x a, b y>": lambda a, x, b, y: omul(x, a) @ omul(b, y),
    "<[a,x,b], y>（結合子）": lambda a, x, b, y: assoc(a, x, b) @ y,
    "<[a,x,y], b>（結合子）": lambda a, x, b, y: assoc(a, x, y) @ b,
    "<[a,b,x], y>（結合子）": lambda a, x, b, y: assoc(a, b, x) @ y,
}
x = unit(rng.normal(size=8)); y = unit(rng.normal(size=8))
print(f"  x, y はランダム（<x,y> = {float(x @ y):.4f}）")
for nm, f in forms.items():
    M = np.array([[float(f(E[i], x, E[j], y)) for j in range(1, 8)] for i in range(1, 8)])
    print(f"  {nm:<24} 特異値 = {np.round(np.linalg.svd(M, compute_uv=False), 4)}")

# ---------------------------------------------------------------------------
print("\n=== 4. 結合子を混ぜると |E| <= 1 を破る ===")
M0 = np.array([[float(forms["<a x, b y>"](E[i], x, E[j], y)) for j in range(1, 8)] for i in range(1, 8)])
M1 = np.array([[float(forms["<[a,x,b], y>（結合子）"](E[i], x, E[j], y)) for j in range(1, 8)] for i in range(1, 8)])
print(f"  {'lambda':>8}{'最大特異値':>14}{'最小特異値':>14}{'1 に張り付く数':>16}{'|E| <= 1':>10}")
for lam in (0, 0.1, 0.25, 0.5, 0.75, 1.0):
    s = np.linalg.svd((1-lam)*M0 + lam*M1, compute_uv=False)
    print(f"  {lam:>8.2f}{s.max():>14.4f}{s.min():>14.4f}{int(np.sum(np.abs(s-1) < 1e-9)):>16}"
          f"{('はい' if s.max() <= 1 + 1e-9 else 'いいえ'):>10}")

# ---------------------------------------------------------------------------
print("\n=== 5. 到達できる相関行列の族の次元と、直積状態への距離 ===")
Ms = np.array([Mmat(unit(rng.normal(size=8)), unit(rng.normal(size=8))).reshape(-1) for _ in range(400)])
Ms -= Ms.mean(0)
s = np.linalg.svd(Ms, compute_uv=False)
print(f"  (x, y) in S^7 x S^7（14 パラメータ）から到達する M の張る次元 = "
      f"{int(np.sum(s > 1e-8*s[0]))}  （全体は 49 次元）")
tgt = np.outer(np.eye(7)[0], np.eye(7)[1])
def f(p): return np.linalg.norm(Mmat(unit(p[:8]), unit(p[8:])) - tgt)
best = min(minimize(f, np.random.default_rng(t).normal(size=16), method="BFGS").fun for t in range(30))
print(f"  階数 1（直積状態）への最小距離 = {best:.4f}（目標のノルム {np.linalg.norm(tgt):.4f}）")
print("  -> 直積状態には届かない。6 個の特異値が 1 に張り付いたままだから。")

# ---------------------------------------------------------------------------
print("\n=== 6. ここまでの計算は H でも成り立っていた（結合子が初めて八元数固有）===")
mx = 0.0
for _ in range(300):
    a = Im(rng.normal(size=8)); a[4:] = 0
    b = Im(rng.normal(size=8)); b[4:] = 0
    c = rng.normal(size=8); c[4:] = 0
    mx = max(mx, np.linalg.norm(assoc(a, c, b)))
print(f"  H の中での [a, x, b] の最大 = {mx:.1e}")
print("  -> 交代性もクリフォード関係も結合的な代数で成り立つ。")
print("     結合子が初めて八元数固有の寄与を入れるが、それでも部分もつれは作れなかった。")
