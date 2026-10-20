#!/usr/bin/env python3
"""
幾何学的位相と八元数の補正（2026-10-20）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか隣のフォルダに置いて実行する（代数の定義だけを読み込む）。

前段：八元数の干渉と二重スリット_20261020.md §7 の検算 3
「ホロノミー＝干渉の順序の記録、の確認」。

この検算で潰れた読み
  「ホロノミーとは干渉の順序の記録である」は外れ。順序の差（結合子）は面積の 2.5 乗、
  ホロノミーは 1 乗で、別の量。ホロノミーと同じ量なのは読み出しの一周（比 1/4）。

代わりに出たこと
  1. 「内側の関係」w_ab = conj(o_a) o_b で作った一周は、結合的な代数では telescoping して
     厳密に実数。位相は出ない。
  2. 「外側の比」s_ab = o_b conj(o_a)（o_a を o_b へ運ぶ輸送子）で作った一周は
     telescoping せず、四元数の中でも枠が 4 x 面積で回る。これは平行移動の
     幾何学的位相（ベリー／パンチャラトナム）で、八元数の性質ではない。
  3. 八元数で新しいのは、輸送が括り方に依ることだけ。しかもその次数は
     素朴な数え上げ（eps^3）より 2 次高い eps^5 = 面積の 5/2 乗。
     三つの一歩が一次で和ゼロになり、結合子が交代形式なので eps^3 が消える。

再現性：乱数シード 20261020 と 7。出力は標準出力のみ。
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
exec(open(SRC).read().split("results = {")[0])          # omul, oconj, oexp, E

def Im(x):
    y = x.copy(); y[0] = 0.0; return y
def nrm(x): return float(np.sqrt(x @ x))
def assoc(a, b, c): return omul(omul(a, b), c) - omul(a, omul(b, c))

def frame_holonomy(oe, of, og):
    """比の一歩 s_ef = o_f conj(o_e) を三角形で一周したときの枠の写像（8x8）"""
    sef = omul(of, oconj(oe)); sfg = omul(og, oconj(of)); sge = omul(oe, oconj(og))
    H = np.zeros((8, 8))
    for k in range(8):
        H[:, k] = omul(sge, omul(sfg, omul(sef, E[k])))
    return H

def sph_area(e, f, g):
    def a(p, q, r):
        tq = q - (p @ q)*p; tr = r - (p @ r)*p
        return np.arccos(np.clip(tq @ tr/np.linalg.norm(tq)/np.linalg.norm(tr), -1, 1))
    return float(a(e, f, g) + a(f, g, e) + a(g, e, f) - np.pi)

def rot_angle(M):
    return float(np.max(np.abs(np.angle(np.linalg.eigvals(M)))))

def readout(a, b):
    w = omul(oconj(a), b); return w[1:]/w[0]

def slope(x, y):
    x = np.asarray(x); y = np.asarray(y); m = (x > 0) & (y > 0)
    return float(np.polyfit(np.log(x[m]), np.log(y[m]), 1)[0])

rng = np.random.default_rng(20261020)

# ---------------------------------------------------------------------------
print("=== 1. 三つの量を同じ三角形で比べる ===")
rows = []
for _ in range(400):
    eps = 10 ** rng.uniform(-2.2, -0.9)
    A = Im(rng.normal(size=8)); A /= np.linalg.norm(A)
    B = Im(rng.normal(size=8)); B /= np.linalg.norm(B)
    Cv = Im(rng.normal(size=8)); Cv /= np.linalg.norm(Cv)
    oe, of, og = oexp(eps*A), oexp(eps*B), oexp(eps*Cv)
    s1 = omul(of, oconj(oe)); s2 = omul(og, oconj(of)); s3 = omul(oe, oconj(og))
    order = np.linalg.norm(assoc(s1, s2, s3))
    H = frame_holonomy(oe, of, og); holo = np.linalg.norm(H - np.eye(8))
    loop = np.linalg.norm(readout(oe, of) + readout(of, og) + readout(og, oe))
    rows.append((sph_area(oe, of, og), order, holo, loop))
R = np.array(rows)
print(f"  {'量':<28}{'面積に対する傾き':>18}{'ホロノミーとの比':>18}")
print(f"  {'枠ホロノミー |H - I|':<28}{slope(R[:,0], R[:,2]):>18.3f}{1.0:>18.4f}")
print(f"  {'読み出しの一周':<28}{slope(R[:,0], R[:,3]):>18.3f}{np.median(R[:,3]/R[:,2]):>18.4f}")
print(f"  {'順序の差（結合子）':<28}{slope(R[:,0], R[:,1]):>18.3f}{np.median(R[:,1]/R[:,2]):>18.2e}")
print("  -> ホロノミーと同じ量なのは読み出しの一周（比 1/4）。順序の差は別の、より高次の量。")

# ---------------------------------------------------------------------------
print("\n=== 2. 内側の関係 conj(o_a) o_b は結合的な代数で telescoping する ===")
mx = 0.0
for _ in range(200):
    x, y, z = rng.normal(size=(3, 8))
    x[4:] = 0; y[4:] = 0; z[4:] = 0          # H（四元数部分代数）の中
    w1 = omul(oconj(x), y); w2 = omul(oconj(y), z); w3 = omul(oconj(z), x)
    P = omul(omul(w1, w2), w3)
    mx = max(mx, np.linalg.norm(Im(P))/np.linalg.norm(P))
print(f"  H の中：|Im(P)|/|P| の最大 = {mx:.1e}  -> 実数。位相は出ない。")
print("  外側の比 s_ab = o_b conj(o_a) で作った一周は telescoping しない。だから枠は回る。")

# ---------------------------------------------------------------------------
print("\n=== 3. 幾何学的位相は四元数の中でも出る（係数 4）===")
qmask = np.array([0, 1, 1, 1, 0, 0, 0, 0], float)
rec = []
for eps in np.logspace(-3.0, -1.5, 10):
    for _ in range(5):
        A = rng.normal(size=8)*qmask; A /= np.linalg.norm(A)
        B = rng.normal(size=8)*qmask; B /= np.linalg.norm(B)
        oe, of, og = E[0], oexp(eps*A), oexp(eps*B)
        H = frame_holonomy(oe, of, og); area = sph_area(oe, of, og)
        rec.append((np.linalg.norm(H[:4, :4] - np.eye(4)),
                    np.linalg.norm(H[:4, 4:]),
                    rot_angle(H[4:, 4:])/area))
rec = np.array(rec)
print(f"  時空ブロックの変化 |Hqq - I| の最大: {rec[:,0].max():.1e}（厳密に戻る）")
print(f"  時空 <-> 内部の混合の最大:          {rec[:,1].max():.1e}")
print(f"  内部の枠の回転角 / 面積（中央値）:   {np.median(rec[:,2]):.4f}  -> 4")

# ---------------------------------------------------------------------------
print("\n=== 4. 三つの一歩は一次で和ゼロ。だから eps^3 が消える ===")
rng2 = np.random.default_rng(7)
mx = 0.0
for _ in range(200):
    eps = 1e-3
    A = Im(rng2.normal(size=8)); A /= np.linalg.norm(A)
    B = Im(rng2.normal(size=8)); B /= np.linalg.norm(B)
    Cv = Im(rng2.normal(size=8)); Cv /= np.linalg.norm(Cv)
    oe, of, og = oexp(eps*A), oexp(eps*B), oexp(eps*Cv)
    d1 = omul(of, oconj(oe)) - E[0]; d2 = omul(og, oconj(of)) - E[0]; d3 = omul(oe, oconj(og)) - E[0]
    mx = max(mx, np.linalg.norm(d1 + d2 + d3)/eps)
print(f"  |d1 + d2 + d3| / eps の最大（eps = 1e-3）: {mx:.3e}  -> 一次で和ゼロ")
print("  結合子は交代形式なので [d1, d2, -d1-d2] = 0。したがって eps^3 の主要項が消える。")

print("\n=== 5. 実際の次数 ===")
eps_list = np.logspace(-3.0, -1.3, 14); med = []
for eps in eps_list:
    vals = []
    for _ in range(80):
        A = Im(rng2.normal(size=8)); A /= np.linalg.norm(A)
        B = Im(rng2.normal(size=8)); B /= np.linalg.norm(B)
        Cv = Im(rng2.normal(size=8)); Cv /= np.linalg.norm(Cv)
        oe, of, og = oexp(eps*A), oexp(eps*B), oexp(eps*Cv)
        s1 = omul(of, oconj(oe)); s2 = omul(og, oconj(of)); s3 = omul(oe, oconj(og))
        vals.append(np.linalg.norm(assoc(s1, s2, s3)))
    med.append(np.median(vals))
med = np.array(med)
p = np.polyfit(np.log(eps_list), np.log(med), 1)[0]
m = eps_list < 2e-2
p2 = np.polyfit(np.log(eps_list[m]), np.log(med[m]), 1)[0]
print(f"  結合子 ∝ eps^{p:.3f}   （eps < 2e-2 だけなら eps^{p2:.3f}）  -> 5 次")
print("  面積は eps^2 なので、面積に対しては 5/2 乗。")

print("\n=== まとめ（面積に対する次数）===")
print("  幾何学的位相（枠ホロノミー、係数 4） : 1 乗、H でも出る")
print("  読み出しの一周（ホロノミーの 1/4）   : 1 乗、H でも出る")
print("  八元数の補正（括り方の曖昧さ）       : 5/2 乗、H では厳密にゼロ")
print("  -> 八元数の効果は幾何学的位相より 3/2 次高く、小さい環では見えない。")
