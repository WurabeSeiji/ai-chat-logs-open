#!/usr/bin/env python3
"""
作用の定式化と、自己無撞着になる配置の条件（2026-10-20）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか、環境変数 OCTO_BASE の指すフォルダに置いて実行する
（八元数の積の定義を読み込む。親フォルダ経由のフォールバックも入れてある）。

問い（著者、2026-10-20）：高対称な八元数の初期値と、それが自己無撞着になりうる作用の
組み合わせと値を定式化できるか。

このプログラムが示すこと
  1. 許される作用の形：o_e の大きさを厳密に保つ作用は o'_e = u_e o_e、u_e in Im(O) に限る。
  2. 定理 A：G2 共変で双線形な作用（3 本しかない）は Sum o = 0 で厳密にゼロになる。
  3. 生き残る最小の作用（第二層）：u_e = (1/2) Sum Im(w_ef^2)、w_ef = o_f conj(o_e)。
     駆動因子は |Im(w^2)| = |o_e|^2 |o_f|^2 |sin 2 theta_ef|。
  4. 定理 B：高対称な配置が静止する二つの仕組み
       (a) 全対の角が 0/90/180 度 -> sin 2theta = 0（正軸体）
       (b) 各要素に直交する鏡映で集合が閉じる -> 60 度の寄与が対で打ち消される（E8、24 胞体）
  5. 自己無撞着の条件：∃ D in g2 で u_e o_e = D o_e（全 e）。D は自己同型なので条件が持続する。
     g2 の重みは正六角形なので三振動数は和ゼロ。ロックには通約可能、最小は 1 : 1 : 2。
  6. 各配置・各作用の残差と振動数の表。

再現性：固定シード。出力は標準出力のみ。
"""
import os, itertools
import numpy as np

BASE = os.environ.get("OCTO_BASE", os.path.dirname(os.path.abspath(__file__)))
src = os.path.join(BASE, "octonion_closure_check_20261019.py")
if not os.path.exists(src):
    src = os.path.join(BASE, "..", "八元数閉塞と関係の分解_20261019", "octonion_closure_check_20261019.py")
exec(open(src).read().split("results = {")[0])

rng = np.random.default_rng(20261020)
def Im(x):
    y = np.asarray(x, float).copy(); y[0] = 0.0; return y
def nullbasis(A, tol=1e-9):
    u, s, vt = np.linalg.svd(A); return vt[int(np.sum(s > tol * s[0])):]

# --- g2 と so(8) -----------------------------------------------------------
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
so8 = []
for i in range(8):
    for j in range(i+1, 8):
        D = np.zeros((8, 8)); D[i, j] = -1; D[j, i] = 1; so8.append(D)
L = E[4]                                   # 参照の向き l

print("=== 1. ノルムを保つ作用の形 ===")
print("  o'_e = u_e o_e、u_e in Im(O) なら d|o_e|^2/dt = 2 Re((u_e o_e) conj(o_e)) = 2 Re(u_e |o_e|^2) = 0")
print(f"  Re((u o) conj(o)) = 0 の誤差: "
      f"{max(abs(float(omul(omul(Im(a), b), oconj(b))[0])) for a, b in rng.normal(size=(500, 2, 8))):.1e}")

print("\n=== 2. G2 共変な双線形写像は 3 本（1x7, 7x1, 7x7）===")
mx = 0.0
for _ in range(500):
    x, y = rng.normal(size=(2, 8))
    mx = max(mx, np.linalg.norm(Im(omul(x, oconj(y)))
                                - (y[0]*Im(x) - x[0]*Im(y) - Im(omul(Im(x), Im(y))))))
print(f"  Im(x conj(y)) = Re(y)Im(x) - Re(x)Im(y) - Im(x)xIm(y) の誤差: {mx:.1e}")

def gen_closure(N, r, k=7):
    one = np.ones(N)/np.sqrt(N); q0 = r.normal(size=N); q0 -= one*(one @ q0)
    P = np.eye(N) - np.outer(one, one) - np.outer(q0, q0)/(q0 @ q0)
    C = np.array([P @ r.normal(size=N) for _ in range(k)]); C *= np.sqrt((q0 @ q0)/np.sum(C*C))
    bs = []
    for _ in range(k):
        v = r.normal(size=8); v[0] = 0
        for b in bs: v -= (b @ v)*b
        bs.append(v/np.linalg.norm(v))
    return np.outer(q0, E[0]) + C.T @ np.array(bs)

print("\n=== 定理 A：双線形の作用は Sum o = 0 で厳密にゼロ ===")
S = gen_closure(9, np.random.default_rng(1))
def u_first(S):
    return np.array([Im(sum(omul(S[f], oconj(S[e])) for f in range(len(S)) if f != e))
                     for e in range(len(S))])
def u_cross(S):
    return np.array([sum(Im(omul(Im(S[e]), Im(S[f]))) for f in range(len(S)) if f != e)
                     for e in range(len(S))])
print(f"  試験配置 N=9: |Sum o| = {np.linalg.norm(S.sum(0)):.0e}, "
      f"|Sum o^2| = {np.linalg.norm(sum(ojordan(v, v) for v in S)):.0e}")
print(f"  Sum Im(w)             の最大ノルム: {max(np.linalg.norm(u) for u in u_first(S)):.1e}")
print(f"  Sum Im(o_e) x Im(o_f) の最大ノルム: {max(np.linalg.norm(u) for u in u_cross(S)):.1e}")
print("  理由: Sum_{f != e} o_f = -o_e。Im(-o_e conj(o_e)) = 0、Im(o_e) x Im(o_e) = 0。")

print("\n=== 3. 生き残る最小の作用（第二層）と駆動因子 ===")
print("  u_e = (1/2) Sum_{f != e} Im(w_ef^2) = Sum Re(w_ef) Im(w_ef)、w_ef = o_f conj(o_e)")
mx = 0.0
for _ in range(500):
    a, b = rng.normal(size=(2, 8))
    w = omul(b, oconj(a)); na, nb = np.linalg.norm(a), np.linalg.norm(b)
    th = np.arccos(np.clip((a @ b)/(na*nb), -1, 1))
    mx = max(mx, abs(np.linalg.norm(Im(omul(w, w))) - na*na*nb*nb*abs(np.sin(2*th))))
print(f"  |Im(w^2)| = |o_e|^2 |o_f|^2 |sin 2theta| の誤差: {mx:.1e}")
print("  -> 0, 90, 180 度で厳密にゼロ、45 度で最大")

def u_second(S):
    N = len(S); U = np.zeros((N, 8))
    for e in range(N):
        s = np.zeros(8); oe = oconj(S[e])
        for f in range(N):
            if f == e: continue
            w = omul(S[f], oe); s += w[0]*Im(w)
        U[e] = s
    return U

# --- 配置 -------------------------------------------------------------------
cfg = {}
cfg["ORTHO16"] = np.vstack([np.vstack([E[i], -E[i]]) for i in range(8)])
nl = []
for i in range(1, 8):
    for s in (1, -1):
        v = (E[0] + s*E[i])/np.sqrt(2); nl.append(v); nl.append(-v)
cfg["NULL28"] = np.array(nl)
c16 = np.zeros((8, 8)); sc = 1/np.sqrt(3); c16[0, 0] = 1; c16[1, 0] = -1
for i in range(1, 4): c16[2*i, i] = sc; c16[2*i+1, i] = -sc
cfg["16-cell v/sqrt3"] = c16
roots = []
for i in range(8):
    for j in range(i+1, 8):
        for si in (1, -1):
            for sj in (1, -1):
                v = np.zeros(8); v[i] = si; v[j] = sj; roots.append(v)
for s in itertools.product([0.5, -0.5], repeat=8):
    if np.prod([1 if x > 0 else -1 for x in s]) > 0: roots.append(np.array(s))
cfg["E8 roots 240"] = np.array(roots)/np.sqrt(2)
d24 = []
for i in range(4):
    for s in (1, -1):
        v = np.zeros(8); v[i] = s; d24.append(v)
for sg in itertools.product([0.5, -0.5], repeat=4):
    v = np.zeros(8); v[:4] = sg; d24.append(v)
cfg["24-cell in H"] = np.array(d24)
cfg["random closure N=9"] = gen_closure(9, np.random.default_rng(3))

def angles(S):
    out = set()
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            c = (S[i] @ S[j])/(np.linalg.norm(S[i])*np.linalg.norm(S[j]))
            out.add(round(float(np.degrees(np.arccos(np.clip(c, -1, 1)))), 1))
    return sorted(out)
def refl_closed(S, tol=1e-9, cap=12):
    for e in range(min(len(S), cap)):
        a = S[e]; n2 = a @ a
        for f in range(len(S)):
            y = S[f] - 2*(S[f] @ a)/n2*a
            if min(np.linalg.norm(y - S[g]) for g in range(len(S))) > tol: return False
    return True

print("\n=== 定理 B：高対称な配置で第二層の作用が消える二つの仕組み ===")
print(f"  {'配置':<20}{'N':>5}  {'角度':<26}{'鏡映で閉じる':>12}{'max |u|':>12}")
for nm, S in cfg.items():
    U = u_second(S)
    print(f"  {nm:<20}{len(S):>5}  {str(angles(S))[:26]:<26}"
          f"{str(refl_closed(S)):>12}{max(np.linalg.norm(u) for u in U):>12.3e}")
print("  (a) 全対が 0/90/180 度 -> sin 2theta = 0（ORTHO16）")
print("  (b) 各要素に直交する鏡映で閉じる -> 60 度の寄与が対で打ち消される（E8、24-cell）")
print("  NULL28 だけが鏡映で閉じず、u != 0 になる。")

# --- 自己無撞着の判定 --------------------------------------------------------
def resid(S, V, rig):
    A = np.column_stack([(S @ D.T).reshape(-1) for D in rig]); b = V.reshape(-1)
    nb = np.linalg.norm(b)
    if nb < 1e-10: return 0.0, None, True
    c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.linalg.norm(A @ c - b)/nb), sum(ci*D for ci, D in zip(c, rig)), False
def freqs(Om):
    w = np.sort(np.abs(np.linalg.eigvals(Om).imag))[::-1]; out = []
    for x in w:
        if x > 1e-9 and (not out or abs(x - out[-1]) > 1e-6*max(1, x)): out.append(x)
    return np.round(out, 5)

def actions(S):
    N = len(S); W = [[omul(S[f], oconj(S[e])) for f in range(N)] for e in range(N)]
    def acc(fn):
        U = np.zeros((N, 8))
        for e in range(N):
            s = np.zeros(8)
            for f in range(N):
                if f == e: continue
                s += fn(e, f, W[e][f])
            U[e] = s
        return U
    a = {}
    a["A1 Sum Im(w)"]            = acc(lambda e, f, w: Im(w))
    a["A2 Sum Re(w)Im(w)"]       = acc(lambda e, f, w: w[0]*Im(w))
    a["A3 Sum Re(w)^2 Im(w)"]    = acc(lambda e, f, w: w[0]**2*Im(w))
    a["A6 Sum Im(oe)xIm(of)"]    = acc(lambda e, f, w: Im(omul(Im(S[e]), Im(S[f]))))
    a["B1 Sum <l oe,of> Im(w)"]  = acc(lambda e, f, w: float(omul(L, S[e]) @ S[f])*Im(w))
    a["B2 u_e = l (uniform)"]    = np.tile(L, (N, 1))
    return a
def label_motion(S):
    N = len(S)
    K = np.array([[float(omul(L, S[e]) @ S[f]) for f in range(N)] for e in range(N)])
    return K @ S

print("\n=== 4. 自己無撞着の判定：∃ D で u_e o_e = D o_e ===")
print("  g2 の残差が 0 なら条件は時間発展で持続する（D は自己同型）。")
print("  so(8) だけ 0 のときは形は保つが、持続は別に確かめる必要がある。")
for nm in ["ORTHO16", "NULL28", "16-cell v/sqrt3", "random closure N=9"]:
    S = cfg[nm]
    print(f"\n  --- {nm} (N = {len(S)}) ---")
    print(f"      {'作用':<26}{'g2 残差':>11}{'so(8) 残差':>12}   振動数(so(8))")
    items = [(k, np.array([omul(v[e], S[e]) for e in range(len(S))]))
             for k, v in actions(S).items()]
    items.append(("C1 label side K=<l oe,of>", label_motion(S)))
    for k, V in items:
        r2, _, z2 = resid(S, V, g2); r8, Om8, z8 = resid(S, V, so8)
        tag = "static" if z8 else (str(freqs(Om8)[:4]) if r8 < 1e-8 else "-")
        print(f"      {k:<26}{r2:>11.2e}{r8:>12.2e}   {tag}")

print("\n=== 5. g2 の中の回転は三振動数が和ゼロ ===")
X = sum(c*D for c, D in zip(rng.normal(size=14), g2))
ev = np.linalg.eigvals(X).imag; w = np.sort(ev[ev > 1e-9])
print(f"  一般の D in g2 の正の振動数: {np.round(w, 6)}  "
      f"最大 - 他の二つの和 = {w.max() - (w.sum() - w.max()):.1e}")
print("  -> omega1 + omega2 = omega3。ロック U^n = I には通約可能が要り、最小は 1 : 1 : 2。")
