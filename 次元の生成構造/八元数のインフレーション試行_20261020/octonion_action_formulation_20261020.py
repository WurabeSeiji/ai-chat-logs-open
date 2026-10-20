#!/usr/bin/env python3
"""
八元数の作用の定式化と、自己無撞着になる配置の探索空間（2026-10-20）

`次元の生成構造/八元数閉塞と関係の分解_20261019/octonion_closure_check_20261019.py`
を同じフォルダか一つ上の隣のフォルダに置いて実行する（代数の定義だけを読み込む）。

問い（著者）：高対称な八元数の初期値と、それが自己無撞着になりうる作用の組み合わせ・値を
定式化せよ。

定式化（本文は 作用の定式化と探索空間_20261020.md）
  1. 無名性・ノルム保存・積からの構成を課すと、作用は o_e' = u_e o_e、u_e in Im(O) に限る。
  2. 定理 A：G2 共変な双線形の作用は Sum o = 0 で厳密にゼロ（3 本すべて）。
  3. 生き残る最小の作用（第二層）：u_e = (1/2) Sum_f Im(w_ef^2)、w_ef = o_f conj(o_e)。
     駆動因子は |Im(w^2)| = |o_e|^2 |o_f|^2 |sin 2 theta_ef|。
  4. 定理 B：全対が 0/90/180 度、または各要素に直交する鏡映で集合が閉じるとき、
     この作用は厳密にゼロ（正軸体・E8 の根・24 胞体がこれに当たる）。
  5. 自己無撞着（持続する相対平衡）の条件：ある D in g2 で u_e o_e = D o_e（全 e）。
     そのとき o_e(tau) = exp(tau D) o_e(0)。D の虚 7 次元での三つの振動数は
     omega_1 + omega_2 + omega_3 = 0（符号付き）を満たし、
     ロック（U^n = I）には通約可能であることが要る。最小の比は 1 : 1 : -2。

出力：標準出力と results/action_formulation.json
再現性：乱数シード 20261020
"""
import json, os, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CAND = [os.path.join(HERE, "octonion_closure_check_20261019.py"),
        os.path.join(HERE, "..", "八元数閉塞と関係の分解_20261019", "octonion_closure_check_20261019.py")]
SRC = next((p for p in CAND if os.path.exists(p)), None)
if SRC is None:
    raise SystemExit("octonion_closure_check_20261019.py が見つかりません（同じフォルダに置いてください）")
exec(open(SRC).read().split("results = {")[0])          # omul, oconj, ojordan, E

SEED = 20261020
rng = np.random.default_rng(SEED)
OUT = os.path.join(HERE, "results"); os.makedirs(OUT, exist_ok=True)
out = {"seed": SEED}

def Im(x):
    y = np.asarray(x, float).copy(); y[0] = 0.0; return y

def nullbasis(A, tol=1e-9):
    u, s, vt = np.linalg.svd(A)
    return vt[int(np.sum(s > tol * s[0])):]

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
L = E[4]                                                 # 参照の向き ell
print(f"dim g2 = {len(g2)}, dim so(8) = {len(so8)}")

# ---------------------------------------------------------------------------
# 1. 双線形の分類と、ノルム保存
# ---------------------------------------------------------------------------
print("\n=== 1. G2 共変・双線形・純虚への写像は 3 本で尽きる ===")
e1 = max(np.linalg.norm(Im(omul(x, oconj(y)))
                        - (y[0]*Im(x) - x[0]*Im(y) - Im(omul(Im(x), Im(y)))))
         for x, y in rng.normal(size=(500, 2, 8)))
print(f"  Im(x conj(y)) = Re(y)Im(x) - Re(x)Im(y) - Im(x)xIm(y) の誤差: {e1:.1e}")
e2 = max(abs(float(omul(omul(Im(u), o), oconj(o))[0])) for u, o in rng.normal(size=(500, 2, 8)))
print(f"  Re((u o) conj(o)) = 0（u 純虚）の誤差: {e2:.1e}  -> |o_e| は厳密に保存")
out["bilinear_identity_err"] = float(e1); out["norm_preserving_err"] = float(e2)

# ---------------------------------------------------------------------------
# 2. 作用の基底
# ---------------------------------------------------------------------------
def actions(S):
    """各作用の u_e（N x 8）。すべて純虚。"""
    N = len(S)
    W = [[omul(S[f], oconj(S[e])) for f in range(N)] for e in range(N)]
    def acc(fn):
        U = np.zeros((N, 8))
        for e in range(N):
            s = np.zeros(8)
            for f in range(N):
                if f == e: continue
                s += fn(e, f, W[e][f])
            U[e] = s
        return U
    g = {}
    g["A1 Sum Im(w)"]            = acc(lambda e, f, w: Im(w))
    g["A2 Sum Re(w)Im(w)"]       = acc(lambda e, f, w: w[0]*Im(w))      # = (1/2) Sum Im(w^2)
    g["A3 Sum Re(w)^2 Im(w)"]    = acc(lambda e, f, w: w[0]**2*Im(w))
    g["A6 Sum Im(o_e)xIm(o_f)"]  = acc(lambda e, f, w: Im(omul(Im(S[e]), Im(S[f]))))
    g["B1 Sum<l o_e,o_f>Im(w)"]  = acc(lambda e, f, w: float(omul(L, S[e]) @ S[f])*Im(w))
    g["B2 u_e = l (uniform)"]    = np.tile(L, (N, 1))
    return g

def vel_elem(S, U):
    return np.array([omul(U[e], S[e]) for e in range(len(S))])

def vel_label(S):
    """番号の側：dot o_e = Sum_f K_ef o_f、K_ef = <l o_e, o_f>（反対称）"""
    N = len(S)
    K = np.array([[float(omul(L, S[e]) @ S[f]) for f in range(N)] for e in range(N)])
    return K @ S

def resid(S, V, rig):
    """min_Omega |V - Omega o| / |V|。Omega は rig（g2 か so(8)）の張る空間。"""
    A = np.column_stack([(S @ D.T).reshape(-1) for D in rig]); b = V.reshape(-1)
    nb = np.linalg.norm(b)
    if nb < 1e-10: return 0.0, None, True
    c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.linalg.norm(A @ c - b)/nb), sum(ci*D for ci, D in zip(c, rig)), False

def freqs(Om):
    w = np.sort(np.abs(np.linalg.eigvals(Om).imag))[::-1]; o = []
    for x in w:
        if x > 1e-9 and (not o or abs(x - o[-1]) > 1e-6*max(1, x)): o.append(x)
    return np.round(o, 6)

# ---------------------------------------------------------------------------
# 3. 駆動因子 sin 2theta
# ---------------------------------------------------------------------------
print("\n=== 3. 第二層の作用の駆動因子 ===")
mx = 0.0
for a, b in rng.normal(size=(500, 2, 8)):
    w = omul(b, oconj(a)); na, nb = np.linalg.norm(a), np.linalg.norm(b)
    th = np.arccos(np.clip((a @ b)/(na*nb), -1, 1))
    mx = max(mx, abs(np.linalg.norm(Im(omul(w, w))) - na*na*nb*nb*abs(np.sin(2*th))))
print(f"  |Im(w^2)| = |o_e|^2|o_f|^2 |sin 2theta| の誤差: {mx:.1e}")
print("  -> 0, 90, 180 度で厳密にゼロ、45 度で最大")
out["sin2theta_err"] = float(mx)

# ---------------------------------------------------------------------------
# 4. 配置
# ---------------------------------------------------------------------------
def angles_of(S):
    s = set()
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            c = (S[i] @ S[j])/(np.linalg.norm(S[i])*np.linalg.norm(S[j]))
            s.add(round(float(np.degrees(np.arccos(np.clip(c, -1, 1)))), 2))
    return sorted(s)

def refl_closed(S, tol=1e-9, probe=12):
    """各 o_e について、o_e に直交する鏡映で集合が閉じるか。"""
    for e in range(min(len(S), probe)):
        a = S[e]; n2 = a @ a
        for f in range(len(S)):
            y = S[f] - 2*(S[f] @ a)/n2*a
            if min(np.linalg.norm(y - S[g]) for g in range(len(S))) > tol:
                return False
    return True

cfgs = {}
cfgs["ORTHO16"] = np.vstack([np.vstack([E[i], -E[i]]) for i in range(8)])
nl = []
for i in range(1, 8):
    for s in (1, -1):
        v = (E[0] + s*E[i])/np.sqrt(2); nl.append(v); nl.append(-v)
cfgs["NULL28"] = np.array(nl)
c16 = np.zeros((8, 8)); sc = 1/np.sqrt(3); c16[0, 0] = 1; c16[1, 0] = -1
for i in range(1, 4): c16[2*i, i] = sc; c16[2*i+1, i] = -sc
cfgs["16-cell v/sqrt3"] = c16
d24 = []
for i in range(4):
    for s in (1, -1):
        v = np.zeros(8); v[i] = s; d24.append(v)
for sg in itertools.product([0.5, -0.5], repeat=4):
    v = np.zeros(8); v[:4] = sg; d24.append(v)
cfgs["24-cell (in H)"] = np.array(d24)
roots = []
for i in range(8):
    for j in range(i+1, 8):
        for si in (1, -1):
            for sj in (1, -1):
                v = np.zeros(8); v[i] = si; v[j] = sj; roots.append(v)
for s in itertools.product([0.5, -0.5], repeat=8):
    if np.prod([1 if x > 0 else -1 for x in s]) > 0: roots.append(np.array(s))
cfgs["E8 roots 240"] = np.array(roots)/np.sqrt(2)
def gen_closure(N, r, k=7):
    one = np.ones(N)/np.sqrt(N); q0 = r.normal(size=N); q0 -= one*(one @ q0)
    P = np.eye(N) - np.outer(one, one) - np.outer(q0, q0)/(q0 @ q0)
    C = np.array([P @ r.normal(size=N) for _ in range(k)]); C *= np.sqrt((q0 @ q0)/np.sum(C*C))
    bs = []
    for i in range(k):
        v = r.normal(size=8); v[0] = 0
        for b in bs: v -= (b @ v)*b
        bs.append(v/np.linalg.norm(v))
    return np.outer(q0, E[0]) + C.T @ np.array(bs)
cfgs["random closure N=9"] = gen_closure(9, np.random.default_rng(3))

# ---------------------------------------------------------------------------
# 5. 定理 B：第二層の作用が消える条件
# ---------------------------------------------------------------------------
print("\n=== 5. 定理 B：第二層の作用が消える条件 ===")
print(f"  {'config':<20}{'N':>5}  {'angles':<28}{'refl.closed':>12}{'max|u|':>12}")
tabB = []
for nm, S in cfgs.items():
    U = actions(S)["A2 Sum Re(w)Im(w)"]
    mu = max(np.linalg.norm(u) for u in U)
    rc = refl_closed(S); ang = str(angles_of(S))
    tabB.append({"name": nm, "N": len(S), "angles": angles_of(S),
                 "reflection_closed": bool(rc), "max_u_second_layer": float(mu)})
    print(f"  {nm:<20}{len(S):>5}  {ang[:28]:<28}{str(rc):>12}{mu:>12.2e}")
out["theorem_B"] = tabB
print("  -> 全対が 0/90/180 度、または鏡映で閉じる配置では厳密にゼロ。")

# ---------------------------------------------------------------------------
# 6. 各作用の自己無撞着性
# ---------------------------------------------------------------------------
print("\n=== 6. 各作用の相対平衡の残差と振動数 ===")
tabC = []
for nm, S in cfgs.items():
    if len(S) > 60: continue                              # E8 は重いので省く
    print(f"\n  --- {nm}（N = {len(S)}、|Sum o| = {np.linalg.norm(S.sum(0)):.0e}、"
          f"|Sum o^2| = {np.linalg.norm(sum(ojordan(v, v) for v in S)):.0e}）")
    print(f"    {'action':<26}{'res g2':>11}{'res so(8)':>12}   freqs(so(8))")
    G = actions(S)
    items = [(k, vel_elem(S, v)) for k, v in G.items()] + [("C1 label K=<l o_e,o_f>", vel_label(S))]
    for k, V in items:
        r2, _, z2 = resid(S, V, g2); r8, Om8, z8 = resid(S, V, so8)
        tag = "static" if z8 else (str(freqs(Om8)[:4]) if r8 < 1e-8 else "-")
        tabC.append({"config": nm, "action": k, "res_g2": r2, "res_so8": r8,
                     "static": bool(z8), "freqs": [] if z8 or r8 >= 1e-8 else [float(x) for x in freqs(Om8)[:4]]})
        print(f"    {k:<26}{r2:>11.2e}{r8:>12.2e}   {tag}")
out["relative_equilibria"] = tabC

# ---------------------------------------------------------------------------
# 7. g2 の振動数の拘束（符号付きで和ゼロ）
# ---------------------------------------------------------------------------
print("\n=== 7. D in g2 の虚 7 次元での三つの振動数は符号付きで和ゼロ ===")
print("  （固有値からは絶対値しか出ないので、最大 = 他の二つの和 で検定する）")
chk = []
for t in range(6):
    D = sum(c*M for c, M in zip(np.random.default_rng(SEED + t).normal(size=14), g2))
    im = np.sort(np.abs(np.linalg.eigvals(D).imag))[::-1]
    a = []
    for x in im:
        if x > 1e-9 and (not a or abs(x - a[-1]) > 1e-7*max(1, x)): a.append(x)
    a = a[:3]
    d = abs(a[0] - (a[1] + a[2])) / a[0]
    chk.append(d)
    print(f"  freqs {np.round(a, 6)}  |max - (sum of other two)| / max = {d:.1e}")
print("  -> omega_1 + omega_2 + omega_3 = 0（符号付き）。G2 の 7 次元表現の重みが")
print("     正六角形の一つおきの三つ組で和ゼロになることの帰結。")
print("     ロック（U^n = I）には三つが通約可能であることが要り、最小の比は 1 : 1 : -2。")
out["g2_frequency_sum_test_max_rel_dev"] = float(max(chk))

with open(os.path.join(OUT, "action_formulation.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("\nsaved:", os.path.join(OUT, "action_formulation.json"))
