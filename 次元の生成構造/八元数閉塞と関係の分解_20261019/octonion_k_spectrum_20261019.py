#!/usr/bin/env python3
"""
k の取りうる値の数え上げ（2026-10-19）

`octonion_closure_check_20261019.py` と同じフォルダに置いて実行する。

問い：粒子とは「ロックした配置」＝ k（配置が使っている虚方向の本数）が一般値より下に
落ちた配置ではないか。各 N でどの k が実現できるかを数え上げれば、粒子のスペクトルの
候補になる。

前提（これまでの結果）
  一般の閉塞解は k = min(7, N - 2)。
  k <= 2 ... 必ず四元数部分代数に収まる（内部ゼロ、電荷を持てない）。
  k = 3 ... 結合的（AxB が平面内）なら収まる。一般には収まらない。
  k >= 4 ... 次元で収まらず、内部が必ず開く。

このプログラムの中身
  1. 各 (N, k) で閉塞解が作れるか（構成）。
  2. 各 (N, k) の族の次元と G2xR 軌道・形の空間。
  3. k = 3 の族を結合的／非結合的に分ける。
  4. 等モジュラー（|o_n| 一定）を各 (N, k) に課したときに解が残るか。
  5. 既知の対称配置がどの (N, k) に落ちるか。
  6. k ごとに残る対称性（G2 の安定化群）。

再現性：乱数シード 20261019。出力は標準出力と results/k_spectrum.json。
"""
import json, os
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "octonion_closure_check_20261019.py")).read().split("results = {")[0])

SEED = 20261019
rng = np.random.default_rng(SEED)
OUT = os.path.join(HERE, "results"); os.makedirs(OUT, exist_ok=True)
out = {"seed": SEED}

def rank_(A, tol=1e-9):
    s = np.linalg.svd(A, compute_uv=False)
    return 0 if s.size == 0 or s[0] < 1e-300 else int(np.sum(s > tol * s[0]))

def closure_res(S):
    return (float(np.linalg.norm(S.sum(0))),
            float(np.linalg.norm(sum(ojordan(x, x) for x in S))))

def k_of(S, tol=1e-9):
    return rank_(S[:, 1:].T, tol)

def assoc_defect(S, tol=1e-9):
    """k = 3 のとき、3-平面が四元数部分代数の虚部か（AxB が平面内か）。"""
    U, sv, _ = np.linalg.svd(S[:, 1:].T)
    k = rank_(S[:, 1:].T, tol)
    if k <= 2: return 0.0
    if k >= 4: return 1.0
    B = np.array([np.concatenate([[0], U[:, i]]) for i in range(3)])
    w = omul(B[0], B[1]); w = w - w[0] * E[0]
    P = np.vstack(B)
    return float(np.linalg.norm(w - P.T @ (P @ w)) / np.linalg.norm(w))

# ---------------------------------------------------------------------------
# 1. 構成：虚方向をちょうど k 本に制限した閉塞解
#    k <= min(7, N - 2) は定理：Sigma o = 0 で虚部が 1 に直交、Sigma o^2 の虚部で q0 にも直交。
# ---------------------------------------------------------------------------
def gen_k(N, k, r, basis=None):
    """虚部がちょうど k 次元に収まる閉塞解。存在しないときは None。"""
    if N < 2 or k < 1 or k > min(7, N - 2):
        return None
    one = np.ones(N) / np.sqrt(N)
    q0 = r.normal(size=N); q0 -= one * (one @ q0)
    P = np.eye(N) - np.outer(one, one) - np.outer(q0, q0) / (q0 @ q0)
    C = np.array([P @ r.normal(size=N) for _ in range(k)])
    if rank_(C) < k: return None
    C *= np.sqrt((q0 @ q0) / np.sum(C * C))
    if basis is None:
        basis = []
        for _ in range(k):
            v = r.normal(size=8); v[0] = 0
            for b in basis: v -= (b @ v) * b
            n = np.linalg.norm(v)
            if n < 1e-9: return None
            basis.append(v / n)
        basis = np.array(basis)
    return np.outer(q0, E[0]) + C.T @ basis

def quaternion_basis(r, k):
    """四元数部分代数 span(1, A, B, AB) の虚部から k 本（k <= 3）。結合的 3-平面。"""
    A = r.normal(size=8); A[0] = 0; A /= np.linalg.norm(A)
    B = r.normal(size=8); B[0] = 0; B -= (A @ B) * A; B /= np.linalg.norm(B)
    C = omul(A, B); C[0] = 0; C /= np.linalg.norm(C)
    return np.array([A, B, C])[:k]

print("=== 1. 各 (N, k) で閉塞解が作れるか（構成による存在の確認）===")
print("  N \\ k    1    2    3    4    5    6    7     generic min(7,N-2)")
exist = {}
for N in range(2, 13):
    row = []
    for k in range(1, 8):
        ok = False
        for t in range(8):
            S = gen_k(N, k, np.random.default_rng(SEED + 31*N + 7*k + t))
            if S is None: continue
            rl, rs = closure_res(S)
            if rl < 1e-10 and rs < 1e-10 and k_of(S) == k: ok = True; break
        exist[(N, k)] = ok
        row.append(" o " if ok else " . ")
    print(f"  {N:>2}    " + "  ".join(row) + f"      {min(7, max(0, N-2))}")
out["exist"] = {f"{N},{k}": bool(v) for (N, k), v in exist.items()}
print("  o = 作れる（閉塞残差 1e-10 未満、虚部の階数がちょうど k）")
print("  -> k は 1 から min(7, N - 2) まで何でも取れる。N と k は独立なラベル。")

# ---------------------------------------------------------------------------
# 2. 各 (N, k) の族の次元
# ---------------------------------------------------------------------------
rows = []
for i in range(8):
    for j in range(8):
        prod = omul(E[i], E[j])
        for c in range(8):
            row = np.zeros(64)
            for m in range(8): row[c*8 + m] += prod[m]
            for m in range(8): row[m*8 + i] -= omul(E[m], E[j])[c]
            for m in range(8): row[m*8 + j] -= omul(E[i], E[m])[c]
            rows.append(row)
u, s, vt = np.linalg.svd(np.array(rows))
g2 = [v.reshape(8, 8) for v in vt[int(np.sum(s > 1e-9 * s[0])):]]
so7 = []
for i in range(1, 8):
    for j in range(i + 1, 8):
        D = np.zeros((8, 8)); D[i, j] = -1; D[j, i] = 1; so7.append(D)

def orbit_dim(S, alg):
    return rank_(np.array([S.reshape(-1)] + [(S @ D.T).reshape(-1) for D in alg]))

def family_dim(N, k):
    """q0：和零で N - 1。係数：和零かつ q0 直交で k(N - 2)、スケール 1 本が閉塞で固定。
       虚方向の基底：Stiefel V_k(R^7) で 7k - k(k+1)/2。"""
    return (N - 1) + k * (N - 2) - 1 + (7*k - k*(k+1)//2)

print("\n=== 2. 各 (N, k) の族の次元 ===")
print("  N   k   族の次元  G2xR 軌道  形（積も見る）")
dims = []
for N in [4, 5, 6, 7, 8]:
    for k in range(1, min(7, N-2) + 1):
        S = gen_k(N, k, np.random.default_rng(SEED + 91*N + 13*k))
        if S is None: continue
        o2 = orbit_dim(S, g2); o7 = orbit_dim(S, so7); fd = family_dim(N, k)
        dims.append({"N": N, "k": k, "family_dim_formula": fd, "orbit_G2": o2,
                     "shape_product": fd - o2, "orbit_SO7": o7, "shape_closure": fd - o7})
        print(f"  {N}   {k}      {fd:>3}        {o2:>3}          {fd-o2:>3}")
out["dims"] = dims

# ---------------------------------------------------------------------------
# 3. k = 3 を結合的／非結合的に分ける
# ---------------------------------------------------------------------------
print("\n=== 3. k = 3 の二つの族 ===")
k3 = []
for N in [5, 6, 7, 8]:
    r = np.random.default_rng(SEED + N)
    Sa = gen_k(N, 3, r, basis=quaternion_basis(r, 3))
    Sg = gen_k(N, 3, np.random.default_rng(SEED + 1000 + N))
    da, dg = assoc_defect(Sa), assoc_defect(Sg)
    ra, rg = closure_res(Sa), closure_res(Sg)
    k3.append({"N": N, "assoc_defect_associative": da, "assoc_defect_generic": dg})
    print(f"  N = {N}: associative {da:.1e}（閉塞 {ra[0]:.0e},{ra[1]:.0e}）/ "
          f"generic {dg:.3f}（閉塞 {rg[0]:.0e},{rg[1]:.0e}）")
out["k3"] = k3
print("  結合的 3-平面は Gr(3,7) 12 次元の中の G2/SO(4) 8 次元 -> 余次元 4。")

# ---------------------------------------------------------------------------
# 4. 等モジュラーを各 (N, k) に課す
# ---------------------------------------------------------------------------
def equal_modulus_in_k(N, k, seed, restarts=8):
    """|o_n| = 1、Sigma o = 0、Sigma o^2 = 0、かつ虚部が k 次元に収まる解を最小二乗で探す。"""
    def f(x):
        q0 = x[:N]; C = x[N:N + k*N].reshape(k, N); Bv = x[N + k*N:].reshape(k, 8)
        Bv = Bv.copy(); Bv[:, 0] = 0
        Q, _ = np.linalg.qr(Bv[:, 1:].T)
        Bo = np.zeros((k, 8)); Bo[:, 1:] = Q.T[:k]
        S = np.outer(q0, E[0]) + C.T @ Bo
        sq = sum(ojordan(v, v) for v in S)
        return (np.sum((np.sum(S**2, axis=1) - 1.0) ** 2) + np.sum(S.sum(0) ** 2) + float(sq @ sq))
    best = None
    for t in range(restarts):
        x0 = np.random.default_rng(seed * 100 + t).normal(size=N + k*N + k*8)
        res = minimize(f, x0, method="L-BFGS-B",
                       options=dict(maxiter=20000, maxfun=60000, ftol=1e-18, gtol=1e-14))
        if best is None or res.fun < best.fun: best = res
    return best.fun

print("\n=== 4. 等モジュラー（|o_n| 一定）を (N, k) に課す ===")
print("  N \\ k      1        2        3        4        5")
em = {}
for N in range(3, 9):
    cells = []
    for k in range(1, 6):
        if k > min(7, N - 2): cells.append("   --   "); continue
        err = equal_modulus_in_k(N, k, SEED + 7*N + k)
        em[(N, k)] = err
        cells.append(f" {err:.0e} " if err > 1e-12 else "  solved ")
    print(f"  {N:>2}    " + " ".join(cells))
out["equal_modulus_residual"] = {f"{N},{k}": float(v) for (N, k), v in em.items()}
print("  残差が 1e-12 未満なら等モジュラー解が存在する。")
print("  注：k を拘束せずに等モジュラーを解くと一般値 k = min(7, N-2) に落ちる。")
print("      k <= N - 2 は定理（Sigma o = 0 で虚部が 1 に直交、Sigma o^2 の虚部で q0 にも直交）。")

# ---------------------------------------------------------------------------
# 5. 既知の対称配置
# ---------------------------------------------------------------------------
print("\n=== 5. 既知の対称配置 ===")
known = []
def reg(name, S):
    k = k_of(S); rl, rs = closure_res(S); d = assoc_defect(S)
    gk = min(7, max(0, len(S) - 2))
    known.append({"name": name, "N": len(S), "k": k, "k_generic": gk,
                  "assoc_defect": d, "res_lin": rl, "res_sq": rs})
    print(f"  {name}: N = {len(S)}, k = {k}（generic {gk}）, assoc defect = {d:.1e}, "
          f"closure {rl:.0e},{rs:.0e}")

cross = np.zeros((4, 8)); cross[0,0]=1; cross[1,0]=-1; cross[2,1]=1; cross[3,1]=-1
reg("cross (1,-1,i,-i)    ", cross)
w = np.exp(2j*np.pi/3); tri = np.zeros((3, 8))
for n in range(3):
    z = w**n; tri[n,0] = z.real; tri[n,1] = z.imag
reg("three-phase (1,w,w^2)", tri)
c16 = np.zeros((8, 8)); sc = 1/np.sqrt(3)
c16[0,0]=1; c16[1,0]=-1
for i in range(1, 4):
    c16[2*i, i] = sc; c16[2*i+1, i] = -sc
reg("16-cell (v/sqrt3)    ", c16)
two = np.zeros((2, 8)); two[0,0]=1; two[1,1]=1
rl, rs = closure_res(two)
print(f"  two-wave (1, i)      : N = 2, k = {k_of(two)}, square closure {rs:.0e}, "
      f"linear closure {rl:.1f}（満たさない）")
out["known"] = known

# ---------------------------------------------------------------------------
# 6. k ごとに残る対称性（G2 の安定化群）
# ---------------------------------------------------------------------------
def nulldim(A, tol=1e-9):
    sv = np.linalg.svd(A, compute_uv=False)
    return int(np.sum(sv < tol * max(1.0, sv[0]))) + (A.shape[1] - len(sv))

def stab_dim(dirs):
    C = []
    for d in dirs:
        for c in range(8):
            C.append([(D @ d)[c] for D in g2])
    return nulldim(np.array(C))

r = np.random.default_rng(5)
A = r.normal(size=8); A[0] = 0; A /= np.linalg.norm(A)
B = r.normal(size=8); B[0] = 0; B -= (A @ B) * A; B /= np.linalg.norm(B)
Cq = omul(A, B); Cq[0] = 0; Cq /= np.linalg.norm(Cq)
Dg = r.normal(size=8); Dg[0] = 0
for v in (A, B, Cq): Dg -= (v @ Dg) * v
Dg /= np.linalg.norm(Dg)

print("\n=== 6. k ごとに残る対称性（G2 の安定化群の次元）===")
lad = {"k=1": stab_dim([A]), "k=2": stab_dim([A, B]),
       "k=3 associative": stab_dim([A, B, Cq]), "k=3 generic": stab_dim([A, B, Dg])}
print(f"  k = 1（1 本固定）                     : {lad['k=1']}   su(3)")
print(f"  k = 2（2 本固定）                     : {lad['k=2']}   su(2)")
print(f"  k = 3 associative（A, B, AxB を固定） : {lad['k=3 associative']}   su(2)（3 本目は A, B から決まる）")
print(f"  k = 3 generic（A, B, 別方向を固定）   : {lad['k=3 generic']}   自明")
print("  k >= 4 も自明（3 本で自明になった後は減らない）")
out["residual_symmetry"] = lad

with open(os.path.join(OUT, "k_spectrum.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("\nsaved:", os.path.join(OUT, "k_spectrum.json"))
