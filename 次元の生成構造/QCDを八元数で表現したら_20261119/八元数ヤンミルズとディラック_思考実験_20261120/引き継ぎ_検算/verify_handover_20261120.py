"""思考実験(発言録 第1〜15部)に出てくる式の検算。2026-11-20
発言録の番号 (## n) を [n] として示す。乱数の八元数で恒等式を確かめる。
"""
import itertools
import numpy as np

rng = np.random.default_rng(20261120)
TRIPLES = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
f = np.zeros((8, 8, 8))
for (i, j, k) in TRIPLES:
    for (a, b, c) in [(i, j, k), (j, k, i), (k, i, j)]:
        f[a, b, c] = 1
        f[b, a, c] = -1
M = np.zeros((8, 8, 8))  # (x y)_k = M[i,j,k] x_i y_j
for i in range(8):
    M[0, i, i] = 1
    M[i, 0, i] = 1
for i in range(1, 8):
    M[i, i, 0] = -1
M[1:, 1:, 1:] += f[1:, 1:, 1:]


def mul(x, y):
    return np.einsum('ijk,i,j->k', M, x, y)


def conj(x):
    y = -x.copy()
    y[0] = x[0]
    return y


def rnd(unit=False):
    x = rng.normal(size=8)
    return x / np.linalg.norm(x) if unit else x


def cross(b, c):
    return np.einsum('ijk,i,j->k', f[1:, 1:, 1:], b, c)


def phi(b1, b2, b3):
    return np.einsum('ijk,i,j,k->', f[1:, 1:, 1:], b1, b2, b3)


def power(x, m):
    y = np.zeros(8)
    y[0] = 1
    for _ in range(m):
        y = mul(y, x)
    return y


results = []


def check(tag, name, err, tol=1e-9):
    ok = err < tol
    results.append(ok)
    print(f"[{tag}] {'成立' if ok else '不成立'}  {name}   最大誤差 = {err:.2e}")


# ---- [7][8] 三つの積の実部
e1 = e2 = e3 = e4 = e5 = 0.0
for _ in range(500):
    x, y, z = rnd(), rnd(), rnd()
    r1 = mul(mul(x, y), z)[0]
    e1 = max(e1, abs(r1 - mul(x, mul(y, z))[0]))
    e2 = max(e2, abs(r1 - mul(mul(y, z), x)[0]))
    a1, a2, a3 = x[0], y[0], z[0]
    b1, b2, b3 = x[1:], y[1:], z[1:]
    S = a1 * a2 * a3 - a1 * (b2 @ b3) - a2 * (b1 @ b3) - a3 * (b1 @ b2) - phi(b1, b2, b3)
    e3 = max(e3, abs(r1 - S))
    # 二つの積の成分表示
    xy = mul(x, y)
    e4 = max(e4, abs(xy[0] - (a1 * a2 - b1 @ b2)) + np.abs(xy[1:] - (a1 * b2 + a2 * b1 + cross(b1, b2))).max())
    # ∂S/∂o1 = conj(o2 o3):  Re(x w) = <x, conj(w)>
    w = mul(y, z)
    e5 = max(e5, abs(mul(x, w)[0] - x @ conj(w)))
check(7, "Re((xy)z) = Re(x(yz))  括弧によらない", e1)
check(7, "Re(xyz) = Re(yzx)  巡回的", e2)
check(8, "Re(o1o2o3) = a1a2a3 − a1(b2·b3) − a2(b1·b3) − a3(b1·b2) − φ(b1,b2,b3)", e3)
check(8, "xy = (a_x a_y − b_x·b_y) + (a_x b_y + a_y b_x + b_x×b_y)", e4)
check(8, "Re(x w) = <x, conj(w)>  (変分が conj(o2o3) になる)", e5)

# ---- [8] 三つ組が一つだと一手で止まる
e1 = e2 = 0.0
for _ in range(200):
    o2, o3 = rnd(True), rnd(True)
    o1n = conj(mul(o2, o3))
    e1 = max(e1, abs(mul(mul(o1n, o2), o3)[0] - 1))
    o2n = conj(mul(o3, o1n))
    e2 = max(e2, np.abs(o2n - o2).max())
check(8, "o1' = conj(o2o3) とすると Re(o1'o2o3) = 1", e1)
check(8, "続けて o2' = conj(o3 o1') は o2 に等しい(固定点)", e2)

# ---- [6] 二元が生成する部分代数は結合的 (アルティン)
e1 = 0.0
for _ in range(200):
    x, y = rnd(), rnd()
    e1 = max(e1, np.abs(mul(mul(x, y), x) - mul(x, mul(y, x))).max(),
             np.abs(mul(mul(x, x), y) - mul(x, mul(x, y))).max(),
             np.abs(mul(mul(y, x), x) - mul(y, mul(x, x))).max())
check(6, "交代性 (xx)y = x(xy), (xy)x = x(yx), (yx)x = y(xx)", e1)

# ---- [5] 𝔉² = −ΣF_iF_i + ½ f_ijk [F_i,F_j] e_k  (F_i は 3×3 行列)
F = [rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)) for _ in range(7)]
sq = [np.zeros((3, 3), complex) for _ in range(8)]
for i in range(1, 8):
    for j in range(1, 8):
        for k in range(8):
            if M[i, j, k] != 0:
                sq[k] += M[i, j, k] * F[i - 1] @ F[j - 1]
err = np.abs(sq[0] + sum(Fi @ Fi for Fi in F)).max()
for k in range(1, 8):
    t = np.zeros((3, 3), complex)
    for i in range(1, 8):
        for j in range(1, 8):
            t += 0.5 * f[i, j, k] * (F[i - 1] @ F[j - 1] - F[j - 1] @ F[i - 1])
    err = max(err, np.abs(sq[k] - t).max())
check(5, "𝔉² = −ΣF_iF_i + ½ f_ijk[F_i,F_j]e_k", err)

# ---- [12] 左掛け算はクリフォード代数
L = [np.einsum('ijk,i->kj', M, np.eye(8)[i]) for i in range(8)]
err = 0.0
for i in range(1, 8):
    for j in range(1, 8):
        err = max(err, np.abs(L[i] @ L[j] + L[j] @ L[i] + 2 * (i == j) * np.eye(8)).max())
check(12, "{L_ei, L_ej} = −2δ_ij", err)
err = max(np.abs(L[1] @ L[2] - np.einsum('ijk,i->kj', M, mul(np.eye(8)[1], np.eye(8)[2]))).max(), 0)
print(f"      (参考) L_e1 L_e2 と L_(e1e2) の差の最大 = {err:.2f}  → 0 でない = L_x L_y ≠ L_xy")

# ---- [3][4] e_ρ(ē_μ e_ν) の分解と ½Φ の固有値
E = np.eye(8)
chi = np.zeros((8, 8, 8, 8))
err = 0.0
for r, m, n in itertools.product(range(8), repeat=3):
    v = mul(E[r], mul(conj(E[m]), E[n]))
    base = (r == m) * E[n] - (r == n) * E[m] + (m == n) * E[r]
    chi[r, m, n] = v - base
for perm_sign, perm in [(-1, (1, 0, 2, 3)), (-1, (0, 2, 1, 3)), (-1, (0, 1, 3, 2))]:
    err = max(err, np.abs(chi - perm_sign * np.transpose(chi, perm)).max())
check(3, "e_ρ(ē_μ e_ν) = δ_ρμ e_ν − δ_ρν e_μ + δ_μν e_ρ + χ_ρμνσ e_σ、χ は完全反対称", err)
pairs = [(m, n) for m in range(8) for n in range(m + 1, 8)]
T = np.zeros((28, 28))
for a, (m, n) in enumerate(pairs):
    for b, (r, s) in enumerate(pairs):
        T[a, b] = chi[m, n, r, s]  # ½Φ_μνρσ F_ρσ を ρ<σ の和で書いたもの
ev = np.round(np.linalg.eigvalsh(T), 6)
vals, cnt = np.unique(ev, return_counts=True)
print(f"      ½Φ· の固有値と重複度 : {dict(zip(vals.tolist(), cnt.tolist()))}")
ok = sorted(cnt.tolist()) == [7, 21] and abs(abs(vals[np.argmin(cnt)]) - 3) < 1e-6 and abs(abs(vals[np.argmax(cnt)]) - 1) < 1e-6
results.append(ok)
print(f"[4] {'成立' if ok else '不成立'}  28 = 7 + 21、固有値の大きさは 7 成分で 3、21 成分で 1 (符号は規約による)")

# ---- [15] 無名な和は Re(P³) に潰れる
N = 7
os_ = [rnd(True) for _ in range(N)]
S = sum(mul(mul(a, b), c)[0] for a in os_ for b in os_ for c in os_)
P = sum(os_)
a, b = P[0], P[1:]
check(15, "Σ_{A,B,C} Re(o_A o_B o_C) = Re(P³) = a³ − 3a|b|²",
      max(abs(S - mul(mul(P, P), P)[0]), abs(S - (a ** 3 - 3 * a * (b @ b)))))
Rr = np.linalg.norm(P)
th = np.arccos(a / Rr)
check(15, "Re(P³) = R³ cos 3θ", abs(S - Rr ** 3 * np.cos(3 * th)))
S2 = sum(mul(mul(a_, b_), c_)[0] ** 2 for a_ in os_ for b_ in os_ for c_ in os_)
print(f"      (参考) S₂ = Σ[Re(ooo)]² = {S2:.4f}  と  [Re(P³)]² = {S**2:.4f}  → 一致しない = P の関数に潰れない")

# ---- [16] 三角法: |o_A o_B − conj(o_C)|² = 2 − 2 Re(o_A o_B o_C)
e1 = e2 = 0.0
for _ in range(300):
    x, y, z = rnd(True), rnd(True), rnd(True)
    e1 = max(e1, abs(np.sum((x - y) ** 2) - (2 - 2 * mul(x, conj(y))[0])))
    e2 = max(e2, abs(np.sum((mul(x, y) - conj(z)) ** 2) - (2 - 2 * mul(mul(x, y), z)[0])))
check(16, "|x − y|² = 2 − 2 Re(x ȳ)", e1)
check(16, "|o_A o_B − ō_C|² = 2 − 2 Re(o_A o_B o_C)", e2)

# ---- [26][30][32] 倍音は冪、O² の実部と虚部
e1 = e2 = e3 = 0.0
for _ in range(200):
    n = rng.normal(size=7)
    n /= np.linalg.norm(n)
    A, th = rng.uniform(0.5, 2), rng.uniform(0, np.pi)
    O = np.concatenate([[A * np.cos(th)], A * np.sin(th) * n])
    for m in (2, 3, 5):
        Om = power(O, m)
        e1 = max(e1, np.abs(Om - A ** m * np.concatenate([[np.cos(m * th)], np.sin(m * th) * n])).max())
    O2 = mul(O, O)
    e2 = max(e2, abs(O2[0] - (O[0] ** 2 - O[1:] @ O[1:])) + np.abs(O2[1:] - 2 * O[0] * O[1:]).max())
    # 向きの違う二つの m 次倍音の干渉の実部
    n2 = rng.normal(size=7)
    n2 /= np.linalg.norm(n2)
    th2 = rng.uniform(0, np.pi)
    psi = np.arccos(n @ n2)
    for m in (1, 2, 4):
        Am = np.concatenate([[np.cos(m * th)], np.sin(m * th) * n])
        Bm = np.concatenate([[np.cos(m * th2)], np.sin(m * th2) * n2])
        val = mul(Am, conj(Bm))[0]
        e3 = max(e3, abs(val - (np.cos(psi / 2) ** 2 * np.cos(m * (th - th2)) + np.sin(psi / 2) ** 2 * np.cos(m * (th + th2)))))
check(26, "Oᵐ = Aᵐ(cos mθ + n sin mθ)", e1)
check(32, "O² = (a² − |b|²) + 2ab", e2)
check(30, "向きの違う倍音の干渉の実部 = cos²(ψ/2)cos m(θA−θB) + sin²(ψ/2)cos m(θA+θB)", e3)

# ---- [26][33] 倍音の和の公式
e1 = e2 = 0.0
for _ in range(200):
    Nn = int(rng.integers(2, 30))
    th = rng.uniform(0.01, np.pi - 0.01)
    s = sum(np.cos(k * th) for k in range(1, Nn + 1))
    e1 = max(e1, abs(s - np.cos((Nn + 1) * th / 2) * np.sin(Nn * th / 2) / np.sin(th / 2)))
    re = sum(np.cos(2 * k * th) for k in range(1, Nn + 1))
    im = sum(np.sin(2 * k * th) for k in range(1, Nn + 1))
    pre = np.sin(Nn * th) / np.sin(th)
    e2 = max(e2, abs(re - pre * np.cos((Nn + 1) * th)) + abs(im - pre * np.sin((Nn + 1) * th)))
check(26, "Σ_{n=1}^N cos nθ = cos((N+1)θ/2)·sin(Nθ/2)/sin(θ/2)", e1)
check(33, "Σ_n o_kn² = a²·sin(Nθ)/sinθ·[cos((N+1)θ) + n_k sin((N+1)θ)]", e2)

# ---- [36] 二次までの展開
Nn, th = 200, 1e-4
re = sum(np.cos(2 * k * th) for k in range(1, Nn + 1)) / Nn
im = sum(np.sin(2 * k * th) for k in range(1, Nn + 1)) / Nn
print(f"[36] N={Nn}, θ={th}:  虚部/P = {im:.6e}  (式 (N+1)θ = {(Nn+1)*th:.6e})")
print(f"      1−実部/P = {1-re:.6e}  (式 (2/3)N²θ² = {2/3*Nn**2*th**2:.6e}、N が大きいときの値)")
ok = abs(im / ((Nn + 1) * th) - 1) < 1e-3 and abs((1 - re) / (2 / 3 * Nn ** 2 * th ** 2) - 1) < 2e-2
results.append(ok)
print(f"[36] {'成立' if ok else '不成立'}  R_k² ≈ P_k[1 + n_k(N+1)θ − (2/3)N²θ²]")

# ---- [40][41] 四元数二つに分けた積 (ℍ = <1,e1,e2,e3>、ℓ = e7)
def split(o):  # o = q1 + q2 ℓ となる q1, q2 (どちらも ℍ の元、8 成分で持つ)
    q1 = np.zeros(8)
    q1[:4] = o[:4]
    rest = o - q1
    q2 = -mul(rest, E[7])  # (q2 ℓ) ℓ = −q2
    return q1, q2


e0 = e1 = e2 = 0.0
for _ in range(300):
    x, y = rnd(), rnd()
    q1, q2 = split(x)
    p1, p2 = split(y)
    e0 = max(e0, np.abs(q2[4:]).max(), np.abs(mul(q2, E[7]) - (x - q1)).max())
    lhs = mul(x, y)
    first = mul(q1, p1) - mul(conj(p2), q2)
    second = mul(mul(p2, q1) + mul(q2, conj(p1)), E[7])
    e1 = max(e1, np.abs(lhs - first - second).max())
    xx = mul(x, x)
    e2 = max(e2, np.abs(xx - (mul(q1, q1) - (q2 @ q2) * E[0]) - 2 * q1[0] * mul(q2, E[7])).max())
check(40, "o = q1 + q2ℓ (q1,q2 ∈ ℍ=<1,e1,e2,e3>、ℓ=e7) と分けられる", e0)
check(40, "(q1+q2ℓ)(p1+p2ℓ) = (q1p1 − p̄2q2) + (p2q1 + q2p̄1)ℓ", e1)
check(40, "O² = (q1² − |q2|²) + 2Re(q1) q2ℓ", e2)

# ---- [11] 複素平面に制限すると三体の蔵本モデル
e1 = 0.0
for _ in range(100):
    tA, tB, tC = rng.uniform(0, 2 * np.pi, 3)
    def z(t):
        o = np.zeros(8)
        o[0], o[7] = np.cos(t), np.sin(t)
        return o
    oA = z(tA)
    V = conj(mul(z(tB), z(tC)))
    G = V - (oA @ V) * oA
    tangent = np.zeros(8)
    tangent[0], tangent[7] = -np.sin(tA), np.cos(tA)
    e1 = max(e1, abs(G @ tangent + np.sin(tA + tB + tC)))
check(11, "一つの複素平面では dθ_A/dt = −sin(θ_A+θ_B+θ_C)", e1)

# ---- [9][19] 勾配流で S は増える / 距離の重みつきの流れで S_f は増える
def flow_test(weighted, Nn=6, steps=400, eps=0.01):
    os2 = np.array([rnd(True) for _ in range(Nn)])
    hist = []
    for _ in range(steps):
        Q = np.einsum('ijk,bi,cj->bck', M, os2, os2)        # o_B o_C
        c = np.einsum('ak,bck->abc', os2, np.einsum('bck,k->bck', Q, np.array([1] + [-1] * 7)))
        # c_ABC = Re(o_A o_B o_C) = <o_A, conj(o_B o_C)>
        Qc = Q * np.array([1] + [-1] * 7)
        if weighted:
            hist.append(np.sum((1 + c) ** 2 / 4))             # S_f, f = (1+c)²/4, f' = (1+c)/2
            w = (1 + c) / 2
        else:
            hist.append(np.sum(c))
            w = np.ones_like(c)
        V = np.einsum('abc,bck->ak', w, Qc)
        G = V - np.einsum('abc,abc->a', w, c)[:, None] * os2
        os2 = os2 + eps * G
        os2 /= np.linalg.norm(os2, axis=1, keepdims=True)
    d = np.diff(hist)
    return hist[0], hist[-1], d.min()


for weighted, tag, name in [(False, 9, "重み一定の勾配流で S = ΣRe(ooo) は単調に増える"),
                            (True, 17, "重み w=(1+c)/2 の流れで S_f = Σ(1+c)²/4 は単調に増える")]:
    s0, s1, dmin = flow_test(weighted)
    ok = dmin > -1e-9
    results.append(ok)
    print(f"[{tag}] {'成立' if ok else '不成立'}  {name}   {s0:.4f} → {s1:.4f}、一歩の増分の最小 = {dmin:.2e}")

# ---- [24][56] 単位系の数値
alpha = 1 / 137.035999
print(f"[24] ħ = 9 q0²/α = {9/alpha:.1f} q0²   [56] √ħ = 3 q0/√α = {3/np.sqrt(alpha):.2f} q0")

print(f"\n合計 {len(results)} 項目、成立 {sum(results)}、不成立 {len(results)-sum(results)}")
