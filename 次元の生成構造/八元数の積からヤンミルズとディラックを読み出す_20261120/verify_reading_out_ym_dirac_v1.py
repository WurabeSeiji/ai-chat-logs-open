"""検算プログラム: 八元数の積から、ヤン＝ミルズとディラックを読み出す (v1)
Python 3 / NumPy。乱数の種は固定。各項目の [§n] は論文の節番号。
実行: python3 verify_reading_out_ym_dirac_v1.py > results/verify_reading_out_ym_dirac_v1.txt
"""
import itertools
import numpy as np

rng = np.random.default_rng(20261120)

# ---------------------------------------------------------------- 八元数
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
E = np.eye(8)
F7 = f[1:, 1:, 1:]


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
    return np.einsum('ijk,i,j->k', F7, b, c)


def phi(b1, b2, b3):
    return np.einsum('ijk,i,j,k->', F7, b1, b2, b3)


def power(x, m):
    y = E[0].copy()
    for _ in range(m):
        y = mul(y, x)
    return y


def oexp(v):
    t = np.linalg.norm(v)
    return np.concatenate([[np.cos(t)], np.sin(t) * v / t])


def unit_im():
    n = rng.normal(size=7)
    return n / np.linalg.norm(n)


def expm(A, terms=40):
    k = max(0, int(np.ceil(np.log2(max(np.abs(A).sum(axis=1).max(), 1e-16)))) + 1)
    B = A / 2 ** k
    X = np.eye(A.shape[0], dtype=A.dtype)
    T = np.eye(A.shape[0], dtype=A.dtype)
    for n in range(1, terms):
        T = T @ B / n
        X = X + T
    for _ in range(k):
        X = X @ X
    return X


results = []


def check(sec, name, err, tol=1e-9):
    ok = bool(err < tol)
    results.append(ok)
    print(f"[§{sec}] {'成立' if ok else '不成立'}  {name}   最大誤差 = {err:.2e}")


def note(sec, name, ok, detail):
    results.append(bool(ok))
    print(f"[§{sec}] {'成立' if ok else '不成立'}  {name}   {detail}")


print("=" * 20, "事前の検討 (§2)", "=" * 20)
# ---- §2.2 全員と一度に作用するとパワーが保存される
N = 4
z = rng.normal(size=N) + 1j * rng.normal(size=N)
K = np.imag(np.conj(z)[:, None] * z[None, :])          # K_ef = Im(conj(z_e) z_f)、反対称
P0 = np.sum(np.abs(z) ** 2)
U = expm(0.3 * K)
z_all = U @ z
z_add = z.copy()                                         # 対ごとに別々に計算して足す
for e in range(N):
    for g in range(e + 1, N):
        Kp = np.zeros((N, N))
        Kp[e, g], Kp[g, e] = K[e, g], K[g, e]
        z_add = z_add + (expm(0.3 * Kp) @ z - z)
def seq(order):
    w = z.copy()
    for (e, g) in order:
        Kp = np.zeros((N, N))
        Kp[e, g], Kp[g, e] = K[e, g], K[g, e]
        w = expm(0.3 * Kp) @ w
    return w
pairs = [(e, g) for e in range(N) for g in range(e + 1, N)]
z_s1, z_s2 = seq(pairs), seq(pairs[::-1])
check("2.2", "一つの N×N の回転 exp(K) でパワーが不変", abs(np.sum(np.abs(z_all) ** 2) - P0))
check("2.2", "回転の各行で |U_ij|² の和が 1", np.abs((np.abs(U) ** 2).sum(axis=1) - 1).max())
note("2.2", "対ごとに別々に計算して足すとパワーが変わる", abs(np.sum(np.abs(z_add) ** 2) - P0) > 1e-3,
     f"元 {P0:.6f} → {np.sum(np.abs(z_add)**2):.6f}")
note("2.2", "対ごとに順番に施すと、順番で結果が変わる", np.abs(z_s1 - z_s2).max() > 1e-3,
     f"パワーは {np.sum(np.abs(z_s1)**2):.6f} で不変、順番を逆にした差 {np.abs(z_s1-z_s2).max():.4f}")

# ---- §2.4 倍音を下に積む
psi = np.linspace(0, 2 * np.pi, 36000, endpoint=False)


def n_top(y, frac=0.999):
    lm = np.where((y > np.roll(y, 1)) & (y >= np.roll(y, -1)))[0]
    return int(np.sum(y[lm] > frac * y.max()))


Mh = 12
single = np.cos(Mh * psi)
up = sum(np.cos(k * Mh * psi) for k in range(1, 5))
down = sum(np.cos(m * psi) for m in range(1, Mh + 1))
note("2.4", "巻き数 12 の波だけ → 山 12 個、上に積んでも 12 個、下に 1 次まで積むと 1 個",
     (n_top(single), n_top(up), n_top(down)) == (12, 12, 1),
     f"山の数 = {n_top(single)}, {n_top(up)}, {n_top(down)}")
psi2 = np.linspace(0, 2 * np.pi, 720000, endpoint=False)
rows = []
for Mx in (6, 12, 24, 48):
    y = sum(np.cos(m * (psi2 - np.radians(137.0))) for m in range(1, Mx + 1))
    pos = np.degrees(psi2[np.argmax(y)])
    yy = y / y.max()
    width = np.degrees(np.sum(yy > 0.5) * (psi2[1] - psi2[0]))
    rows.append((Mx, pos, width))
note("2.4", "中心は次数によらずベース波の中心位相、幅は最高次数にほぼ反比例",
     all(abs(p - 137.0) < 1e-2 for _, p, _ in rows) and abs(rows[0][2] / rows[-1][2] - 8) < 1.5,
     " / ".join(f"M={m}: 位置 {p:.3f}°, 半値幅 {w:.2f}°" for m, p, w in rows))
ac = {}
for name, orders in [("奇偶両方", range(1, 13)), ("奇数だけ", range(1, 13, 2)), ("偶数だけ", range(2, 13, 2))]:
    ac[name] = sum(np.cos(m * np.pi) for m in orders) / len(list(orders))
note("2.4", "自己相関の 180° の値: 奇偶両方 0、奇数だけ −1、偶数だけ ＋1",
     abs(ac["奇偶両方"]) < 1e-12 and abs(ac["奇数だけ"] + 1) < 1e-12 and abs(ac["偶数だけ"] - 1) < 1e-12,
     str({k: round(float(v), 3) for k, v in ac.items()}))
e1 = 0.0
for _ in range(300):
    a, b = rnd(), rnd()
    cosang = a @ b / np.linalg.norm(a) / np.linalg.norm(b)
    ab = mul(conj(a), b)
    e1 = max(e1, abs(ab[0] - np.linalg.norm(a) * np.linalg.norm(b) * cosang),
             abs(np.linalg.norm(ab[1:]) - np.linalg.norm(a) * np.linalg.norm(b) * np.sqrt(1 - cosang ** 2)))
check("2.5", "Re(āb) = |a||b|cos、|Im(āb)| = |a||b|sin (二つが張る面積)", e1)

print("=" * 20, "思考実験 (§3〜§18)", "=" * 20)
# ---- §4 そのままの二乗
Fm = [rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)) for _ in range(7)]
sq = [np.zeros((3, 3), complex) for _ in range(8)]
for i in range(1, 8):
    for j in range(1, 8):
        for k in range(8):
            if M[i, j, k] != 0:
                sq[k] += M[i, j, k] * Fm[i - 1] @ Fm[j - 1]
err = np.abs(sq[0] + sum(Fi @ Fi for Fi in Fm)).max()
for k in range(1, 8):
    t = sum(0.5 * f[i, j, k] * (Fm[i - 1] @ Fm[j - 1] - Fm[j - 1] @ Fm[i - 1]) for i in range(1, 8) for j in range(1, 8))
    err = max(err, np.abs(sq[k] - t).max())
check(4, "𝔉² = −ΣF_iF_i + ½ f_ijk[F_i,F_j]e_k", err)
e1 = 0.0
for _ in range(200):
    x, y = rnd(), rnd()
    e1 = max(e1, np.abs(mul(mul(x, y), x) - mul(x, mul(y, x))).max(), np.abs(mul(mul(x, x), y) - mul(x, mul(x, y))).max())
check(4, "交代性: 同じ二つの元の積は括弧によらない", e1)

# ---- §5 三つの積
e1 = e2 = e3 = e4 = e5 = e6 = 0.0
for _ in range(400):
    x, y, zz = rnd(), rnd(), rnd()
    r1 = mul(mul(x, y), zz)[0]
    e1 = max(e1, abs(r1 - mul(x, mul(y, zz))[0]))
    e2 = max(e2, abs(r1 - mul(mul(y, zz), x)[0]))
    a1, a2, a3, b1, b2, b3 = x[0], y[0], zz[0], x[1:], y[1:], zz[1:]
    e3 = max(e3, abs(r1 - (a1 * a2 * a3 - a1 * (b2 @ b3) - a2 * (b1 @ b3) - a3 * (b1 @ b2) - phi(b1, b2, b3))))
    xy = mul(x, y)
    e4 = max(e4, abs(xy[0] - (a1 * a2 - b1 @ b2)) + np.abs(xy[1:] - (a1 * b2 + a2 * b1 + cross(b1, b2))).max())
    e5 = max(e5, abs(np.linalg.norm(xy) - np.linalg.norm(x) * np.linalg.norm(y)))
    w = [rnd() for _ in range(4)]
    # ((ab)c)d と (a(bc))d の実部の差 = 結合子 (a,b,c) と d の積の実部
    e6 = max(e6, abs(mul(mul(mul(w[0], w[1]), w[2]), w[3])[0] - mul(mul(w[0], mul(w[1], w[2])), w[3])[0]))
check(5, "Re((xy)z) = Re(x(yz))", e1)
check(5, "Re(xyz) = Re(yzx)  巡回的", e2)
check(5, "Re(o1o2o3) = a1a2a3 − a1(b2·b3) − a2(b1·b3) − a3(b1·b2) − φ(b1,b2,b3)", e3)
check(5, "xy = (a_x a_y − b_x·b_y) + (a_x b_y + a_y b_x + b_x×b_y)", e4)
check(5, "|xy| = |x||y|  ノルムは積で保たれる", e5)
note(5, "四つの積の実部は括弧の付け方による", e6 > 1e-3, f"差の最大 = {e6:.3f}")
e1 = e2 = 0.0
for _ in range(200):
    o2, o3 = rnd(True), rnd(True)
    o1n = conj(mul(o2, o3))
    e1 = max(e1, abs(mul(mul(o1n, o2), o3)[0] - 1))
    e2 = max(e2, np.abs(conj(mul(o3, o1n)) - o2).max())
check(5, "o1' = conj(o2o3) とすると Re(o1'o2o3) = 1、続く o2' は o2 のまま", max(e1, e2))

# ---- §6 ヤン＝ミルズの形
a, b = rng.normal(size=7), rng.normal(size=7)
ratios = []
for eps in (0.1, 0.05, 0.025, 0.0125):
    S = mul(mul(oexp(eps * a), oexp(eps * b)), oexp(eps * (-a - b)))[0]
    ratios.append((1 - S) / (0.5 * eps ** 4 * cross(a, b) @ cross(a, b)))
note(6, "一次で閉じた三辺: 1 − Re(ooo) → (ε⁴/2)|a×b|²", abs(ratios[-1] - 1) < 2e-3 and abs(ratios[-1] - 1) < abs(ratios[0] - 1),
     "比 = " + ", ".join(f"{r:.5f}" for r in ratios) + "  (ε = 0.1, 0.05, 0.025, 0.0125)")
c = rng.normal(size=7)
ratios = []
for eps in (0.05, 0.025, 0.0125, 0.00625):
    S = mul(mul(oexp(eps * a), oexp(eps * b)), oexp(eps * c))[0]
    Fv = eps * (a + b + c) + eps ** 2 * (cross(a, b) + cross(a, c) + cross(b, c))
    ratios.append((1 - S) / (0.5 * Fv @ Fv))
note(6, "一般の三辺: 1 − Re(ooo) → ½|𝔉|²、𝔉 = ε(a+b+c) + ε²(a×b+a×c+b×c)", abs(ratios[-1] - 1) < 2e-3,
     "比 = " + ", ".join(f"{r:.5f}" for r in ratios))
S = mul(mul(oexp(0.1 * a), oexp(0.07 * a)), oexp(-0.17 * a))[0]
check(6, "同じ複素平面の三辺では、一次で閉じれば閉じ残りは 0", abs(1 - S))

# ---- §7 ディラックの形
L = [np.einsum('ijk,i->kj', M, E[i]) for i in range(8)]
err = max(np.abs(L[i] @ L[j] + L[j] @ L[i] + 2 * (i == j) * np.eye(8)).max() for i in range(1, 8) for j in range(1, 8))
check(7, "{L_ei, L_ej} = −2δ_ij  左からの掛け算はクリフォード代数", err)
note(7, "L_x L_y ≠ L_xy  左からの掛け算は群の表現にならない",
     np.abs(L[1] @ L[2] - np.einsum('ijk,i->kj', M, mul(E[1], E[2]))).max() > 1, "L_e1 L_e2 と L_(e1e2) の差の最大 = 2")
e1 = e2 = e3 = 0.0
for _ in range(300):
    x, p = rnd(), rnd()
    e1 = max(e1, np.abs(mul(p, mul(conj(p), x)) - (p @ p) * x).max())
    oA, oB, oC = rnd(), rnd(), rnd()
    LB = np.einsum('ijk,i->kj', M, oB)
    e2 = max(e2, abs(mul(mul(oA, oB), oC)[0] - conj(oA) @ (LB @ oC)))
    u = rnd(True)
    Lu = np.einsum('ijk,i->kj', M, u)
    e3 = max(e3, np.abs(Lu @ Lu.T - np.eye(8)).max())
check(7, "𝔇(𝔇̄ x) = |p|² x  二乗するとラプラシアン", e1)
check(7, "Re(o_A o_B o_C) = ψ̄ᵀ L(o_B) ψ  (ψ̄ = conj(o_A), ψ = o_C)", e2)
check(7, "|o_B| = 1 のとき L(o_B) は回転 (直交行列)", e3)
# 一本の式 𝔇𝔉 = 𝔍 の代数 (論文では §3)
chi = np.zeros((8, 8, 8, 8))
for r, m, n in itertools.product(range(8), repeat=3):
    chi[r, m, n] = mul(E[r], mul(conj(E[m]), E[n])) - ((r == m) * E[n] - (r == n) * E[m] + (m == n) * E[r])
err = max(np.abs(chi + np.transpose(chi, p)).max() for p in [(1, 0, 2, 3), (0, 2, 1, 3), (0, 1, 3, 2)])
check(3, "e_ρ(ē_μ e_ν) = δ_ρμ e_ν − δ_ρν e_μ + δ_μν e_ρ + χ_ρμνσ e_σ、χ は完全反対称", err)
T = rng.normal(size=(8, 8, 8))
T = T - np.transpose(T, (0, 2, 1))                               # T_ρμν = D_ρ F_μν (μν で反対称)
cyc = (T + np.transpose(T, (1, 2, 0)) + np.transpose(T, (2, 0, 1))) / 3
T = T - cyc                                                       # ビアンキ恒等式 D_[ρ F_μν] = 0
lhs = np.zeros(8)
for r, m, n in itertools.product(range(8), repeat=3):
    lhs += 0.5 * T[r, m, n] * mul(E[r], mul(conj(E[m]), E[n]))
rhs = np.einsum('mmn->n', T)
check(3, "ビアンキ恒等式のもとで 𝔇𝔉 の e_ν 成分 = D_μF_μν", np.abs(lhs - rhs).max())
err = max(abs(mul(conj(E[m]), E[n])[0] - (m == n)) for m in range(8) for n in range(8))
check(3, "Re(ē_μ e_ν) = δ_μν  (𝔉 = ½Σ ē_μ e_ν F_μν は虚部だけ、7 成分)", err)
pairs8 = [(m, n) for m in range(8) for n in range(m + 1, 8)]
Top = np.array([[chi[m, n, r, s] for (r, s) in pairs8] for (m, n) in pairs8])
ev = np.round(np.linalg.eigvalsh(Top), 6)
vals, cnt = np.unique(ev, return_counts=True)
note(3, "28 成分は 7 ＋ 21 に分かれる (固有値 3 が 7 個、−1 が 21 個)",
     dict(zip(vals.tolist(), cnt.tolist())) == {-1.0: 21, 3.0: 7}, str(dict(zip(vals.tolist(), cnt.tolist()))))

# ---- §9 無名な和
Nn = 7
os_ = [rnd(True) for _ in range(Nn)]
S = sum(mul(mul(x, y), w)[0] for x in os_ for y in os_ for w in os_)
P = sum(os_)
Rr, th = np.linalg.norm(P), np.arccos(P[0] / np.linalg.norm(P))
check(9, "Σ_{A,B,C} Re(o_A o_B o_C) = Re(P³) = a³ − 3a|b|² = R³cos3θ",
      max(abs(S - mul(mul(P, P), P)[0]), abs(S - (P[0] ** 3 - 3 * P[0] * (P[1:] @ P[1:]))), abs(S - Rr ** 3 * np.cos(3 * th))))
S2 = sum(mul(mul(x, y), w)[0] ** 2 for x in os_ for y in os_ for w in os_)
note(9, "二乗して足す和は総和 P の関数に潰れない", abs(S2 - S ** 2) > 1, f"Σ[Re(ooo)]² = {S2:.3f}、[Re(P³)]² = {S**2:.3f}")

# ---- §10 三角法
e1 = e2 = 0.0
for _ in range(300):
    x, y, w = rnd(True), rnd(True), rnd(True)
    e1 = max(e1, abs(np.sum((x - y) ** 2) - (2 - 2 * mul(x, conj(y))[0])))
    e2 = max(e2, abs(np.sum((mul(x, y) - conj(w)) ** 2) - (2 - 2 * mul(mul(x, y), w)[0])))
check(10, "|x − y|² = 2 − 2Re(x ȳ)", e1)
check(10, "|o_A o_B − ō_C|² = 2 − 2Re(o_A o_B o_C)", e2)

# ---- §11 一つの複素平面では三体の蔵本型
e1 = 0.0
for _ in range(100):
    tA, tB, tC = rng.uniform(0, 2 * np.pi, 3)
    zc = lambda t: np.concatenate([[np.cos(t)], np.zeros(6), [np.sin(t)]])
    V = conj(mul(zc(tB), zc(tC)))
    G = V - (zc(tA) @ V) * zc(tA)
    tang = np.concatenate([[-np.sin(tA)], np.zeros(6), [np.cos(tA)]])
    e1 = max(e1, abs(G @ tang + np.sin(tA + tB + tC)))
check(11, "一つの複素平面では、接線方向の成分 = −sin(θ_A+θ_B+θ_C)", e1)

# ---- §11 参考: §5 の連続の流れ (§17 で退けられる更新) を小さな系で一例だけ走らせる
def flow_once(weighted, Nn=6, steps=400, eps=0.01):
    os2 = np.array([rnd(True) for _ in range(Nn)])
    sign = np.array([1] + [-1] * 7)
    for _ in range(steps):
        Qc = np.einsum('ijk,bi,cj->bck', M, os2, os2) * sign          # conj(o_B o_C)
        c = np.einsum('ak,bck->abc', os2, Qc)                          # c_ABC = Re(o_A o_B o_C)
        w = (1 + c) / 2 if weighted else np.ones_like(c)
        G = np.einsum('abc,bck->ak', w, Qc) - np.einsum('abc,abc->a', w, c)[:, None] * os2
        os2 = os2 + eps * G
        os2 /= np.linalg.norm(os2, axis=1, keepdims=True)
    Qc = np.einsum('ijk,bi,cj->bck', M, os2, os2) * sign
    return np.einsum('ak,bck->abc', os2, Qc).min()


print(f"[§11] 参考  連続の流れを N = 6 で一例: 重み一定で c の最小 = {flow_once(False):.4f}、"
      f"距離の重み (1+c)/2 で c の最小 = {flow_once(True):.4f}  (1 なら全員が同じ値に揃った)")

# ---- §12 倍音の和、干渉項、視野角
e1 = 0.0
for _ in range(200):
    Nh = int(rng.integers(2, 40))
    t = rng.uniform(0.01, np.pi - 0.01)
    s = sum(np.cos(k * t) for k in range(1, Nh + 1))
    e1 = max(e1, abs(s - np.cos((Nh + 1) * t / 2) * np.sin(Nh * t / 2) / np.sin(t / 2)),
             abs(s - (np.sin((Nh + 0.5) * t) / (2 * np.sin(t / 2)) - 0.5)))
check(12, "Σ_{n=1}^N cos nθ = cos((N+1)θ/2)sin(Nθ/2)/sin(θ/2) = sin((N+½)θ)/(2sin(θ/2)) − ½", e1)
Nh = 400
rows = []
for D in (0.3, 0.6, 1.2, 2.4, np.pi):
    per = 2 * np.pi / (Nh + 0.5)                       # 速い振動の周期。その整数倍の幅で平均する
    ds = D + np.linspace(-per, per, 2000, endpoint=False)   # 幅が広いと 1/sin² 自体の変化が混ざる
    s2 = np.mean([sum(np.cos(k * d) for k in range(1, Nh + 1)) ** 2 for d in ds])
    rows.append((D, s2, 1 / (8 * np.sin(D / 2) ** 2) + 0.25))
note(12, "幅の外: 干渉項の二乗の平均 ≈ 1/(8sin²(Δ/2)) ＋ ¼ (小さい Δ で 1/(2Δ²))",
     all(abs(a_ / b_ - 1) < 0.02 for _, a_, b_ in rows),
     " / ".join(f"Δ={D:.2f}: {a_:.4f} と {b_:.4f}" for D, a_, b_ in rows))
note(12, "幅の内: Δ ≪ 幅 では和は次数 N で一定", abs(sum(np.cos(k * 1e-5) for k in range(1, Nh + 1)) / Nh - 1) < 1e-4,
     f"N = {Nh}、Δ = 1e-5 で 和/N = {sum(np.cos(k * 1e-5) for k in range(1, Nh + 1)) / Nh:.6f}")
e1 = 0.0
for _ in range(200):
    d = rng.uniform(0.3, np.pi - 0.3)
    rho = rng.uniform(0.02, min(d, np.pi - d) * 0.9)
    cen = np.array([np.sin(d), 0, np.cos(d)])
    u1 = np.array([np.cos(d), 0, -np.sin(d)])
    u2 = np.array([0, 1.0, 0])
    al = np.linspace(0, 2 * np.pi, 20001)
    pts = np.cos(rho) * cen[None, :] + np.sin(rho) * (np.cos(al)[:, None] * u1 + np.sin(al)[:, None] * u2)
    az = np.arctan2(pts[:, 1], pts[:, 0])
    e1 = max(e1, abs(np.sin(az.max()) - np.sin(rho) / np.sin(d)))
check(12, "球面の上で、角半径 ρ の相手を角距離 d から見込む角 ψ: sinψ = sinρ/sin d", e1, tol=1e-6)
alpha = 1 / 137.035999
print(f"[§12] 参考  c = G = 1、q0 = e/3 の単位系: ħ = 9q0²/α = {9/alpha:.1f} q0²、√ħ = 3q0/√α = {3/np.sqrt(alpha):.2f} q0")

# ---- §13 中心位相、O²、R
e1 = e2 = e3 = e4 = 0.0
for _ in range(200):
    n1, n2 = unit_im(), unit_im()
    A, t1, t2 = rng.uniform(0.5, 2), rng.uniform(0, np.pi), rng.uniform(0, np.pi)
    O = np.concatenate([[A * np.cos(t1)], A * np.sin(t1) * n1])
    for m in (2, 3, 5):
        e1 = max(e1, np.abs(power(O, m) - A ** m * np.concatenate([[np.cos(m * t1)], np.sin(m * t1) * n1])).max())
    O2 = mul(O, O)
    e2 = max(e2, abs(O2[0] - (O[0] ** 2 - O[1:] @ O[1:])) + np.abs(O2[1:] - 2 * O[0] * O[1:]).max())
    ps = np.arccos(n1 @ n2)
    for m in (1, 2, 4):
        Am = np.concatenate([[np.cos(m * t1)], np.sin(m * t1) * n1])
        Bm = np.concatenate([[np.cos(m * t2)], np.sin(m * t2) * n2])
        e3 = max(e3, abs(mul(Am, conj(Bm))[0] - (np.cos(ps / 2) ** 2 * np.cos(m * (t1 - t2)) + np.sin(ps / 2) ** 2 * np.cos(m * (t1 + t2)))))
    Nh = int(rng.integers(2, 30))
    t = rng.uniform(0.01, np.pi - 0.01)
    re = sum(np.cos(2 * k * t) for k in range(1, Nh + 1))
    im = sum(np.sin(2 * k * t) for k in range(1, Nh + 1))
    pre = np.sin(Nh * t) / np.sin(t)
    e4 = max(e4, abs(re - pre * np.cos((Nh + 1) * t)) + abs(im - pre * np.sin((Nh + 1) * t)))
check(13, "Oᵐ = Aᵐ(cos mθ + n sin mθ)  倍音は冪", e1)
check(13, "O² = (a² − |b|²) + 2ab", e2)
check(13, "向きの違う倍音の干渉の実部 = cos²(ψ/2)cos m(θ_A−θ_B) + sin²(ψ/2)cos m(θ_A+θ_B)", e3)
check(13, "Σ_n o_kn² = a²·sin(Nθ)/sinθ·[cos((N+1)θ) + n_k sin((N+1)θ)]", e4)

# ---- §14 二次の展開、グノモン投影
Nh, t = 200, 1e-4
re = sum(np.cos(2 * k * t) for k in range(1, Nh + 1)) / Nh
im = sum(np.sin(2 * k * t) for k in range(1, Nh + 1)) / Nh
note(14, "R_k² ≈ P_k[1 + n_k(N+1)θ − (2/3)N²θ²]  (N が大きく Nθ が小さいとき)",
     abs(im / ((Nh + 1) * t) - 1) < 1e-3 and abs((1 - re) / (2 / 3 * Nh ** 2 * t ** 2) - 1) < 2e-2,
     f"一次 {im:.6e} と {(Nh+1)*t:.6e}、二次 {1-re:.6e} と {2/3*Nh**2*t**2:.6e}")
e1 = 0.0
for _ in range(300):
    x = rng.normal(size=3)
    dx = 1e-5 * rng.normal(size=3)
    def to_sphere(v):
        w = np.concatenate([v, [1.0]])
        return w / np.linalg.norm(w)
    ds2 = np.arccos(np.clip(to_sphere(x) @ to_sphere(x + dx), -1, 1)) ** 2
    r2 = x @ x
    form = dx @ dx / (1 + r2) - (x @ dx) ** 2 / (1 + r2) ** 2
    e1 = max(e1, abs(ds2 / form - 1))
check(14, "グノモン投影の計量 ds² = dx·dx/(1+r²) − (x·dx)²/(1+r²)²", e1, tol=1e-3)

# ---- §16 二つの四元数として読む
def split(o):
    q1 = np.zeros(8)
    q1[:4] = o[:4]
    return q1, -mul(o - q1, E[7])


e0 = e1 = e2 = 0.0
for _ in range(300):
    x, y = rnd(), rnd()
    q1, q2 = split(x)
    p1, p2 = split(y)
    e0 = max(e0, np.abs(q2[4:]).max(), np.abs(mul(q2, E[7]) - (x - q1)).max())
    e1 = max(e1, np.abs(mul(x, y) - (mul(q1, p1) - mul(conj(p2), q2)) - mul(mul(p2, q1) + mul(q2, conj(p1)), E[7])).max())
    e2 = max(e2, np.abs(mul(x, x) - (mul(q1, q1) - (q2 @ q2) * E[0]) - 2 * q1[0] * mul(q2, E[7])).max())
check(16, "o = q1 + q2ℓ  (q1, q2 は四元数 <1,e1,e2,e3>、ℓ = e7)", e0)
check(16, "(q1+q2ℓ)(p1+p2ℓ) = (q1p1 − p̄2q2) + (p2q1 + q2p̄1)ℓ", e1)
check(16, "O² = (q1² − |q2|²) + 2Re(q1) q2ℓ", e2)
e1 = 0.0
for _ in range(200):
    s, x, y = rng.normal(), rnd(), rnd()
    sc = s * E[0]
    e1 = max(e1, np.abs(mul(sc, x) - mul(x, sc)).max(), np.abs(mul(mul(sc, x), y) - mul(sc, mul(x, y))).max(),
             np.abs(mul(mul(x, sc), y) - mul(x, mul(sc, y))).max())
check(16, "実数 (1＋0) は、どの八元数とも交換し、括弧の付け方によらない", e1)
note(16, "虚数単位どうしは交換しない (交換するのは実数だけ)", np.abs(mul(E[1], E[2]) - mul(E[2], E[1])).max() > 1, "e1e2 − e2e1 = 2e3")

# ---- §17 有限位数
n = 12
ks = rng.integers(0, n, size=5)
Ud = np.exp(2j * np.pi * ks / n)
psi0 = rng.normal(size=5) + 1j * rng.normal(size=5)
fobs = lambda v: np.real(v[0] * np.conj(v[1])) + np.abs(v[2] + v[3]) ** 2 + np.real(v[4] ** 2)
orb = lambda v: sum(fobs(Ud ** t * v) for t in range(n))
single_change = abs(fobs(Ud * psi0) - fobs(psi0))
check(17, "Uⁿ = I のとき、一周ぶんの和は ψ → Uψ で変わらない", abs(orb(Ud * psi0) - orb(psi0)))
note(17, "一回ごとの値は変わる (一回だけ見ると保存されないように見える)", single_change > 1e-3, f"一回での変化 = {single_change:.4f}")
check(17, "1 の n 乗根の総和は 0、m 乗の総和も m が n の倍数でなければ 0",
      max(abs(np.sum(np.exp(2j * np.pi * np.arange(n) * m / n))) for m in range(1, n)))

print(f"\n合計 {len(results)} 項目、成立 {sum(results)}、不成立 {len(results)-sum(results)}")
