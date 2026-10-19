#!/usr/bin/env python3
"""
八元数閉塞と関係の分解（検討 2026-10-19）

目的
  1. 八元数 O の代数（ケーリー＝ディクソン構成）を実装し、基本性質を数値で確認する。
  2. 「本当の二乗」の閉塞 Σo = 0, Σo² = 0 を保つ連続対称性の次元を数える
     （無限小写像の核の次元）。比較：積を保つ微分（G2, 14 次元）、ℓ を固定する部分（su(3), 8 次元）。
  3. 閉塞の厳密解を C, H, O で生成し、解多様体の次元・対称性の軌道の次元・形の空間の次元を
     N = 3..10 について数値ヤコビアンの階数から求める。
  4. 関係 ō o′ の分解：実部 = ⟨q,q′⟩ + ⟨p,p′⟩（時空側と内部側の二つの内積）、ℓ 部 = 混合項。
  5. 閉塞が内部半分のスケール κ（内と外の換算比の居場所）を決めること。
  6. 三角形の枠ホロノミー：比の一歩 s_ef = o_f ō_e を三角形で一周したとき、
     H に乗る枠の回転角。H では恒等、O では非結合性により非自明。

再現性
  乱数シード 20261019。出力は results/ 以下（JSON, md）、図は figures/ 以下（PNG, SVG）。
"""
import json, os, sys
import numpy as np

SEED = 20261019
rng = np.random.default_rng(SEED)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------------------
# 四元数（[w, x, y, z]）
# ---------------------------------------------------------------------------
def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])
def qconj(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])

# ---------------------------------------------------------------------------
# 八元数 o = (a, b), a, b ∈ H。ケーリー＝ディクソン：(a,b)(c,d) = (ac − d̄b, da + bc̄)
# 基底：e0=1, e1..e3 = i,j,k, e4 = ℓ, e5..e7 = iℓ, jℓ, kℓ
# ---------------------------------------------------------------------------
def omul(o1, o2):
    a, b = o1[:4], o1[4:]; c, d = o2[:4], o2[4:]
    return np.concatenate([qmul(a, c) - qmul(qconj(d), b), qmul(d, a) + qmul(b, qconj(c))])
def oconj(o):
    return np.concatenate([qconj(o[:4]), -o[4:]])
def onorm2(o):
    return float(np.dot(o, o))
def oinv(o):
    return oconj(o) / onorm2(o)
def ojordan(x, y):
    return 0.5 * (omul(x, y) + omul(y, x))
def oexp(A):
    """純虚八元数 A の指数 exp(A) = cos|A| + (A/|A|) sin|A|"""
    t = np.linalg.norm(A)
    if t < 1e-300: return np.eye(8)[0]
    return np.cos(t) * np.eye(8)[0] + np.sin(t) * A / t

E = np.eye(8)

results = {"seed": SEED}

# ---------------------------------------------------------------------------
# 1. 代数の基本性質
# ---------------------------------------------------------------------------
def alg_checks(trials=2000):
    mx_norm = mx_alt = 0.0; mn_assoc = np.inf; mx_assoc = 0.0
    for _ in range(trials):
        x, y, z = rng.normal(size=(3, 8))
        mx_norm = max(mx_norm, abs(onorm2(omul(x, y)) - onorm2(x) * onorm2(y)) / (onorm2(x) * onorm2(y)))
        mx_alt = max(mx_alt, np.linalg.norm(omul(omul(x, x), y) - omul(x, omul(x, y))),
                     np.linalg.norm(omul(omul(y, x), x) - omul(y, omul(x, x))))
        assoc = np.linalg.norm(omul(omul(x, y), z) - omul(x, omul(y, z)))
        mn_assoc = min(mn_assoc, assoc); mx_assoc = max(mx_assoc, assoc)
    # 四元数部分代数 (b = 0) では結合的
    xs = [np.concatenate([rng.normal(size=4), np.zeros(4)]) for _ in range(3)]
    assoc_H = np.linalg.norm(omul(omul(xs[0], xs[1]), xs[2]) - omul(xs[0], omul(xs[1], xs[2])))
    # 基底積の表
    table = [["" for _ in range(8)] for _ in range(8)]
    names = ["1", "e1", "e2", "e3", "e4", "e5", "e6", "e7"]
    for i in range(8):
        for j in range(8):
            r = omul(E[i], E[j]); k = int(np.argmax(np.abs(r))); s = "-" if r[k] < 0 else ""
            table[i][j] = s + names[k]
    return {"max_rel_norm_mult_err": mx_norm, "max_alternativity_err": mx_alt,
            "associator_norm_min_max_random": [mn_assoc, mx_assoc], "associator_in_quaternion_subalgebra": assoc_H,
            "basis_table": table}
results["algebra"] = alg_checks()
print("algebra:", {k: v for k, v in results["algebra"].items() if k != "basis_table"})

# ---------------------------------------------------------------------------
# 2. 無限小対称性の次元（線形方程式の核の次元）
# ---------------------------------------------------------------------------
def nullspace_dim(M, tol=1e-9):
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s < tol * max(1.0, s[0]))) + (M.shape[1] - len(s))

def derivation_dims():
    # D ∈ R^{8x8}（列ベクトル作用）。条件 D(e_i e_j) = D(e_i) e_j + e_i D(e_j)
    rows = []
    for i in range(8):
        for j in range(8):
            # 各成分 c について、D の 64 成分に対する線形式
            prod = omul(E[i], E[j])
            for c in range(8):
                row = np.zeros(64)
                # D(prod)_c = Σ_m D[c,m] prod_m
                for m in range(8): row[c*8 + m] += prod[m]
                # −(D(e_i) e_j)_c = −Σ_m D[m,i] (e_m e_j)_c
                for m in range(8): row[m*8 + i] -= omul(E[m], E[j])[c]
                # −(e_i D(e_j))_c = −Σ_m D[m,j] (e_i e_m)_c
                for m in range(8): row[m*8 + j] -= omul(E[i], E[m])[c]
                rows.append(row)
    M = np.array(rows)
    dim_g2 = nullspace_dim(M)
    # ℓ = e4 を固定：D e4 = 0 → D[:,4] = 0 を追加
    extra = []
    for c in range(8):
        row = np.zeros(64); row[c*8 + 4] = 1.0; extra.append(row)
    dim_stab = nullspace_dim(np.vstack([M, np.array(extra)]))
    # 本当の二乗を保つ無限小写像：D(x∘y) = D(x)∘y + x∘D(y)（ジョルダン積、対称化）
    rows = []
    for i in range(8):
        for j in range(i, 8):
            prod = ojordan(E[i], E[j])
            for c in range(8):
                row = np.zeros(64)
                for m in range(8): row[c*8 + m] += prod[m]
                for m in range(8): row[m*8 + i] -= ojordan(E[m], E[j])[c]
                for m in range(8): row[m*8 + j] -= ojordan(E[i], E[m])[c]
                rows.append(row)
    dim_sq = nullspace_dim(np.array(rows))
    # 四元数で同じ計算（比較）：4x4
    E4 = np.eye(4); rows = []
    for i in range(4):
        for j in range(4):
            prod = qmul(E4[i], E4[j])
            for c in range(4):
                row = np.zeros(16)
                for m in range(4): row[c*4 + m] += prod[m]
                for m in range(4): row[m*4 + i] -= qmul(E4[m], E4[j])[c]
                for m in range(4): row[m*4 + j] -= qmul(E4[i], E4[m])[c]
                rows.append(row)
    dim_aut_H = nullspace_dim(np.array(rows))
    rows = []
    for i in range(4):
        for j in range(i, 4):
            prod = 0.5 * (qmul(E4[i], E4[j]) + qmul(E4[j], E4[i]))
            for c in range(4):
                row = np.zeros(16)
                for m in range(4): row[c*4 + m] += prod[m]
                for m in range(4): row[m*4 + i] -= 0.5 * (qmul(E4[m], E4[j]) + qmul(E4[j], E4[m]))[c]
                for m in range(4): row[m*4 + j] -= 0.5 * (qmul(E4[i], E4[m]) + qmul(E4[m], E4[i]))[c]
                rows.append(row)
    dim_sq_H = nullspace_dim(np.array(rows))
    return {"O_derivations_dim(G2)": dim_g2, "O_derivations_fixing_l_dim(su3)": dim_stab,
            "O_square_preserving_dim(so7)": dim_sq,
            "H_derivations_dim(so3)": dim_aut_H, "H_square_preserving_dim(so3)": dim_sq_H}
results["symmetry_dims"] = derivation_dims()
print("symmetry dims:", results["symmetry_dims"])

# ---------------------------------------------------------------------------
# 3. 閉塞の厳密解と次元の数え上げ（C, H, O を同じ枠で）
#    要素 x = x0 + X（X は虚部、dim d−1）。Σx = 0, Σx² = 0 ⇔
#    Σx0 = 0, ΣX_i = 0 (各虚成分), Σx0² = Σ|X|², Σ x0·X_i = 0 (各 i)
# ---------------------------------------------------------------------------
def gen_closure_solution(d, N):
    """d = 2, 4, 8。戻り値 (N, d) の実配列。厳密解。"""
    one = np.ones(N) / np.sqrt(N)
    x0 = rng.normal(size=N); x0 -= one * (one @ x0)              # 和零
    P = np.eye(N) - np.outer(one, one) - np.outer(x0, x0) / (x0 @ x0)  # 和零 ∩ x0⊥ への射影
    X = np.array([P @ rng.normal(size=N) for _ in range(d - 1)])  # (d−1, N)
    s = np.sqrt((x0 @ x0) / np.sum(X * X))
    X *= s
    return np.column_stack([x0, X.T])

def jordan_d(d, x, y):
    if d == 2:
        return np.array([x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0]])
    if d == 4:
        return 0.5 * (qmul(x, y) + qmul(y, x))
    return ojordan(x, y)

def square_d(d, x):
    return jordan_d(d, x, x)

def closure_residual(d, S):
    lin = np.sum(S, axis=0)
    sq = np.sum([square_d(d, x) for x in S], axis=0)
    return float(np.linalg.norm(lin)), float(np.linalg.norm(sq))

def constraint_jacobian(d, S):
    """F(S) = (Σx, Σx²) ∈ R^{2d}。J は (2d, N d)。"""
    N = S.shape[0]; J = np.zeros((2 * d, N * d)); Ed = np.eye(d)
    for n in range(N):
        for k in range(d):
            J[:d, n*d + k] = Ed[k]
            J[d:, n*d + k] = 2.0 * jordan_d(d, S[n], Ed[k])
    return J

def symmetry_tangents(d, S):
    """閉塞を保つ連続対称性の無限小作用（配置の接ベクトル）。
       scaling；虚部の回転 so(d−1)；C のみ位相回転（可換なので一方掛けが対称）。"""
    N = S.shape[0]; vecs = [S.reshape(-1)]  # scaling
    if d == 2:
        rot = np.column_stack([-S[:, 1], S[:, 0]])  # z → i z
        vecs.append(rot.reshape(-1))
    else:
        for i in range(1, d):
            for j in range(i + 1, d):
                T = np.zeros_like(S); T[:, i] = -S[:, j]; T[:, j] = S[:, i]
                vecs.append(T.reshape(-1))
    return np.array(vecs)

def rank(M, tol=1e-9):
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > tol * s[0]))

def dimension_table():
    table = []
    for d, name in [(2, "C"), (4, "H"), (8, "O")]:
        for N in range(3, 11):
            S = gen_closure_solution(d, N)
            rl, rs = closure_residual(d, S)
            J = constraint_jacobian(d, S); rJ = rank(J)
            sol_dim = N * d - rJ
            T = symmetry_tangents(d, S); rT = rank(T)
            # 対称性の接ベクトルが解多様体の接空間に入ることの確認
            inJ = float(np.max(np.abs(J @ T.T)))
            img_rank = rank(S[:, 1:].T)  # 開いている虚方向の数
            table.append({"algebra": name, "d": d, "N": N, "residual_lin": rl, "residual_sq": rs,
                          "constraint_rank": rJ, "solution_dim": sol_dim, "symmetry_orbit_dim": rT,
                          "shape_dim": sol_dim - rT, "open_imaginary_directions": img_rank,
                          "symmetry_tangent_in_solution_tangent_maxabs": inJ})
    return table
results["dimension_table"] = dimension_table()
print("\nalgebra  N  constr_rank  sol_dim  orbit_dim  shape_dim  open_dirs")
for r in results["dimension_table"]:
    print(f"{r['algebra']:>3}  {r['N']:>2}  {r['constraint_rank']:>11}  {r['solution_dim']:>7}  {r['symmetry_orbit_dim']:>9}  {r['shape_dim']:>9}  {r['open_imaginary_directions']:>9}   res={r['residual_lin']:.1e},{r['residual_sq']:.1e}")

# ---------------------------------------------------------------------------
# 4. 関係 ō o′ の分解
# ---------------------------------------------------------------------------
def relation_decomposition(trials=1000):
    mx = [0.0, 0.0, 0.0, 0.0]
    for _ in range(trials):
        o1, o2 = rng.normal(size=(2, 8))
        a, b, c, d = o1[:4], o1[4:], o2[:4], o2[4:]
        r = omul(oconj(o1), o2)
        # 実部 = ⟨a,c⟩ + ⟨b,d⟩
        mx[0] = max(mx[0], abs(r[0] - (a @ c + b @ d)))
        # 四元数虚部 = Vec(ā c) + Vec(d̄ b)
        v = qmul(qconj(a), c)[1:] + qmul(qconj(d), b)[1:]
        mx[1] = max(mx[1], np.linalg.norm(r[1:4] - v))
        # ℓ 部 = d ā − b c̄
        l = qmul(d, qconj(a)) - qmul(b, qconj(c))
        mx[2] = max(mx[2], np.linalg.norm(r[4:] - l))
        # 自己関係の ℓ 部と虚部は 0
        rs = omul(oconj(o1), o1)
        mx[3] = max(mx[3], np.linalg.norm(rs[1:]))
    return {"realpart_eq_qdot_plus_pdot_err": mx[0], "quat_imag_eq_Vec(ac)+Vec(db)_err": mx[1],
            "l_part_eq_da-bc_err": mx[2], "self_relation_imag_err": mx[3]}
results["relation_decomposition"] = relation_decomposition()
print("\nrelation decomposition:", results["relation_decomposition"])

# ---------------------------------------------------------------------------
# 5. 閉塞が内部半分のスケール κ を決める
# ---------------------------------------------------------------------------
def conversion_ratio(N=12, trials=200):
    out = []
    for _ in range(trials):
        S = gen_closure_solution(8, N)             # (N, 8): 列 0..3 = q, 4..7 = p
        q, p = S[:, :4], S[:, 4:]
        # p を κ 倍：線形閉塞は保たれる。二乗閉塞の実部 Σ(q0² − |v|² − κ²|p|²) = 0 から κ² が決まる
        num = np.sum(q[:, 0]**2 - np.sum(q[:, 1:]**2, axis=1)); den = np.sum(p**2)
        kappa2 = num / den
        out.append(kappa2)
        # 確認：κ=1 でないスケールでは閉塞が破れること、κ² の式で戻ること
    out = np.array(out)
    # 具体例：ある解で p を 2 倍にして残差を見る → 式で直す
    S = gen_closure_solution(8, N); S2 = S.copy(); S2[:, 4:] *= 2.0
    rl, rs = closure_residual(8, S2)
    num = np.sum(S2[:, 0]**2 - np.sum(S2[:, 1:4]**2, axis=1)); den = np.sum(S2[:, 4:]**2)
    S3 = S2.copy(); S3[:, 4:] *= np.sqrt(num / den)
    rl3, rs3 = closure_residual(8, S3)
    return {"kappa2_of_generated_solutions_min_max": [float(out.min()), float(out.max())],
            "note": "生成解は κ=1 で作っているので κ²=1 が戻る。p を 2 倍にした配置の残差とその修復：",
            "residual_after_p_x2": [rl, rs], "residual_after_rescale_by_formula": [rl3, rs3]}
results["conversion_ratio"] = conversion_ratio()
print("\nconversion ratio:", results["conversion_ratio"])

# ---------------------------------------------------------------------------
# 6. 三角形の枠ホロノミー
#    比の一歩 s_ef = o_f ō_e（単位元）。枠 w ↦ s_ge(s_fg(s_ef w))。
#    H（結合的）では q 半分の枠は必ず戻る（telescoping）。O では非結合性のため戻らない。
#    どの三角形も自分の四元数部分代数 H_T = span(1, A, B, AB) を持ち、その中では枠は厳密に戻り、
#    補空間（4 次元）が回る。小さな三角形では回転角 = 4 × 面積（可換子 [A,B] = 2 A×B から）。
# ---------------------------------------------------------------------------
def unit(o): return o / np.sqrt(onorm2(o))

def frame_holonomy(oe, of, og):
    sef, sfg, sge = omul(of, oconj(oe)), omul(og, oconj(of)), omul(oe, oconj(og))
    H = np.zeros((8, 8))
    for k in range(8):
        H[:, k] = omul(sge, omul(sfg, omul(sef, E[k])))
    return H

def rot_angle(M):
    ev = np.linalg.eigvals(M); return float(np.max(np.abs(np.angle(ev))))

def sph_area(e, f, g):
    """三点を含む大 2 球面上の測地三角形の面積（球面過剰）。"""
    def ang(p, q, r):
        tq = q - (p @ q) * p; tr = r - (p @ r) * p
        return np.arccos(np.clip(tq @ tr / np.linalg.norm(tq) / np.linalg.norm(tr), -1, 1))
    return float(ang(e, f, g) + ang(f, g, e) + ang(g, e, f) - np.pi)

def im_unit(mask=None):
    v = rng.normal(size=8); v[0] = 0.0
    if mask is not None: v = v * mask
    return v / np.linalg.norm(v)

def triangle_subalgebra_basis(A, B):
    Bp = B - (A @ B) * A; Bp /= np.linalg.norm(Bp)
    C = omul(A, Bp); C[0] = 0.0; C /= np.linalg.norm(C)
    S = np.array([E[0], A, Bp, C]).T                      # 8x4
    Q, _ = np.linalg.qr(np.hstack([S, rng.normal(size=(8, 8))]))
    return S, Q[:, 4:]                                    # 部分代数の基底と補空間の基底

def holonomy_study():
    out = {}
    # (a) H に埋め込んだ三角形（p = 0）：q 半分は戻る、ℓ 半分（内部枠）が回る。角/面積 → 4
    rec = []
    qmask = np.array([0, 1, 1, 1, 0, 0, 0, 0], float)
    for eps in np.logspace(-2.5, 0.2, 14):
        for _ in range(6):
            A, B = im_unit(qmask), im_unit(qmask)
            oe, of, og = E[0], oexp(eps * A), oexp(eps * B)
            H = frame_holonomy(oe, of, og); area = sph_area(oe, of, og)
            rec.append({"eps": float(eps), "area": area, "Hqq_dev": float(np.linalg.norm(H[:4, :4] - np.eye(4))),
                        "mix": float(np.linalg.norm(H[:4, 4:]) + np.linalg.norm(H[4:, :4])),
                        "angle_ll": rot_angle(H[4:, 4:])})
    out["H_embedded"] = rec
    small = [r for r in rec if r["area"] < 1e-3]
    out["H_embedded_small_triangle_angle_over_area_mean"] = float(np.mean([r["angle_ll"] / r["area"] for r in small]))
    out["H_embedded_max_Hqq_dev"] = float(max(r["Hqq_dev"] for r in rec))
    out["H_embedded_max_mix"] = float(max(r["mix"] for r in rec))
    # (b) O の一般の小さな三角形：自分の部分代数の中では恒等、補空間は 4×面積で回る。
    #     観測者の時空半分 q の枠がどれだけ回るかは、q 半分と補空間の重なり（不整合）で決まる
    rec = []
    for _ in range(300):
        eps = 10 ** rng.uniform(-2.5, -1)
        A, B = im_unit(), im_unit()
        S, Cm = triangle_subalgebra_basis(A, B)
        oe, of, og = E[0], oexp(eps * A), oexp(eps * B)
        H = frame_holonomy(oe, of, og); area = sph_area(oe, of, og)
        Pq = np.zeros((8, 4)); Pq[:4, :4] = np.eye(4)        # 観測者の時空半分
        misalign = float(np.linalg.norm(Cm.T @ Pq) ** 2)       # q 半分のうち補空間に乗る分（0..3）
        rec.append({"eps": eps, "area": area,
                    "Hsub_dev": float(np.linalg.norm(S.T @ H @ S - np.eye(4))),
                    "mix_sub_comp": float(np.linalg.norm(S.T @ H @ Cm)),
                    "angle_comp_over_area": rot_angle(Cm.T @ H @ Cm) / area,
                    "Hqq_dev_over_area": float(np.linalg.norm(H[:4, :4] - np.eye(4)) / area),
                    "Hll_dev_over_area": float(np.linalg.norm(H[4:, 4:] - np.eye(4)) / area),
                    "qp_mix_over_area": float(np.linalg.norm(H[:4, 4:]) / area),
                    "misalignment": misalign})
    out["O_small_triangles"] = rec
    out["O_max_Hsub_dev"] = float(max(r["Hsub_dev"] for r in rec))
    out["O_max_mix_sub_comp"] = float(max(r["mix_sub_comp"] for r in rec))
    out["O_angle_comp_over_area_min_max"] = [float(min(r["angle_comp_over_area"] for r in rec)), float(max(r["angle_comp_over_area"] for r in rec))]
    # (c) 大きな三角形（単位八元数をランダムに）：枠が大きく回る。基点は戻る（逆元性）
    angles, base_err = [], []
    for _ in range(400):
        oe, of, og = [unit(rng.normal(size=8)) for _ in range(3)]
        H = frame_holonomy(oe, of, og)
        angles.append(rot_angle(H)); base_err.append(float(np.linalg.norm(H @ oe - oe)))
    out["O_random_unit_triangles_angle_median_max"] = [float(np.median(angles)), float(np.max(angles))]
    out["O_random_basepoint_return_max_err"] = float(np.max(base_err))
    # (d) 共通の四元数部分代数 span(1, e1, e4, e5) の三角形：部分代数の枠は戻る
    Bq = np.array([E[0], E[1], E[4], E[5]])
    dev = []
    for _ in range(100):
        os_ = [unit(rng.normal(size=4) @ Bq) for _ in range(3)]
        H = frame_holonomy(*os_); dev.append(float(np.linalg.norm(Bq @ H @ Bq.T - np.eye(4))))
    out["O_common_subalgebra_frame_return_max_dev"] = float(max(dev))
    return out
results["holonomy"] = holonomy_study()
h = results["holonomy"]
print("\nholonomy:")
print("  H-embedded: max|Hqq-I| = %.1e, max mix = %.1e, small-triangle angle(Hll)/area = %.4f"
      % (h["H_embedded_max_Hqq_dev"], h["H_embedded_max_mix"], h["H_embedded_small_triangle_angle_over_area_mean"]))
print("  O small triangles: max|H_sub-I| = %.1e, max mix(sub,comp) = %.1e, angle(comp)/area in [%.4f, %.4f]"
      % (h["O_max_Hsub_dev"], h["O_max_mix_sub_comp"], *h["O_angle_comp_over_area_min_max"]))
print("  O random unit triangles: angle median/max = %.3f / %.3f, basepoint return err = %.1e; common-subalgebra frame return dev = %.1e"
      % (*h["O_random_unit_triangles_angle_median_max"], h["O_random_basepoint_return_max_err"], h["O_common_subalgebra_frame_return_max_dev"]))

# ---------------------------------------------------------------------------
# 保存
# ---------------------------------------------------------------------------
with open(os.path.join(OUT, "results.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=1)

# 図
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = ["DejaVu Sans"]

# 図1：形の空間の次元
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for name, color in [("C", "#2457e6"), ("H", "#d9800f"), ("O", "#d1247a")]:
    rows = [r for r in results["dimension_table"] if r["algebra"] == name]
    ax.plot([r["N"] for r in rows], [r["shape_dim"] for r in rows], "o-", color=color, label=f"{name}: shape dim")
Ns = np.arange(3, 11)
ax.plot(Ns, 2 * (Ns - 3), "--", color="#2457e6", alpha=0.5, label="2(N−3)")
ax.plot(Ns, 4 * (Ns - 3), "--", color="#d9800f", alpha=0.5, label="4(N−3)")
ax.plot(Ns, np.maximum(0, 8 * Ns - 38), "--", color="#d1247a", alpha=0.5, label="8N−38 (N≥9)")
ax.set_xlabel("N"); ax.set_ylabel("dim (shape space)"); ax.set_title("Zero-closure shape-space dimension: C, H, O")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_shape_dimension.png"), dpi=160); fig.savefig(os.path.join(FIG, "fig1_shape_dimension.svg"))

# 図2：ホロノミー
fig, axs = plt.subplots(1, 3, figsize=(13, 3.9))
rec = h["H_embedded"]
axs[0].semilogx([r["area"] for r in rec], [r["angle_ll"] / r["area"] for r in rec], "o", ms=3, color="#d9800f")
axs[0].axhline(4.0, color="#888", ls="--", lw=1)
axs[0].set_xlabel("spherical area of triangle (spacetime half only)"); axs[0].set_ylabel("angle(internal frame rotation) / area")
axs[0].set_title("Triangle inside H: q-frame returns, internal frame rotates\n(ratio → 4 as area → 0)", fontsize=9)
rec = h["O_small_triangles"]
axs[1].plot([r["misalignment"] for r in rec], [r["Hqq_dev_over_area"] for r in rec], "o", ms=3, color="#d1247a", label="observer spacetime frame |Hqq−I|/area")
axs[1].plot([r["misalignment"] for r in rec], [r["qp_mix_over_area"] for r in rec], "o", ms=3, color="#2457e6", alpha=0.6, label="spacetime↔internal mixing /area")
axs[1].set_xlabel("misalignment of observer's spacetime half with triangle's own H_T (0..3)"); axs[1].set_ylabel("per unit area")
axs[1].set_title("Small triangles in O: holonomy seen by the observer\n(zero when aligned; complement always rotates by 4×area)", fontsize=9)
axs[1].legend(fontsize=7)
axs[2].plot([r["eps"] for r in rec], [r["angle_comp_over_area"] for r in rec], "o", ms=3, color="#444")
axs[2].set_xscale("log"); axs[2].set_ylim(3.9, 4.1)
axs[2].set_xlabel("ε (triangle size)"); axs[2].set_ylabel("angle(complement) / area")
axs[2].set_title("Universal: complement of the triangle's own H_T\nrotates by 4 × area", fontsize=9)
for a in axs: a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_frame_holonomy.png"), dpi=160); fig.savefig(os.path.join(FIG, "fig2_frame_holonomy.svg"))
print("\nsaved results and figures to", OUT, FIG)

# ---------------------------------------------------------------------------
# コンパクト SVG（Drive 保存用、依存なし・小サイズ。内容は上の図と同じデータ）
# ---------------------------------------------------------------------------
def compact_svg(path, panels, W=1080, Hh=330):
    """panels: list of dict(title, xlabel, ylabel, xlog, ylog, xlim, ylim, series=[(x,y,color,kind,label)], hlines=[...])"""
    n = len(panels); pw = W // n; m = dict(l=58, r=14, t=34, b=44)
    def tf(v, lo, hi, log):
        if log: v, lo, hi = np.log10(v), np.log10(lo), np.log10(hi)
        return (v - lo) / (hi - lo)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{Hh}" viewBox="0 0 {W} {Hh}" font-family="Helvetica, Arial, sans-serif" font-size="11">',
           f'<rect width="{W}" height="{Hh}" fill="white"/>']
    for i, p in enumerate(panels):
        x0 = i * pw + m["l"]; x1 = (i + 1) * pw - m["r"]; y0 = m["t"]; y1 = Hh - m["b"]
        out.append(f'<rect x="{x0}" y="{y0}" width="{x1-x0}" height="{y1-y0}" fill="none" stroke="#888"/>')
        out.append(f'<text x="{(x0+x1)/2:.0f}" y="{y0-16}" text-anchor="middle" font-size="12">{p["title"]}</text>')
        out.append(f'<text x="{(x0+x1)/2:.0f}" y="{Hh-10}" text-anchor="middle">{p["xlabel"]}</text>')
        out.append(f'<text transform="translate({i*pw+16},{(y0+y1)/2:.0f}) rotate(-90)" text-anchor="middle">{p["ylabel"]}</text>')
        xl, yl = p["xlim"], p["ylim"]
        # ticks
        xt = [10**k for k in range(int(np.floor(np.log10(xl[0]))), int(np.ceil(np.log10(xl[1])))+1)] if p.get("xlog") else np.linspace(xl[0], xl[1], 5)
        for v in xt:
            if v < xl[0] - 1e-12 or v > xl[1] + 1e-12: continue
            xx = x0 + tf(v, xl[0], xl[1], p.get("xlog")) * (x1 - x0)
            lab = (f'1e{int(np.log10(v))}' if p.get("xlog") else f'{v:g}')
            out.append(f'<line x1="{xx:.1f}" y1="{y1}" x2="{xx:.1f}" y2="{y1+4}" stroke="#888"/><text x="{xx:.1f}" y="{y1+16}" text-anchor="middle" font-size="10">{lab}</text>')
        for v in np.linspace(yl[0], yl[1], 5):
            yy = y1 - tf(v, yl[0], yl[1], False) * (y1 - y0)
            out.append(f'<line x1="{x0-4}" y1="{yy:.1f}" x2="{x0}" y2="{yy:.1f}" stroke="#888"/><text x="{x0-6}" y="{yy+3:.1f}" text-anchor="end" font-size="10">{v:g}</text>')
        for hv in p.get("hlines", []):
            yy = y1 - tf(hv, yl[0], yl[1], False) * (y1 - y0)
            out.append(f'<line x1="{x0}" y1="{yy:.1f}" x2="{x1}" y2="{yy:.1f}" stroke="#999" stroke-dasharray="4 3"/>')
        ly = y0 + 14
        for (xs, ys, color, kind, label) in p["series"]:
            pts = [(x0 + tf(x, xl[0], xl[1], p.get("xlog")) * (x1 - x0), y1 - tf(min(max(y, yl[0]), yl[1]), yl[0], yl[1], False) * (y1 - y0)) for x, y in zip(xs, ys)]
            if kind in ("line", "dash"):
                d = " ".join(f'{a:.1f},{b:.1f}' for a, b in pts)
                out.append(f'<polyline points="{d}" fill="none" stroke="{color}" stroke-width="1.6" {"stroke-dasharray=\"5 4\"" if kind=="dash" else ""}/>')
            if kind in ("line", "dots"):
                r = 3 if kind == "line" else 2.4
                out.append(f'<g fill="{color}" fill-opacity="0.85">' + "".join(f'<circle cx="{a:.1f}" cy="{b:.1f}" r="{r}"/>' for a, b in pts) + '</g>')
            if label:
                out.append(f'<circle cx="{x0+10}" cy="{ly}" r="3" fill="{color}"/><text x="{x0+17}" y="{ly+4}" font-size="10">{label}</text>'); ly += 14
    out.append('</svg>')
    with open(path, "w", encoding="utf-8") as f: f.write("\n".join(out))

Ns = list(range(3, 11))
sd = lambda name: [r["shape_dim"] for r in results["dimension_table"] if r["algebra"] == name]
compact_svg(os.path.join(FIG, "fig1_shape_dimension_compact.svg"), [dict(
    title="Zero-closure shape-space dimension", xlabel="N", ylabel="dim(shape space)", xlim=(3, 10), ylim=(0, 45),
    series=[(Ns, sd("C"), "#2457e6", "line", "C (= 2(N-3))"), (Ns, sd("H"), "#d9800f", "line", "H (= 4(N-3), N>=4)"),
            (Ns, sd("O"), "#d1247a", "line", "O (= 8N-38, N>=9)"), (Ns, [max(0, 8*n-38) for n in Ns], "#d1247a", "dash", "")])], W=560, Hh=360)

rngc = np.random.default_rng(0)
recH = h["H_embedded"]; recO = h["O_small_triangles"]
idx = rngc.choice(len(recO), 120, replace=False); recOs = [recO[i] for i in idx]
compact_svg(os.path.join(FIG, "fig2_frame_holonomy_compact.svg"), [
    dict(title="Triangle in H: internal frame rotation angle / area", xlabel="spherical area of triangle", ylabel="angle / area",
         xlog=True, xlim=(1e-6, 10), ylim=(0, 4.5), hlines=[4.0],
         series=[([r["area"] for r in recH], [r["angle_ll"]/r["area"] for r in recH], "#d9800f", "dots", "q-frame returns exactly; ratio -> 4 for small triangles")]),
    dict(title="Small triangles in O: what the observer sees (per unit area)", xlabel="misalignment of observer's spacetime half with triangle's H_T", ylabel="per unit area",
         xlim=(0, 3), ylim=(0, 6),
         series=[([r["misalignment"] for r in recOs], [r["Hqq_dev_over_area"] for r in recOs], "#d1247a", "dots", "|Hqq - I| / area (spacetime frame)"),
                 ([r["misalignment"] for r in recOs], [r["qp_mix_over_area"] for r in recOs], "#2457e6", "dots", "spacetime-internal mixing / area")]),
    dict(title="Complement of the triangle's own H_T rotates by 4 x area", xlabel="epsilon (triangle size)", ylabel="angle(complement) / area",
         xlog=True, xlim=(1e-3, 0.2), ylim=(3.9, 4.1), hlines=[4.0],
         series=[([r["eps"] for r in recOs], [r["angle_comp_over_area"] for r in recOs], "#444444", "dots", "")])])
print("compact SVGs written")
