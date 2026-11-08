# -*- coding: utf-8 -*-
"""実装篇の中核関数（論文原案 v0.7 の第 II 部に、この関数の本文をそのまま載せる）。

関係は N × N = 16 の行列 Y_ij,m = ⟨q̄_{i,m} q_{j,m}⟩_ψ（左共役、倍音 m ごと、四元数値、i = j の自己関係を含む）。左だけの無名性（全体を左から g 倍）で Y は 4 成分とも不変。
対角から R_i²・偶奇・局在・N_eff、非対角から 縞の切片 Φ_ij・傾き Δ_ij・対の枠 (û_ij, ŵ_ij)・軸の角 n_i·n_j を読む。
作用は 16 の関係から生成子 A（4×4 四元数値、反エルミート）を作り U = exp(A) を右乗の 16×16 実行列として一度に掛ける（順序なし）。
実験の走行は含まない。"""
import math
import numpy as np
from scipy.linalg import expm

DEG = math.pi / 180.0

# ---------- 四元数（最後の軸が [w, x, y, z]）----------
def qmul(p, q):
    pw, px, py, pz = p[..., 0], p[..., 1], p[..., 2], p[..., 3]
    qw, qx, qy, qz = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([pw * qw - px * qx - py * qy - pz * qz,
                     pw * qx + px * qw + py * qz - pz * qy,
                     pw * qy - px * qz + py * qw + pz * qx,
                     pw * qz + px * qy - py * qx + pz * qw], axis=-1)

def qconj(q):
    return q * np.array([1.0, -1.0, -1.0, -1.0])

def unit(v):
    L = float(np.linalg.norm(v))
    return np.asarray(v, float) / L if L > 1e-15 else np.array([0.0, 0.0, 1.0])

def Rmat(h):
    """右乗の 4×4 実行列：vec(q h) = Rmat(h) @ vec(q)。列 k は e_k h"""
    E = np.eye(4)
    return np.stack([qmul(E[k], np.asarray(h, float)) for k in range(4)], axis=1)

# ---------- 状態：index_v9q4 のパルス ----------
def coefficients(phi_deg, M, odd=True, even=True):
    ms = [m for m in range(1, M + 1) if (m % 2 == 1 and odd) or (m % 2 == 0 and even)]
    cs = [(2.0 / math.pi) * math.sin(m * math.pi * phi_deg / 360.0) / m for m in ms]
    return np.array(ms), np.array(cs), phi_deg / 360.0

def make_pulse(psi, theta_deg, phi0_deg, phi_deg, M, A, n, odd=True, even=True):
    """q_i(ψ) = A [c₀ e^{n φ₀} + Σ_m c_m e^{n(φ₀ + m(ψ − θ))}]、n は単位純虚四元数（index_v9q4 と同じ式）"""
    ms, cs, c0 = coefficients(phi_deg, M, odd, even)
    beta = phi0_deg * DEG + np.outer(ms, psi - theta_deg * DEG)
    re = c0 * math.cos(phi0_deg * DEG) + np.sum(cs[:, None] * np.cos(beta), axis=0)
    im = c0 * math.sin(phi0_deg * DEG) + np.sum(cs[:, None] * np.sin(beta), axis=0)
    q = np.zeros((len(psi), 4)); q[:, 0] = A * re; q[:, 1:] = A * np.outer(im, unit(n))
    return q

# ---------- 読出し：倍音係数と 16 の関係 ----------
def harmonic_coeffs(q):
    """各成分の実フーリエ係数。q_m(ψ) = a_m cos mψ + b_m sin mψ（a_m, b_m は四元数、形 (NG/2+1, 4)）"""
    NG = q.shape[0]; F = np.fft.rfft(q, axis=0) / NG
    a = 2 * F.real; b = -2 * F.imag; a[0] /= 2
    return a, b

def relation_matrix(Q, Mmax):
    """Y[i, j, m] = ⟨q̄_{i,m} q_{j,m}⟩_ψ = (ā_{i,m} a_{j,m} + b̄_{i,m} b_{j,m})/2。形 (N, N, Mmax+1, 4)。i = j を含む 16 本。Y_ji = conj(Y_ij)"""
    N = len(Q); ab = [harmonic_coeffs(q) for q in Q]
    Y = np.zeros((N, N, Mmax + 1, 4))
    for i in range(N):
        ai, bi = ab[i]
        for j in range(N):
            aj, bj = ab[j]
            Y[i, j] = ((qmul(qconj(ai), aj) + qmul(qconj(bi), bj)) / 2)[:Mmax + 1]
            Y[i, j, 0] *= 2                                   # m = 0（直流）は ψ 平均が ā₀ a₀ で 1/2 が付かない
    return Y

def read_self(Yii, q, dpsi):
    """自己関係から読むもの：R² = Σ_m |X_ii,m|、偶奇の重み（対蹠像との重なり）、局在 L、有効倍音次数 N_eff"""
    pw = Yii[:, 0].copy(); R2 = float(pw.sum())
    dist = pw / R2 if R2 > 0 else pw
    n_eff = float(np.sum(np.arange(len(dist)) * dist))
    P = np.roll(q, q.shape[0] // 2, axis=0)
    eta = float(np.sum(qmul(qconj(q), P)[:, 0]) * dpsi / (np.sum(q * q) * dpsi))
    p = np.sum(q * q, axis=1); p = p / p.sum(); L = float(np.sum(p * p))
    return {"R2": R2, "N_eff": n_eff, "w_e": (1 + eta) / 2, "w_o": (1 - eta) / 2, "eta": eta, "L": L}

def pair_frame(Yij):
    """非対角関係のベクトル部の列 Vec_m（m ≥ 1）から対の枠を読む。
    Y_m = (c_i c_j/2)[(1 + c) cos ϑ_m − sin ϑ_m (n_i + n_j) − cos ϑ_m (n_i × n_j)]、c = n_i·n_j、ϑ_m = Φ − mΔ なので、
    列は (n_i + n_j) と (n_i × n_j) の張る面にある。主成分の第 1 軸 û ∝ n_i + n_j（符号は未定、orient_axes で全体の向きから決める）、
    第 2 軸 ŵ ∝ n_i × n_j。ŵ の向きは Vec·ŵ と Sc の符号が逆になるように取る（cos ϑ が両方に入る）。
    軸の角 c は (Vec·ŵ)/Sc = −√((1 − c)/(1 + c)) から読める（m によらない）。平面性 = 上位 2 固有値の割合（1 なら式どおり）"""
    V = Yij[1:, 1:]; Sc = Yij[1:, 0]
    C = V.T @ V; w, U = np.linalg.eigh(C)
    planarity = float((w[-1] + w[-2]) / max(np.trace(C), 1e-300))
    # ŵ：Vec_m·ŵ = −tan(γ/2) Sc_m が全 m で成り立つ唯一の方向。w_vec = argmin Σ_m (V_m·w + Sc_m)² = −C⁺ Σ_m Sc_m V_m は ŵ/tan(γ/2) なので
    # |w_vec| = 1/tan(γ/2)、c = cos γ = cos(2 arctan(1/|w_vec|))
    wvec = -np.linalg.pinv(C, rcond=1e-10) @ (V.T @ Sc)
    k = float(np.linalg.norm(wvec))
    if k > 1e-12:
        wv = wvec / k; c = float(np.cos(2 * math.atan(1.0 / k)))
        u = U[:, -1] - np.dot(U[:, -1], wv) * wv            # 面の中で ŵ に直交する方向が û（符号は未定）
        if np.linalg.norm(u) < 1e-9: u = U[:, -2] - np.dot(U[:, -2], wv) * wv
    else:                                                     # 軸共有：Vec_m ∝ sin ϑ_m n の一本線、ŵ なし、c = 1
        wv = np.zeros(3); c = 1.0; u = U[:, -1]
    return unit(u), (unit(wv) if k > 1e-12 else wv), c, planarity

def orient_axes(Uax):
    """対の軸 û_ij の符号を全体の向きで揃える。û は線（±）としてしか読めないので、全対の û ûᵀ の和（符号によらない）の主固有ベクトル r̂ を基準にし、
    û·r̂ > 0 に取る。r̂ 自体の符号は全体の鏡像の約束（4 体で決まらない唯一の離散量）。順序 i↔j で同じ û になる"""
    N = len(Uax); T = np.zeros((3, 3))
    for i in range(N):
        for j in range(N):
            if i != j and Uax[i][j] is not None: T += np.outer(Uax[i][j], Uax[i][j])
    r = np.linalg.eigh(T)[1][:, -1]
    if r[np.argmax(np.abs(r))] < 0: r = -r
    out = [[None] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            if i != j and Uax[i][j] is not None: out[i][j] = Uax[i][j] if np.dot(Uax[i][j], r) >= 0 else -Uax[i][j]
    return out, r

def read_pair(Yij, u, c):
    """縞の位相 ϑ_m = Φ − mΔ を m にわたって直線に当て、切片 Φ_ij と傾き −Δ_ij を読む。左共役の関係では tan ϑ = −s (Vec·û)/Sc、s = √((1 + c)/2)。
    û の符号を反転すると Φ, Δ の符号が反転する（対の向きと時計の向きは一組）"""
    Mmax = Yij.shape[0] - 1; ms = np.arange(1, Mmax + 1)
    Sc = Yij[1:, 0]; Vu = Yij[1:, 1:] @ u; s = math.sqrt(max((1 + c) / 2, 1e-12))
    th = np.unwrap(np.arctan2(-s * Vu, Sc)); w = np.sqrt(Sc * Sc + Vu * Vu); w = w / (w.sum() + 1e-300)
    keep = w > 1e-3 * w.max(); A_ = np.vstack([ms[keep], np.ones(int(keep.sum()))]).T; W = np.sqrt(w[keep])
    slope, icpt = np.linalg.lstsq(A_ * W[:, None], th[keep] * W, rcond=None)[0]
    resid = float(np.sqrt(np.sum(w[keep] * (th[keep] - (slope * ms[keep] + icpt)) ** 2)))
    Phi = (icpt + math.pi) % (2 * math.pi) - math.pi
    return {"Phi": float(Phi), "Delta": float(-slope), "resid": resid}

# ---------- 作用：16 の関係から生成子、同時ユニタリ ----------
def generator(N, Phi, U_axis):
    """A_ij = −(π + |Φ_ij|) û_ij/2（i ≠ j）、A_ii = −Σ_{j≠i} A_ij。|Φ| を使うのは順序 i↔j で A_ji = A_ij となるため（名前を使わない）。
    右乗の実表現で 4N×4N の反対称実行列 G を返す。各対の生成子は φ P_D^{(ij)} ⊗ R(û)、P_D^{(ij)} = (|i⟩−|j⟩)(⟨i|−⟨j|)/2 で、
    exp(φ P_D ⊗ R(û)) = P_S ⊗ I + P_D ⊗ R(e^{ûφ})、φ = π + Δ_F が 7/13・7/18 の二チャネル交換写像"""
    G = np.zeros((4 * N, 4 * N))
    for i in range(N):
        for j in range(i + 1, N):
            phi = math.pi + abs(Phi[i][j]); Ru = Rmat(np.concatenate([[0.0], U_axis[i][j]]))
            PD = np.zeros((N, N)); PD[i, i] = PD[j, j] = 0.5; PD[i, j] = PD[j, i] = -0.5
            G += phi * np.kron(PD, Ru)
    return G

def simultaneous_unitary(G):
    """U = exp(G)。G が反対称実行列なので U は直交（全ノルム保存）。順序は入らない"""
    return expm(G)

def apply_unitary(U, Q):
    """標本ごとに 4N ベクトルへ U を掛ける（ψ に依らない一つの作用素）"""
    N = len(Q); V = np.concatenate(Q, axis=1)           # (NG, 4N)
    W = V @ U.T
    return [W[:, 4 * i:4 * i + 4].copy() for i in range(N)]

def blocks_as_quaternions(U, N):
    """U の (i, j) ブロックは右乗 R(h_ij) になる（右乗の和と積は右乗）。h_ij = ブロックの第 1 列。検査用に残差も返す"""
    H = np.zeros((N, N, 4)); err = 0.0
    for i in range(N):
        for j in range(N):
            B = U[4 * i:4 * i + 4, 4 * j:4 * j + 4]; h = B[:, 0]; H[i, j] = h
            err = max(err, float(np.max(np.abs(B - Rmat(h)))))
    return H, err

def handshake_prediction(Ytot, H, i):
    """観測者 i の位相の回転の予測：⟨q_i, q_i′⟩ = ∫ q̄_i q_i′ = Σ_j ∫ q̄_i (q_j h_ij) = Σ_j Y_ij h_ij。
    Y_ii = R_i² は実数、Y_ij（j ≠ i）の大きさ |Y_ij|/R_i² が視野角（相手の大きさ × 重なり）。回転角は arg、軸は Vec の向き。
    16 の関係と作用素だけから書けるので、視野角が位相の回転にどう入るかがこの式で読める"""
    acc = np.zeros(4)
    for j in range(Ytot.shape[1]): acc += qmul(Ytot[i, j], H[i, j])
    ang = math.atan2(float(np.linalg.norm(acc[1:])), float(acc[0]))
    return ang, unit(acc[1:]) if np.linalg.norm(acc[1:]) > 1e-15 else np.zeros(3)

def rotation_measured(q0, q1, dpsi):
    ip = np.sum(qmul(qconj(q0), q1), axis=0) * dpsi
    ang = math.atan2(float(np.linalg.norm(ip[1:])), float(ip[0]))
    return ang, unit(ip[1:]) if np.linalg.norm(ip[1:]) > 1e-15 else np.zeros(3)

def gram_eigs(Ytot):
    """グラム行列 G_ij = Sc(Y_ij)/√(R_i² R_j²) の固有値（昇順）。閾値の下に沈む本数が有効階数の落ち＝凝縮体"""
    N = Ytot.shape[0]; G = np.zeros((N, N))
    for i in range(N):
        for j in range(N): G[i, j] = Ytot[i, j, 0] / math.sqrt(Ytot[i, i, 0] * Ytot[j, j, 0])
    return np.sort(np.linalg.eigvalsh(G))
