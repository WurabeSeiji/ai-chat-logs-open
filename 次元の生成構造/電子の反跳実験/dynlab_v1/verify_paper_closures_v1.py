#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""紙と鉛筆で閉じた3課題の保存検証 v1（DLdyn ノート 命題 Dyn-7〜Dyn-10 の正本）

W1  Dyn-7  θ̇ の閉形式: ẋ_{A,i} = v_i(A)λ̇_i/(2√λ_i) + √λ_i Σ_{j≠i} Ω_ij v_j(A)
           を、固有分解の直接差分と突き合わせる（一次精度）
W2  Dyn-8  N体調和閉鎖の恒等式: |ω_n|Δθ_n = 2πω_1（ω_n=nω_1, Δθ_n=2π/|n|）と
           α_n = RΩ²/Δθ_n² の log-log 勾配が厳密 −2 であること（算術恒等式）
W3  Dyn-9  線形部の変分性: スライス Cayley = 中点則（Crank–Nicolson）であり、
           h=iK（エルミート）の二次不変量 ⟨z,hz⟩・Σ|z|²・Σz² を厳密保存すること
W4  Dyn-10 頂点のハミルトン形式: R 固定の下で δz = i ∂H/∂z̄、
           H = g Σ_{e<e'} R_ee'( |z_e|²|z_e'|² − Re(z_e'² z̄_e²) )。
           勾配一致（Wirtinger 数値微分）と、凍結 R 流での H 保存（RK4）、
           および R を状態依存にしたときの H ドリフト（障害の実証）

全て決定的（固定シード）。結果: result_paper_closures_v1.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(2026)


# ----------------------------------------------------------------------
def w1_theta_dot_closed_form():
    n, k3 = 9, 3
    Q = RNG.normal(size=(n, n)); B = (Q + Q.T) / 2
    # ギャップを開ける（規約1域）
    B += np.diag(np.linspace(0, 4.0, n))
    P = RNG.normal(size=(n, n)); dB = (P + P.T) / 2
    eps = 1e-7

    lam, V = np.linalg.eigh(B)
    order = np.argsort(lam)[::-1]
    lam, V = lam[order], V[:, order]
    lam2, V2 = np.linalg.eigh(B + eps * dB)
    o2 = np.argsort(lam2)[::-1]
    lam2, V2 = lam2[o2], V2[:, o2]
    for i in range(n):  # 符号整列
        if V2[:, i] @ V[:, i] < 0:
            V2[:, i] = -V2[:, i]

    X1 = V[:, :k3] * np.sqrt(lam[:k3])[None, :]
    X2 = V2[:, :k3] * np.sqrt(lam2[:k3])[None, :]
    dX_num = (X2 - X1) / eps

    lam_dot = np.array([V[:, i] @ dB @ V[:, i] for i in range(n)])
    dX_th = np.zeros((n, k3))
    for i in range(k3):
        term = V[:, i] * lam_dot[i] / (2 * np.sqrt(lam[i]))
        for j in range(n):
            if j != i:
                om = (V[:, j] @ dB @ V[:, i]) / (lam[i] - lam[j])
                term = term + np.sqrt(lam[i]) * om * V[:, j]
        dX_th[:, i] = term
    return {"W1_max_err": float(np.max(np.abs(dX_num - dX_th)))}


# ----------------------------------------------------------------------
def w2_harmonic_closure_nbody():
    omega1 = np.pi / 72.0  # 実測の集団時計（S17）
    ns = np.arange(1, 33)
    dth = 2 * np.pi / ns
    om = ns * omega1
    duality = np.abs(om) * dth  # 全 n で 2π ω1 のはず
    Omega = 2 * np.pi * omega1
    Rc = 1.0
    alpha = Rc * om ** 2  # = R Ω² / Δθ²
    slope = np.polyfit(np.log(dth), np.log(alpha), 1)[0]
    return {
        "W2_duality_max_dev": float(np.max(np.abs(duality - Omega))),
        "W2_Omega_N": float(Omega),
        "W2_loglog_slope": float(slope),
        "W2_slope_minus_neg2": float(abs(slope + 2.0)),
    }


# ----------------------------------------------------------------------
def w3_linear_variational():
    n = 8
    A = RNG.normal(size=(n, n))
    K = A - A.T  # 実反対称
    h = 1j * K   # エルミート（検算）
    herm_dev = float(np.max(np.abs(h - h.conj().T)))
    tau = 0.3
    I = np.eye(n)
    O = np.linalg.solve(I - (tau / 2) * K, I + (tau / 2) * K)  # Cayley＝中点則
    orth_dev = float(np.max(np.abs(O.T @ O - I)))

    z = RNG.normal(size=n) + 1j * RNG.normal(size=n)
    invs = {"norm": [], "bilin": [], "energy": []}
    for _ in range(200):
        invs["norm"].append(float(np.vdot(z, z).real))
        invs["bilin"].append(complex(z @ z))
        invs["energy"].append(float((np.conj(z) @ (h @ z)).real))
        z = O @ z
    def drift(xs):
        xs = np.array(xs)
        return float(np.max(np.abs(xs - xs[0])))
    return {
        "W3_h_hermitian_dev": herm_dev,
        "W3_cayley_orthogonality": orth_dev,
        "W3_norm_drift": drift(invs["norm"]),
        "W3_bilinear_drift": drift(np.abs(np.array(invs["bilin"]) - invs["bilin"][0])),
        "W3_energy_drift": drift(invs["energy"]),
    }


# ----------------------------------------------------------------------
def _pairs(n):
    return [(a, b) for a in range(n) for b in range(a + 1, n)]


def w4_vertex_hamiltonian():
    # K_N の辺集合上の頂点（対対称 R_ee'。隣接＝頂点共有）
    nv = 5
    edges = _pairs(nv)
    m = len(edges)
    adj = [[j for j in range(m) if j != i and (set(edges[i]) & set(edges[j]))]
           for i in range(m)]
    Rp = {(i, j): float(RNG.random()) for i in range(m) for j in adj[i] if i < j}
    def Rget(i, j):
        return Rp[(i, j)] if i < j else Rp[(j, i)]

    g = 1.0

    def H(z):
        s = 0.0
        for i in range(m):
            for j in adj[i]:
                if i < j:
                    s += Rget(i, j) * (abs(z[i]) ** 2 * abs(z[j]) ** 2
                                       - np.real(z[j] ** 2 * np.conj(z[i]) ** 2))
        return g * s

    def rate(z):
        dz = np.zeros(m, complex)
        for i in range(m):
            acc = 0.0 + 0.0j
            for j in adj[i]:
                acc += Rget(i, j) * (abs(z[j]) ** 2 * z[i]
                                     - z[j] ** 2 * np.conj(z[i]))
            dz[i] = 1j * g * acc
        return dz

    z0 = RNG.normal(size=m) + 1j * RNG.normal(size=m)

    # (a) 勾配一致: δz_i = i ∂H/∂z̄_i（Wirtinger 数値微分）
    eps = 1e-6
    grad = np.zeros(m, complex)
    for i in range(m):
        zx = z0.copy(); zx[i] += eps
        zy = z0.copy(); zy[i] += 1j * eps
        dHx = (H(zx) - H(z0)) / eps
        dHy = (H(zy) - H(z0)) / eps
        grad[i] = 0.5 * (dHx + 1j * dHy)   # ∂H/∂z̄
    ga = float(np.max(np.abs(rate(z0) - 1j * grad)))

    # (b) 凍結 R 流の保存量（RK4）: H・Σ|z|²・Σz²
    def rk4(z, dt):
        k1 = rate(z); k2 = rate(z + dt / 2 * k1)
        k3 = rate(z + dt / 2 * k2); k4 = rate(z + dt * k3)
        return z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    z = z0.copy()
    H0, n0, c0 = H(z), float(np.vdot(z, z).real), complex(z @ z)
    dt, T = 5e-4, 400
    for _ in range(T):
        z = rk4(z, dt)
    hb = abs(H(z) - H0)
    nb = abs(float(np.vdot(z, z).real) - n0)
    cb = abs(complex(z @ z) - c0)

    # (c) 障害の実証: R を状態依存（例示則）にすると H は保存しない
    def rate_sd(z):
        dz = np.zeros(m, complex)
        for i in range(m):
            acc = 0.0 + 0.0j
            for j in adj[i]:
                Rij = Rget(i, j) * abs(z[i]) ** 2 / (1.0 + abs(z[i]) ** 2)
                acc += Rij * (abs(z[j]) ** 2 * z[i] - z[j] ** 2 * np.conj(z[i]))
            dz[i] = 1j * g * acc
        return dz

    def rk4sd(z, dt):
        k1 = rate_sd(z); k2 = rate_sd(z + dt / 2 * k1)
        k3 = rate_sd(z + dt / 2 * k2); k4 = rate_sd(z + dt * k3)
        return z + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    z = z0.copy(); H0sd = H(z)
    for _ in range(T):
        z = rk4sd(z, dt)
    hsd = abs(H(z) - H0sd)

    return {
        "W4a_rate_vs_igrad": ga,
        "W4b_H_drift_frozenR": hb,
        "W4b_norm_drift": nb,
        "W4b_closure_drift": cb,
        "W4c_H_drift_stateR": hsd,
    }


# ----------------------------------------------------------------------
def main():
    t0 = time.time()
    res = {"W1": w1_theta_dot_closed_form(), "W2": w2_harmonic_closure_nbody(),
           "W3": w3_linear_variational(), "W4": w4_vertex_hamiltonian()}
    res["elapsed_sec"] = time.time() - t0
    (HERE / "result_paper_closures_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False, default=str))
    print(json.dumps(res, indent=1, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
