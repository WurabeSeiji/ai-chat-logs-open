#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
実験 II: 奇数 N の厳密周期結晶 (反対称波 F_{1,2}(u,v) = ω^{u+2v} - ω^{2u+v},
U^n = I かつ ΣT^2 = 0) に対する安定性解析と 1 マス欠陥の非線形発展。

対象:
  - N = 5 (周期 4, 125 成分): 1 周期 (4 ステップ) のモノドロミー行列 (125x125) の
    Floquet 乗数 (固有値) を厳密計算 + 1 マス欠陥の 200 ステップ非線形発展。
  - N = 7 (周期 3, 343 成分): 1 マス欠陥の 200 ステップ非線形発展。

乱数なし。3 次斉次スケーリング分離を用い、||P|| < 1e-12 を厳密ゼロとして停止する。
"""
import os
import numpy as np


def build_antisym_tensor(N, a=1, b=2):
    w = np.exp(2j * np.pi / N)
    T = np.zeros((N, N, N), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                u = (j - i) % N
                v = (k - i) % N
                T[i, j, k] = w ** ((a * u + b * v) % N) - w ** ((a * v + b * u) % N)
    return T


def bm_prod(T):
    return np.einsum("itk,ijt,tjk->ijk", T, T, T)


def bm_jacobian(T):
    N = T.shape[0]
    dim = N ** 3
    J = np.zeros((dim, dim), dtype=np.complex128)
    for idx in range(dim):
        E = np.zeros((N, N, N), dtype=np.complex128)
        E.flat[idx] = 1.0
        dP = (
            np.einsum("itk,ijt,tjk->ijk", E, T, T)
            + np.einsum("itk,ijt,tjk->ijk", T, E, T)
            + np.einsum("itk,ijt,tjk->ijk", T, T, E)
        )
        J[:, idx] = dP.reshape(-1)
    return J


def proj_dist(A, B):
    na = np.linalg.norm(A)
    nb = np.linalg.norm(B)
    if na == 0 or nb == 0:
        return 1.0
    ov = abs(np.vdot(A / na, B / nb))
    return float(np.sqrt(max(0.0, 1.0 - ov * ov)))


def sum_z2_ratio(T):
    denom = float(np.sum(np.abs(T) ** 2))
    return 0.0 if denom == 0 else float(abs(np.sum(T ** 2)) / denom)


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("実験 II: 奇数 N (N=5, 7) 反対称波周期結晶の Floquet 安定性と 1 マス欠陥発展")
    out("=" * 72)

    # --- 1. N = 5 (周期 4) の Floquet 乗数 (125x125 モノドロミー行列) ---
    out("\n--- 1. N = 5 (周期 4, 125 成分) の 1 周期モノドロミー行列の固有値 ---")
    T_cycle = []
    T_cur = build_antisym_tensor(5, 1, 2)
    for st in range(4):
        T_cycle.append(T_cur)
        T_cur = bm_prod(T_cur) / 5.0
    assert np.allclose(T_cur, T_cycle[0], atol=1e-13), "N=5 周期 4 検算失敗"
    out("[検算合格] N = 5 反対称波は 3 次元全成分で厳密に周期 4 (誤差 < 1e-13)。")

    # M = (J_3 / 5) (J_2 / 5) (J_1 / 5) (J_0 / 5)
    M = np.eye(125, dtype=np.complex128)
    for st in range(4):
        J_k = bm_jacobian(T_cycle[st]) / 5.0
        M = J_k @ M
    evals = np.linalg.eigvals(M)
    # 絶対値で分類
    n_gt1 = int(np.sum(np.abs(evals) > 1.0 + 1e-6))
    n_eq1 = int(np.sum(np.abs(np.abs(evals) - 1.0) <= 1e-6))
    n_lt1 = int(np.sum((np.abs(evals) < 1.0 - 1e-6) & (np.abs(evals) > 1e-6)))
    n_zero = int(np.sum(np.abs(evals) <= 1e-6))
    out(f"Floquet 乗数 (全 125 本) の分布:")
    out(f"  |μ| > 1 (不安定): {n_gt1} 本  (うち 1 本はスカラー倍方向の自明な 3^4 = 81)")
    out(f"  |μ| = 1 (中立)  : {n_eq1} 本")
    out(f"  0 < |μ| < 1 (安定): {n_lt1} 本")
    out(f"  |μ| = 0 (超安定): {n_zero} 本")

    # 非ゼロ固有値を絶対値と角度で丸めて集計
    groups = {}
    for ev in evals:
        mag = abs(ev)
        if mag <= 1e-6:
            key = (0.0, 0.0)
        else:
            ang = round(float(np.degrees(np.angle(ev))), 2)
            if abs(ang + 180.0) < 0.05:
                ang = 180.0
            key = (round(mag, 5), ang)
        groups[key] = groups.get(key, 0) + 1
    out("  詳細 (|μ|, 偏角[度], 重複度):")
    for (mag, ang), cnt in sorted(groups.items(), key=lambda x: (-x[0][0], x[0][1])):
        out(f"    |μ| = {mag:9.5f}, 角 = {ang:+7.2f}°  × {cnt} 重")

    # --- 2. N = 5 および N = 7 の 1 マス欠陥 200 ステップ非線形発展 ---
    for N, period in [(5, 4), (7, 3)]:
        out(f"\n--- 2. N = {N} (周期 {period}, {N**3} 成分) の 1 マス欠陥 200 ステップ発展 ---")
        T0 = build_antisym_tensor(N, 1, 2)
        # 1 周期分の正規化基準状態
        cycle_units = []
        tc = T0.copy()
        for _ in range(period):
            cycle_units.append(tc / np.linalg.norm(tc))
            tc = bm_prod(tc) / float(N)

        defects = [
            ("対角 (0,0,0) に +ε", (0, 0, 0), 1e-4),
            ("対角 (0,0,0) に +iε", (0, 0, 0), 1e-4j),
            ("非対角 (0,1,0) に +ε", (0, 1, 0), 1e-4),
            ("非対角 (0,1,0) に +iε", (0, 1, 0), 1e-4j),
            ("R4 (0,1,2) に +ε", (0, 1, 2), 1e-4),
            ("R4 (0,1,2) に +iε", (0, 1, 2), 1e-4j),
        ]
        for label, pos, eps in defects:
            U = T0.copy()
            U[pos] += eps
            U = U / np.linalg.norm(U)
            hist = [U]
            status = "running"
            for step in range(1, 201):
                P = bm_prod(U)
                nrm = float(np.linalg.norm(P))
                if nrm < 1e-12:
                    status = f"zero at step {step}"
                    break
                U_next = P / nrm
                found_p = None
                for prev_i in range(max(0, len(hist) - 16), len(hist)):
                    if proj_dist(hist[prev_i], U_next) < 1e-10:
                        found_p = step - prev_i
                        break
                hist.append(U_next)
                U = U_next
                if found_p is not None and step >= 15:
                    d_cyc = min(proj_dist(cu, U) for cu in cycle_units)
                    sz2 = sum_z2_ratio(U)
                    nz_cnt = int(np.sum(np.abs(U) > 1e-5))
                    status = f"converged period={found_p} at step {step} (dist_cycle={d_cyc:.3e}, |Σz^2|/Σ|z|^2={sz2:.3e}, nz={nz_cnt}/{N**3})"
                    break
            if status == "running":
                d_cyc = min(proj_dist(cu, U) for cu in cycle_units)
                min_rec = min(proj_dist(h, U) for h in hist[-51:-1])
                sz2 = sum_z2_ratio(U)
                nz_cnt = int(np.sum(np.abs(U) > 1e-5))
                status = f"non-converged at step 200 (dist_cycle={d_cyc:.3f}, min_rec_50={min_rec:.3e}, |Σz^2|/Σ|z|^2={sz2:.3e}, nz={nz_cnt}/{N**3})"
            out(f"  [{label:24s}] → {status}")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_antisym_defect_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
