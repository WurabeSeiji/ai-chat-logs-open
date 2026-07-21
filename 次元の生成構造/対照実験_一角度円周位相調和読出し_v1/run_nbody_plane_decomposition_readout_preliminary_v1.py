"""N体固定生成子系・平面分解読出し予備実験 v1（第4論文用）

固定生成子のスペクトル分解で得られる各回転平面について、
前シリーズのAB二体読出し（位置・加速度様・エネルギー様）が
平面ごとに成立するかを検査する。

実験F: 平面射影のAB運動学的同型性
  各平面の符号付き位相進行が単一周波数（1ステップ角 = 2 atan(γσ_j)）で
  進み、平面振幅が一定であること。基底の向きは Kp_j = σ_j q_j で固定し、
  面内 SO(2) ゲージで符号付き位相進行が不変であること（絶対位相は
  読めない＝公理5整合。反射を含む O(2) では符号は反転する）。
  名称置換共変性。周波数非縮退（独立の一般位置仮定）の検査。

実験G: 零次調和二階構造の平面確認
  各平面座標の離散二階差分が Δ²c = -4 sin²(θ_j/2) c を満たすこと。
  最遅平面を基準とする連続スペクトル再表示（比は一般に非整数）で、
  調和係数 g_j のセル幅に対する指数が -2 に一致すること
  （定義に基づく帰結の確認。整数倍音閉鎖や公理17の観測加速度
  射影の検証ではない）。

実験H: エネルギー様分解と閉鎖分解
  H_j と H_ker が個別に保存され Σ H_j + H_ker = H が成立すること。
  複素双線形の閉鎖量 Q_j = (p_j^T X)² + (q_j^T X)² と Q_ker も
  個別に保存され Σ Q_j + Q_ker = R² が成立すること
  （固定生成子では閉鎖量もエネルギー様量も平面間を移動しない）。

実験I: 残余部分空間の読出し制限
  核成分は位相進行を持たず（静的）、進行に基づく読出しが
  成立しないこと。射影子・核二乗量は一意に読め、核内部の
  標準基底は K だけからは選べないこと（帰結16.3の動力学的実証）。
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

TAU = 2.0 * math.pi
BASE_DIR = Path(__file__).resolve().parent
RESULT_DIR = BASE_DIR / "nbody_plane_decomposition_readout_result_v1"

BODY_COUNTS = (3, 4, 5, 6, 7, 8, 9)
TRIAL_COUNT = 32
STEP_COUNT = 720
TARGET_PERIOD_STEPS = 144
RADIUS_SQUARED = 1.0
IMAG_SEED = 0.35
SEED = 20260731
TOL = 1.0e-10
RANK_TOL = 1.0e-10
AMP_FLOOR = 1.0e-8


def relation_pairs(body_count: int) -> List[Tuple[int, int]]:
    return [(i, j) for i in range(body_count) for j in range(i + 1, body_count)]


def relation_adjacency(pairs: Sequence[Tuple[int, int]]) -> np.ndarray:
    count = len(pairs)
    matrix = np.zeros((count, count))
    for a in range(count):
        for b in range(count):
            if a != b and set(pairs[a]) & set(pairs[b]):
                matrix[a, b] = 1.0
    return matrix


def initial_state(relation_count: int, rng: np.random.Generator) -> np.ndarray:
    u = rng.normal(size=relation_count)
    u /= np.linalg.norm(u)
    v = rng.normal(size=relation_count)
    v -= float(np.dot(v, u)) * u
    v /= np.linalg.norm(v)
    return math.sqrt(RADIUS_SQUARED + IMAG_SEED**2) * u + 1j * IMAG_SEED * v


def build_generator(state: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    phases = np.angle(state)
    raw = adjacency * np.sin(phases[np.newaxis, :] - phases[:, np.newaxis])
    raw = 0.5 * (raw - raw.T)
    norm = float(np.linalg.norm(raw, ord=2))
    return raw / norm if norm > 1.0e-14 else np.zeros_like(raw)


def cayley(generator: np.ndarray) -> Tuple[np.ndarray, float]:
    dim = generator.shape[0]
    gamma = math.tan((TAU / TARGET_PERIOD_STEPS) / 2.0)
    identity = np.eye(dim)
    update = np.linalg.solve(identity - gamma * generator, identity + gamma * generator)
    return update, gamma


def plane_decomposition(generator: np.ndarray, rank_tol: float):
    """実反対称生成子の回転平面（正規直交基底 p,q と周波数 σ）と核基底。"""
    eigvals, eigvecs = np.linalg.eig(generator)
    scale = max(1.0, float(np.max(np.abs(eigvals))) if eigvals.size else 0.0)
    planes = []
    used = np.zeros(len(eigvals), dtype=bool)
    for idx in range(len(eigvals)):
        if used[idx]:
            continue
        sigma = float(eigvals[idx].imag)
        if sigma <= rank_tol * scale:
            continue
        w = eigvecs[:, idx]
        p = np.real(w)
        q = np.imag(w)
        # 正規直交化（正規行列なので p⊥q, |p|=|q| が成立するはず）
        p = p / np.linalg.norm(p)
        q = q - float(np.dot(q, p)) * p
        q = q / np.linalg.norm(q)
        # 向き規約: K p = +σ q（本文第3.3節と一致させる）
        if float(np.dot(generator @ p, q)) < 0.0:
            q = -q
        planes.append({"sigma": sigma, "p": p, "q": q})
        used[idx] = True
        # 共役対を消費
        conj_idx = int(np.argmin(np.abs(eigvals - np.conj(eigvals[idx]))))
        used[conj_idx] = True
    planes.sort(key=lambda item: item["sigma"])
    # 核基底
    _, s_values, vh = np.linalg.svd(generator)
    rank = int(np.sum(s_values > rank_tol * max(1.0, s_values[0])))
    kernel_basis = vh[rank:].T
    return planes, kernel_basis


def phase_series(states: np.ndarray, p: np.ndarray, q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """実部の平面座標 ζ = c + i d の系列と1ステップ位相進行。"""
    c = states.real @ p
    d = states.real @ q
    zeta = c + 1j * d
    increments = np.angle(zeta[1:] / zeta[:-1])
    return zeta, increments


def trial_run(body_count: int, rng: np.random.Generator) -> Dict[str, Any]:
    pairs = relation_pairs(body_count)
    relation_count = len(pairs)
    adjacency = relation_adjacency(pairs)
    state0 = initial_state(relation_count, rng)
    generator = build_generator(state0, adjacency)
    update, gamma = cayley(generator)

    states = np.empty((STEP_COUNT + 1, relation_count), dtype=complex)
    states[0] = state0
    for s in range(STEP_COUNT):
        states[s + 1] = update @ states[s]

    planes, kernel_basis = plane_decomposition(generator, RANK_TOL)
    theta_theory = [2.0 * math.atan(gamma * pl["sigma"]) for pl in planes]

    # 向き規約の検査: (Kp)·q / σ = +1
    orientation_dev = max(
        (abs(float(np.dot(generator @ pl["p"], pl["q"])) / pl["sigma"] - 1.0)
         for pl in planes), default=0.0
    )
    # 周波数非縮退（独立の一般位置仮定）の検査
    sigmas = sorted(pl["sigma"] for pl in planes)
    min_freq_gap = min(
        (b - a for a, b in zip(sigmas, sigmas[1:])), default=float("inf")
    )
    frequencies_distinct = bool(min_freq_gap > RANK_TOL)

    # ---------------- 実験F ----------------
    max_amp_drift = 0.0
    max_increment_spread = 0.0
    max_theory_dev = 0.0
    max_gauge_dev = 0.0
    skipped_planes = 0
    plane_H = []
    for pl, theta in zip(planes, theta_theory):
        zeta, increments = phase_series(states, pl["p"], pl["q"])
        amp = np.abs(zeta)
        if float(np.min(amp)) < AMP_FLOOR:
            skipped_planes += 1
            continue
        max_amp_drift = max(max_amp_drift, float(np.max(np.abs(amp - amp[0]))))
        max_increment_spread = max(
            max_increment_spread, float(np.max(np.abs(increments - increments[0])))
        )
        max_theory_dev = max(
            max_theory_dev, abs(float(np.mean(increments)) - theta)
        )
        # 面内ゲージ（絶対位相の不可読・進行の可読）
        beta = rng.uniform(0.0, TAU)
        p_g = math.cos(beta) * pl["p"] + math.sin(beta) * pl["q"]
        q_g = -math.sin(beta) * pl["p"] + math.cos(beta) * pl["q"]
        _, increments_g = phase_series(states, p_g, q_g)
        max_gauge_dev = max(
            max_gauge_dev, float(np.max(np.abs(increments_g - increments)))
        )
        # エネルギー様（複素両成分）
        h_series = np.abs(states @ pl["p"]) ** 2 + np.abs(states @ pl["q"]) ** 2
        plane_H.append(h_series)

    # 名称置換共変性: 周波数集合とエネルギー集合の不変性
    perm = rng.permutation(body_count)
    index_by_pair = {tuple(sorted(pr)): i for i, pr in enumerate(pairs)}
    edge_perm = np.zeros((relation_count, relation_count))
    for old, pr in enumerate(pairs):
        new_pair = tuple(sorted((int(perm[pr[0]]), int(perm[pr[1]]))))
        edge_perm[index_by_pair[new_pair], old] = 1.0
    generator_p = build_generator(edge_perm @ state0, adjacency)
    sig = np.sort([pl["sigma"] for pl in planes])
    eig_p = np.linalg.eigvals(generator_p)
    sig_p = np.sort([v.imag for v in eig_p if v.imag > RANK_TOL])
    perm_freq_dev = float(np.max(np.abs(sig - sig_p))) if len(sig) == len(sig_p) else float("inf")

    planes_p, _ = plane_decomposition(generator_p, RANK_TOL)
    state0_p = edge_perm @ state0
    H0 = np.sort([
        float(abs(np.dot(pl["p"], state0_p.real)) ** 2 + abs(np.dot(pl["p"], state0_p.imag)) ** 2
              + abs(np.dot(pl["q"], state0_p.real)) ** 2 + abs(np.dot(pl["q"], state0_p.imag)) ** 2)
        for pl in planes_p
    ])
    H0_orig = np.sort([float(h[0]) for h in plane_H]) if plane_H else np.array([])
    perm_energy_dev = (
        float(np.max(np.abs(H0 - H0_orig)))
        if len(H0) == len(H0_orig) and len(H0) > 0 else 0.0
    )

    # ---------------- 実験G ----------------
    # 判定は絶対偏差（理論係数 4sin²(θ/2) が微小な最遅平面では、
    # 機械精度の絶対誤差が相対値で増幅されるため）。相対値は参考記録。
    max_slope_abs_dev = 0.0
    max_slope_rel_dev = 0.0
    for pl, theta in zip(planes, theta_theory):
        c_series = states.real @ pl["p"]
        second = c_series[2:] - 2.0 * c_series[1:-1] + c_series[:-2]
        mid = c_series[1:-1]
        denom = float(np.dot(mid, mid))
        if denom < AMP_FLOOR**2:
            continue
        slope = float(np.dot(mid, second)) / denom
        theory = -4.0 * math.sin(theta / 2.0) ** 2
        max_slope_abs_dev = max(max_slope_abs_dev, abs(slope - theory))
        max_slope_rel_dev = max(max_slope_rel_dev, abs(slope - theory) / abs(theory))

    # 調和表示の距離指数（平面が2枚以上のとき）
    exponent = None
    if len(theta_theory) >= 2:
        theta_min = min(theta_theory)
        log_dtheta = [math.log(TAU / (t / theta_min)) for t in theta_theory]
        log_alpha = [math.log(4.0 * math.sin(t / 2.0) ** 2) for t in theta_theory]
        exponent = float(np.polyfit(log_dtheta, log_alpha, 1)[0])

    # ---------------- 実験H ----------------
    # 閉鎖分解 Σ Q_j + Q_ker = R²（複素双線形）
    Q_total = np.sum(states ** 2, axis=1)
    Q_planes = []
    for pl in planes:
        cq = states @ pl["p"]
        dq = states @ pl["q"]
        Q_planes.append(cq ** 2 + dq ** 2)
    if kernel_basis.shape[1] > 0:
        Q_ker = np.sum((states @ kernel_basis) ** 2, axis=1)
    else:
        Q_ker = np.zeros(STEP_COUNT + 1, dtype=complex)
    max_planeQ_drift = max(
        (float(np.max(np.abs(q_series - q_series[0]))) for q_series in Q_planes),
        default=0.0,
    )
    max_kerQ_drift = float(np.max(np.abs(Q_ker - Q_ker[0])))
    Q_sum = (np.sum(Q_planes, axis=0) if Q_planes else 0.0) + Q_ker
    closure_decomposition_err = float(np.max(np.abs(Q_sum - Q_total)))

    H_total = np.sum(np.abs(states) ** 2, axis=1)
    max_planeH_drift = max(
        (float(np.max(np.abs(h - h[0]))) for h in plane_H), default=0.0
    )
    if kernel_basis.shape[1] > 0:
        kernel_proj = states @ kernel_basis
        H_ker = np.sum(np.abs(kernel_proj) ** 2, axis=1)
    else:
        kernel_proj = None
        H_ker = np.zeros(STEP_COUNT + 1)
    if plane_H:
        H_sum = np.sum(plane_H, axis=0) + H_ker
        decomposition_identity_err = float(np.max(np.abs(H_sum - H_total)))
    else:
        decomposition_identity_err = float(np.max(np.abs(H_ker - H_total)))
    max_kerH_drift = float(np.max(np.abs(H_ker - H_ker[0])))

    # ---------------- 実験I ----------------
    if kernel_proj is not None and kernel_basis.shape[1] > 0:
        kernel_static_err = float(np.max(np.abs(kernel_proj - kernel_proj[0])))
        nullity = kernel_basis.shape[1]
        g1 = np.linalg.qr(rng.normal(size=(nullity, nullity)))[0]
        g2 = np.linalg.qr(rng.normal(size=(nullity, nullity)))[0]
        b1 = kernel_basis @ g1
        b2 = kernel_basis @ g2
        kernel_projector_gauge_err = float(
            np.max(np.abs(b1 @ b1.T - b2 @ b2.T))
        )
        kernel_direction_gauge_angle = math.degrees(
            math.acos(min(1.0, abs(float(np.dot(b1[:, 0], b2[:, 0])))))
        ) if nullity >= 2 else 0.0
    else:
        kernel_static_err = 0.0
        kernel_projector_gauge_err = 0.0
        kernel_direction_gauge_angle = 0.0

    return {
        "body_count": body_count,
        "relation_count": relation_count,
        "plane_count": len(planes),
        "kernel_dim": int(kernel_basis.shape[1]),
        "rank_law_holds": bool(
            len(planes) == min(body_count, relation_count // 2)
            and kernel_basis.shape[1] == relation_count - 2 * len(planes)
        ),
        "orientation_dev": orientation_dev,
        "min_freq_gap": min_freq_gap,
        "frequencies_distinct": frequencies_distinct,
        "skipped_low_amp_planes": skipped_planes,
        "f_max_amp_drift": max_amp_drift,
        "f_max_increment_spread": max_increment_spread,
        "f_max_theory_dev": max_theory_dev,
        "f_max_gauge_dev": max_gauge_dev,
        "f_perm_freq_dev": perm_freq_dev,
        "f_perm_energy_dev": perm_energy_dev,
        "g_max_slope_abs_dev": max_slope_abs_dev,
        "g_max_slope_rel_dev": max_slope_rel_dev,
        "g_exponent": exponent,
        "h_max_planeQ_drift": max_planeQ_drift,
        "h_max_kerQ_drift": max_kerQ_drift,
        "h_closure_decomposition_err": closure_decomposition_err,
        "h_max_planeH_drift": max_planeH_drift,
        "h_max_kerH_drift": max_kerH_drift,
        "h_decomposition_identity_err": decomposition_identity_err,
        "i_kernel_static_err": kernel_static_err,
        "i_kernel_projector_gauge_err": kernel_projector_gauge_err,
        "i_kernel_direction_gauge_angle_deg": kernel_direction_gauge_angle,
    }


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    rows: List[Dict[str, Any]] = []
    for body_count in BODY_COUNTS:
        for _ in range(TRIAL_COUNT):
            rows.append(trial_run(body_count, rng))

    with (RESULT_DIR / "plane_decomposition_trials_v1.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summaries = []
    for body_count in BODY_COUNTS:
        sel = [r for r in rows if r["body_count"] == body_count]
        exps = [r["g_exponent"] for r in sel if r["g_exponent"] is not None]
        summary = {
            "body_count": body_count,
            "relation_count": sel[0]["relation_count"],
            "plane_count_values": sorted({r["plane_count"] for r in sel}),
            "kernel_dim_values": sorted({r["kernel_dim"] for r in sel}),
            "skipped_low_amp_total": sum(r["skipped_low_amp_planes"] for r in sel),
            "F_max_amp_drift": max(r["f_max_amp_drift"] for r in sel),
            "F_max_increment_spread": max(r["f_max_increment_spread"] for r in sel),
            "F_max_theory_dev": max(r["f_max_theory_dev"] for r in sel),
            "F_max_gauge_dev": max(r["f_max_gauge_dev"] for r in sel),
            "F_max_orientation_dev": max(r["orientation_dev"] for r in sel),
            "F_all_rank_law": all(r["rank_law_holds"] for r in sel),
            "F_all_frequencies_distinct": all(r["frequencies_distinct"] for r in sel),
            "F_min_freq_gap": min(r["min_freq_gap"] for r in sel),
            "F_max_perm_freq_dev": max(r["f_perm_freq_dev"] for r in sel),
            "F_max_perm_energy_dev": max(r["f_perm_energy_dev"] for r in sel),
            "G_max_slope_abs_dev": max(r["g_max_slope_abs_dev"] for r in sel),
            "G_max_slope_rel_dev": max(r["g_max_slope_rel_dev"] for r in sel),
            "G_exponent_min": min(exps) if exps else None,
            "G_exponent_max": max(exps) if exps else None,
            "H_max_planeQ_drift": max(r["h_max_planeQ_drift"] for r in sel),
            "H_max_kerQ_drift": max(r["h_max_kerQ_drift"] for r in sel),
            "H_max_closure_decomp_err": max(r["h_closure_decomposition_err"] for r in sel),
            "H_max_planeH_drift": max(r["h_max_planeH_drift"] for r in sel),
            "H_max_kerH_drift": max(r["h_max_kerH_drift"] for r in sel),
            "H_max_identity_err": max(r["h_decomposition_identity_err"] for r in sel),
            "I_max_kernel_static_err": max(r["i_kernel_static_err"] for r in sel),
            "I_max_kernel_projector_gauge_err": max(
                r["i_kernel_projector_gauge_err"] for r in sel
            ),
            "I_min_kernel_direction_angle_deg": min(
                (r["i_kernel_direction_gauge_angle_deg"] for r in sel
                 if r["kernel_dim"] >= 2), default=0.0
            ),
        }
        exponent_ok = (
            summary["G_exponent_min"] is None
            or (-2.001 <= summary["G_exponent_min"]
                and summary["G_exponent_max"] <= -1.999)
        )
        summary["passed"] = bool(
            summary["F_max_amp_drift"] <= TOL
            and summary["F_max_increment_spread"] <= TOL
            and summary["F_max_theory_dev"] <= TOL
            and summary["F_max_gauge_dev"] <= TOL
            and summary["F_max_orientation_dev"] <= TOL
            and summary["F_all_rank_law"]
            and summary["F_all_frequencies_distinct"]
            and summary["F_max_perm_freq_dev"] <= TOL
            and summary["F_max_perm_energy_dev"] <= TOL
            and summary["G_max_slope_abs_dev"] <= TOL
            and exponent_ok
            and summary["H_max_planeQ_drift"] <= TOL
            and summary["H_max_kerQ_drift"] <= TOL
            and summary["H_max_closure_decomp_err"] <= TOL
            and summary["H_max_planeH_drift"] <= TOL
            and summary["H_max_kerH_drift"] <= TOL
            and summary["H_max_identity_err"] <= TOL
            and summary["I_max_kernel_static_err"] <= TOL
            and summary["I_max_kernel_projector_gauge_err"] <= TOL
        )
        summaries.append(summary)

    payload = {
        "parameters": {
            "body_counts": BODY_COUNTS,
            "trial_count": TRIAL_COUNT,
            "step_count": STEP_COUNT,
            "target_period_steps": TARGET_PERIOD_STEPS,
            "radius_squared": RADIUS_SQUARED,
            "imag_seed": IMAG_SEED,
            "seed": SEED,
            "tolerance": TOL,
        },
        "summaries": summaries,
        "all_passed": all(s["passed"] for s in summaries),
    }
    with (RESULT_DIR / "plane_decomposition_readout_result_v1.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)

    print("all_passed:", payload["all_passed"])
    for s in summaries:
        print(
            f"N={s['body_count']} planes={s['plane_count_values']} ker={s['kernel_dim_values']} "
            f"F(amp={s['F_max_amp_drift']:.2e} inc={s['F_max_increment_spread']:.2e} "
            f"th={s['F_max_theory_dev']:.2e} gauge={s['F_max_gauge_dev']:.2e} "
            f"orient={s['F_max_orientation_dev']:.2e} gap={s['F_min_freq_gap']:.2e} "
            f"permF={s['F_max_perm_freq_dev']:.2e} permE={s['F_max_perm_energy_dev']:.2e}) "
            f"G(abs={s['G_max_slope_abs_dev']:.2e} rel={s['G_max_slope_rel_dev']:.2e} exp=[{s['G_exponent_min']},{s['G_exponent_max']}]) "
            f"H(pl={s['H_max_planeH_drift']:.2e} ker={s['H_max_kerH_drift']:.2e} "
            f"id={s['H_max_identity_err']:.2e} Qpl={s['H_max_planeQ_drift']:.2e} "
            f"Qid={s['H_max_closure_decomp_err']:.2e}) "
            f"I(static={s['I_max_kernel_static_err']:.2e} "
            f"proj={s['I_max_kernel_projector_gauge_err']:.2e} "
            f"angle>={s['I_min_kernel_direction_angle_deg']:.1f}) "
            f"pass={s['passed']}"
        )


if __name__ == "__main__":
    main()
