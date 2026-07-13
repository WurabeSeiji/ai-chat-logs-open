from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "elastic_collision_eta_resolution_sweep_result_v2"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fermionic_reflection_q(q_initial: float, delta_f: float = math.pi) -> float:
    reflection_rate = math.sin(delta_f / 2.0) ** 2
    transmission_rate = math.cos(delta_f / 2.0) ** 2
    return q_initial * (transmission_rate - reflection_rate)


@dataclass
class Params:
    A_A: float = 1.0
    A_B: float = 1.0
    Nh_chi_A: int = 99
    Nh_chi_B: int = 99
    Nh_tau_A: int = 99
    Nh_tau_B: int = 99
    chi_A0: float = -0.2
    chi_B0: float = 0.2
    tau0: float = 0.2
    q_A0: int = 1
    q_B0: int = -1
    delta_s: float = 0.01
    v_chi: float = 1.0
    omega_A: float = 1.0
    omega_B: float = 1.0
    s_max: int = 10000
    purity_acceptance: float = 0.99
    tie_tol: float = 1e-9


@dataclass
class State:
    chi: float
    tau: float
    q: int
    amplitude: float


def eta_grid(sample_count: int) -> np.ndarray:
    return np.linspace(-math.pi, math.pi, sample_count, endpoint=False)


def eta_overlap(m_particle: int, m_read: int, grid: np.ndarray) -> complex:
    return complex(np.mean(np.exp(1j * m_particle * grid) * np.exp(-1j * m_read * grid)))


def observe_eta_modes(m_particle: int, candidate_modes: Iterable[int], grid: np.ndarray) -> Dict[int, complex]:
    return {mode: eta_overlap(m_particle, mode, grid) for mode in candidate_modes}


def detect_mode(observations: Dict[int, complex], tie_tol: float) -> int | str:
    ranked = sorted(observations.items(), key=lambda item: abs(item[1]), reverse=True)
    if len(ranked) >= 2 and abs(abs(ranked[0][1]) - abs(ranked[1][1])) <= tie_tol:
        return "ambiguous"
    return ranked[0][0]


def purity(observations: Dict[int, complex], target_mode: int) -> float:
    denom = sum(abs(value) for value in observations.values())
    if denom == 0:
        return 0.0
    return float(abs(observations[target_mode]) / denom)


def eps_chi_ab(params: Params) -> float:
    return math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)


def eps_tau_ab(params: Params) -> float:
    return math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)


def run_motion(params: Params) -> Dict[str, object]:
    a = State(params.chi_A0, -params.tau0, params.q_A0, params.A_A)
    b = State(params.chi_B0, -params.tau0, params.q_B0, params.A_B)
    eps_chi = eps_chi_ab(params)
    eps_tau = eps_tau_ab(params)

    step = 0
    collision_cell_reached = False
    while abs(a.chi - b.chi) >= eps_chi or abs(a.tau - b.tau) >= eps_tau:
        if step >= params.s_max:
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        step += 1
    else:
        collision_cell_reached = True

    q_a_before = a.q
    q_b_before = b.q
    if collision_cell_reached:
        a.q = fermionic_reflection_q(a.q)
        b.q = fermionic_reflection_q(b.q)

    post_step = 0
    post_collision_completed = False
    while abs(a.chi - b.chi) <= eps_chi or min(a.tau, b.tau) < params.tau0 - 1e-12:
        if post_step >= params.s_max:
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        post_step += 1
    else:
        post_collision_completed = True

    return {
        "final_A": a,
        "final_B": b,
        "collision_cell_reached": collision_cell_reached,
        "post_collision_propagation_completed": post_collision_completed,
        "collision_step": step,
        "post_collision_steps": post_step,
        "q_reversed": a.q == -q_a_before and b.q == -q_b_before,
        "separated_after_collision": a.chi < b.chi,
    }


def run_case(params: Params, eta_samples: int, m_a: int, m_b: int, motion: Dict[str, object]) -> Dict[str, object]:
    grid = eta_grid(eta_samples)
    modes = [m_a, m_b]

    obs_a0 = observe_eta_modes(m_a, modes, grid)
    obs_b0 = observe_eta_modes(m_b, modes, grid)
    obs_a2 = observe_eta_modes(m_a, modes, grid)
    obs_b2 = observe_eta_modes(m_b, modes, grid)

    detected_a0 = detect_mode(obs_a0, params.tie_tol)
    detected_b0 = detect_mode(obs_b0, params.tie_tol)
    detected_a2 = detect_mode(obs_a2, params.tie_tol)
    detected_b2 = detect_mode(obs_b2, params.tie_tol)

    purity_a0 = purity(obs_a0, m_a)
    purity_b0 = purity(obs_b0, m_b)
    purity_a2 = purity(obs_a2, m_a)
    purity_b2 = purity(obs_b2, m_b)
    min_purity = min(purity_a0, purity_b0, purity_a2, purity_b2)

    mode_difference = abs(m_b - m_a)
    alias_collision = mode_difference % eta_samples == 0
    label_ambiguous = "ambiguous" in {detected_a0, detected_b0, detected_a2, detected_b2}
    label_preserved = detected_a0 == m_a and detected_b0 == m_b and detected_a2 == m_a and detected_b2 == m_b
    purity_pass = min_purity >= params.purity_acceptance
    eta_resolution_valid = not alias_collision and not label_ambiguous and label_preserved and purity_pass

    motion_valid = (
        bool(motion["collision_cell_reached"])
        and bool(motion["post_collision_propagation_completed"])
        and bool(motion["q_reversed"])
        and bool(motion["separated_after_collision"])
    )

    return {
        "eta_samples": eta_samples,
        "m_A": m_a,
        "m_B": m_b,
        "mode_difference": mode_difference,
        "alias_collision": alias_collision,
        "detected_A_initial": detected_a0,
        "detected_B_initial": detected_b0,
        "detected_A_final": detected_a2,
        "detected_B_final": detected_b2,
        "label_ambiguous": label_ambiguous,
        "label_preserved": label_preserved,
        "purity_A_initial": purity_a0,
        "purity_B_initial": purity_b0,
        "purity_A_final": purity_a2,
        "purity_B_final": purity_b2,
        "min_purity": min_purity,
        "purity_pass": purity_pass,
        "eta_resolution_valid": eta_resolution_valid,
        "motion_valid": motion_valid,
        "case_valid": eta_resolution_valid and motion_valid,
    }


def serialise_rows(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for row in rows:
        item = dict(row)
        for key, value in list(item.items()):
            if isinstance(value, np.generic):
                item[key] = value.item()
        out.append(item)
    return out


def write_outputs(result: Dict[str, object]) -> None:
    rows = serialise_rows(result["cases"])
    (OUT_DIR / "eta_resolution_sweep_result_v2.json").write_text(
        json.dumps({**result, "cases": rows}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    csv_path = OUT_DIR / "eta_resolution_sweep_cases_v2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    report_lines = [
        "# Eta Resolution Sweep Result v1",
        "",
        "## Verdict",
        "",
        f"- total_cases: `{result['summary']['total_cases']}`",
        f"- valid_cases: `{result['summary']['valid_cases']}`",
        f"- invalid_cases: `{result['summary']['invalid_cases']}`",
        f"- alias_collision_cases: `{result['summary']['alias_collision_cases']}`",
        f"- non_alias_failures: `{result['summary']['non_alias_failures']}`",
        f"- first_eta_samples_all_pairs_valid: `{result['summary']['first_eta_samples_all_pairs_valid']}`",
        "",
        "## Interpretation",
        "",
        "The internal identification vibration is readable only up to the eta-readout sampling resolution.",
        "When the mode difference is a multiple of the eta sample count, the two label modes alias and the readout becomes ambiguous.",
        "This is a readout-resolution failure, not a failure of the collision map itself.",
        "",
        "## Cases",
        "",
        "| eta samples | m_A | m_B | diff | alias | min purity | detected A0 | detected B0 | valid |",
        "|---:|---:|---:|---:|---|---:|---|---|---|",
    ]
    for row in rows:
        report_lines.append(
            f"| {row['eta_samples']} | {row['m_A']} | {row['m_B']} | {row['mode_difference']} | "
            f"{row['alias_collision']} | {row['min_purity']:.12g} | {row['detected_A_initial']} | "
            f"{row['detected_B_initial']} | {row['case_valid']} |"
        )
    (OUT_DIR / "eta_resolution_sweep_report_v2.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    eta_values = result["eta_sample_counts"]
    pair_labels = [f"{m_a}-{m_b}" for m_a, m_b in result["mode_pairs"]]
    valid_matrix = np.zeros((len(result["mode_pairs"]), len(eta_values)), dtype=float)
    purity_matrix = np.zeros_like(valid_matrix)
    lookup = {(row["m_A"], row["m_B"], row["eta_samples"]): row for row in rows}
    for i, pair in enumerate(result["mode_pairs"]):
        pair_tuple = tuple(pair)
        for j, eta_samples in enumerate(eta_values):
            row = lookup[(pair_tuple[0], pair_tuple[1], eta_samples)]
            valid_matrix[i, j] = 1.0 if row["case_valid"] else 0.0
            purity_matrix[i, j] = float(row["min_purity"])

    fig, ax = plt.subplots(figsize=(10, 6))
    image = ax.imshow(valid_matrix, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks(range(len(eta_values)))
    ax.set_xticklabels([str(v) for v in eta_values])
    ax.set_yticks(range(len(pair_labels)))
    ax.set_yticklabels(pair_labels)
    ax.set_xlabel("eta sample count")
    ax.set_ylabel("label mode pair")
    ax.set_title("Eta readout resolution validity")
    fig.colorbar(image, ax=ax, ticks=[0, 1], label="case valid")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "eta_resolution_validity_v2.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    image = ax.imshow(purity_matrix, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(eta_values)))
    ax.set_xticklabels([str(v) for v in eta_values])
    ax.set_yticks(range(len(pair_labels)))
    ax.set_yticklabels(pair_labels)
    ax.set_xlabel("eta sample count")
    ax.set_ylabel("label mode pair")
    ax.set_title("Minimum identification purity")
    fig.colorbar(image, ax=ax, label="min purity")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "eta_resolution_purity_v2.png", dpi=180)
    plt.close(fig)


def run() -> Dict[str, object]:
    params = Params()
    eta_sample_counts = [4, 6, 8, 12, 16, 24, 32, 64]
    mode_pairs = [
        (1, 2),
        (1, 5),
        (1, 9),
        (1, 17),
        (1, 33),
        (8, 9),
        (15, 31),
        (16, 32),
        (17, 33),
        (5, 29),
        (7, 31),
    ]
    motion = run_motion(params)

    cases: List[Dict[str, object]] = []
    for eta_samples in eta_sample_counts:
        for m_a, m_b in mode_pairs:
            cases.append(run_case(params, eta_samples, m_a, m_b, motion))

    invalid = [row for row in cases if not row["case_valid"]]
    non_alias_failures = [row for row in invalid if not row["alias_collision"]]
    all_valid_eta = None
    for eta_samples in eta_sample_counts:
        eta_rows = [row for row in cases if row["eta_samples"] == eta_samples]
        if all(row["case_valid"] for row in eta_rows):
            all_valid_eta = eta_samples
            break

    return {
        "parameters": asdict(params),
        "eta_sample_counts": eta_sample_counts,
        "mode_pairs": [list(pair) for pair in mode_pairs],
        "motion": {
            key: value
            for key, value in motion.items()
            if key not in {"final_A", "final_B"}
        },
        "summary": {
            "total_cases": len(cases),
            "valid_cases": sum(1 for row in cases if row["case_valid"]),
            "invalid_cases": len(invalid),
            "alias_collision_cases": sum(1 for row in cases if row["alias_collision"]),
            "non_alias_failures": len(non_alias_failures),
            "first_eta_samples_all_pairs_valid": all_valid_eta,
        },
        "cases": cases,
    }


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["summary"], ensure_ascii=False, indent=2))
    print(f"result_dir: {OUT_DIR}")


if __name__ == "__main__":
    main()
