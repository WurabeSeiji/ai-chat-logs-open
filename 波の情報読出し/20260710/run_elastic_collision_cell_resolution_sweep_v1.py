from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Dict, List


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "elastic_collision_cell_resolution_sweep_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class Params:
    A_A: float = 1.0
    A_B: float = 1.0
    A_C: float = 1000.0
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


@dataclass
class State:
    chi: float
    tau: float
    q: int


def run_case(params: Params) -> Dict[str, object]:
    eps_chi_ab = math.pi / (min(params.Nh_chi_A, params.Nh_chi_B) + 1)
    eps_tau_ab = math.pi / (min(params.Nh_tau_A, params.Nh_tau_B) + 1)
    a = State(params.chi_A0, -params.tau0, params.q_A0)
    b = State(params.chi_B0, -params.tau0, params.q_B0)

    step = 0
    collision_cell_reached = False
    crossed_without_detection = False
    min_abs_chi_gap = abs(a.chi - b.chi)
    min_abs_tau_gap = abs(a.tau - b.tau)
    previous_gap = b.chi - a.chi

    while abs(a.chi - b.chi) >= eps_chi_ab or abs(a.tau - b.tau) >= eps_tau_ab:
        if step >= params.s_max:
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        step += 1
        current_gap = b.chi - a.chi
        min_abs_chi_gap = min(min_abs_chi_gap, abs(current_gap))
        min_abs_tau_gap = min(min_abs_tau_gap, abs(a.tau - b.tau))
        if previous_gap * current_gap < 0 and abs(current_gap) >= eps_chi_ab:
            crossed_without_detection = True
            break
        previous_gap = current_gap
    else:
        collision_cell_reached = True

    collision_step = step
    if collision_cell_reached:
        a.q = -a.q
        b.q = -b.q

    post_step = 0
    post_completed = False
    while collision_cell_reached and (abs(a.chi - b.chi) <= eps_chi_ab or min(a.tau, b.tau) < params.tau0 - 1e-12):
        if post_step >= params.s_max:
            break
        a.chi += a.q * params.v_chi * params.delta_s
        b.chi += b.q * params.v_chi * params.delta_s
        a.tau += params.omega_A * params.delta_s
        b.tau += params.omega_B * params.delta_s
        post_step += 1
    else:
        post_completed = collision_cell_reached

    q_reversed = a.q == -params.q_A0 and b.q == -params.q_B0
    separated_after_collision = collision_cell_reached and a.chi < b.chi and abs(a.chi - b.chi) > eps_chi_ab
    sampling_condition = params.delta_s <= eps_chi_ab
    case_valid = collision_cell_reached and post_completed and q_reversed and separated_after_collision

    return {
        "Nh_chi_AB": min(params.Nh_chi_A, params.Nh_chi_B),
        "delta_s": params.delta_s,
        "d0": abs(params.chi_A0),
        "epsilon_chi_AB": eps_chi_ab,
        "relative_step": 2.0 * params.v_chi * params.delta_s,
        "sampling_condition_delta_s_le_epsilon": sampling_condition,
        "min_abs_chi_gap": min_abs_chi_gap,
        "min_abs_tau_gap": min_abs_tau_gap,
        "collision_cell_reached": collision_cell_reached,
        "crossed_without_detection": crossed_without_detection,
        "collision_step": collision_step,
        "post_collision_steps": post_step,
        "post_collision_propagation_completed": post_completed,
        "q_reversed": q_reversed,
        "separated_after_collision": separated_after_collision,
        "case_valid": case_valid,
    }


def write_outputs(result: Dict[str, object]) -> None:
    json_path = OUT_DIR / "cell_resolution_sweep_result_v1.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    rows = result["cases"]
    csv_path = OUT_DIR / "cell_resolution_sweep_cases_v1.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Cell Resolution Sweep Result v1",
        "",
        "## Verdict",
        "",
        f"- total_cases: `{result['summary']['total_cases']}`",
        f"- valid_cases: `{result['summary']['valid_cases']}`",
        f"- invalid_cases: `{result['summary']['invalid_cases']}`",
        f"- offgrid_valid_cases: `{result['summary']['offgrid_valid_cases']}`",
        f"- offgrid_invalid_cases: `{result['summary']['offgrid_invalid_cases']}`",
        "",
        "## Interpretation",
        "",
        "The finite collision cell is reliably detected when the calculation step is not larger than the cell width.",
        "Aligned cases can pass accidentally by landing exactly on the center; off-grid cases expose skipped-cell failures.",
        "",
        "## Cases",
        "",
        "| d0 | Nh_chi_AB | delta_s | epsilon_chi_AB | min_abs_chi_gap | collision_cell_reached | crossed_without_detection | case_valid |",
        "|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['d0']:.12g} | {row['Nh_chi_AB']} | {row['delta_s']:.12g} | "
            f"{row['epsilon_chi_AB']:.12g} | {row['min_abs_chi_gap']:.12g} | "
            f"{row['collision_cell_reached']} | {row['crossed_without_detection']} | {row['case_valid']} |"
        )
    (OUT_DIR / "cell_resolution_sweep_report_v1.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    for d0 in sorted({row["d0"] for row in rows}):
        subset = [row for row in rows if row["d0"] == d0]
        nhs = sorted({row["Nh_chi_AB"] for row in subset})
        steps = sorted({row["delta_s"] for row in subset})
        matrix = []
        for nh in nhs:
            matrix_row = []
            for delta_s in steps:
                case = next(row for row in subset if row["Nh_chi_AB"] == nh and row["delta_s"] == delta_s)
                matrix_row.append(1 if case["case_valid"] else 0)
            matrix.append(matrix_row)

        fig, ax = plt.subplots(figsize=(8, 5))
        image = ax.imshow(matrix, aspect="auto", origin="lower", cmap="RdYlGn", vmin=0, vmax=1)
        ax.set_xticks(range(len(steps)))
        ax.set_xticklabels([f"{value:g}" for value in steps], rotation=45)
        ax.set_yticks(range(len(nhs)))
        ax.set_yticklabels([str(value) for value in nhs])
        ax.set_xlabel("delta_s")
        ax.set_ylabel("Nh_chi_AB")
        ax.set_title(f"Collision cell detection validity, d0={d0:g}")
        fig.colorbar(image, ax=ax, ticks=[0, 1], label="case_valid")
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"cell_resolution_validity_d0_{str(d0).replace('.', '_')}_v1.png", dpi=180)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for d0 in sorted({row["d0"] for row in rows}):
        subset = [row for row in rows if row["d0"] == d0]
        ax.scatter(
            [row["delta_s"] / row["epsilon_chi_AB"] for row in subset],
            [row["min_abs_chi_gap"] / row["epsilon_chi_AB"] for row in subset],
            c=["green" if row["case_valid"] else "red" for row in subset],
            label=f"d0={d0:g}",
            alpha=0.75,
        )
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1, label="delta_s = epsilon")
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=1, label="min gap = epsilon")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("delta_s / epsilon_chi_AB")
    ax.set_ylabel("min_abs_chi_gap / epsilon_chi_AB")
    ax.set_title("Cell sampling condition")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "cell_resolution_sampling_condition_v1.png", dpi=180)
    plt.close(fig)


def run() -> Dict[str, object]:
    base = Params()
    nh_values = [19, 39, 99, 199, 399, 999]
    delta_values = [0.02, 0.01, 0.005, 0.002, 0.001]
    d0_values = [0.2, 0.203]
    cases: List[Dict[str, object]] = []
    for d0 in d0_values:
        for nh in nh_values:
            for delta_s in delta_values:
                params = replace(
                    base,
                    Nh_chi_A=nh,
                    Nh_chi_B=nh,
                    Nh_tau_A=nh,
                    Nh_tau_B=nh,
                    chi_A0=-d0,
                    chi_B0=d0,
                    delta_s=delta_s,
                )
                cases.append(run_case(params))

    offgrid = [row for row in cases if not math.isclose(row["d0"], 0.2)]
    result = {
        "parameters": asdict(base),
        "summary": {
            "total_cases": len(cases),
            "valid_cases": sum(1 for row in cases if row["case_valid"]),
            "invalid_cases": sum(1 for row in cases if not row["case_valid"]),
            "offgrid_valid_cases": sum(1 for row in offgrid if row["case_valid"]),
            "offgrid_invalid_cases": sum(1 for row in offgrid if not row["case_valid"]),
        },
        "cases": cases,
    }
    return result


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["summary"], ensure_ascii=False, indent=2))
    print(f"result_dir: {OUT_DIR}")


if __name__ == "__main__":
    main()
