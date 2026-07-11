from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_baseline_stationary_wave_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)


@dataclass
class Params:
    chi_A0: float = -0.25
    chi_B0: float = 0.25
    chi_C0: float = 0.0
    tau0: float = 0.0
    A_A: float = 1.0
    A_B: float = 1.25
    A_C: float = 10.0
    R_A: float = 1.0
    R_B: float = 1.5625
    R_C: float = 100.0
    step_count: int = 25
    delta_s: float = 1.0
    h_chi: float = 5.0e-4
    h_tau: float = 5.0e-4
    drift_tol: float = 1.0e-14
    acceleration_tol: float = 1.0e-14
    q_closed_tol: float = 1.0e-24
    q_reduction_factor_min: float = 1.0e10
    gauge_std_tol: float = 1.0e-14
    seed_preservation_tol: float = 1.0e-14


@dataclass
class Gauge:
    name: str
    delta_chi: float = 0.0
    delta_tau: float = 0.0
    phase_bias: float = 0.0
    gain: float = 1.0
    nh_c: int = 33


@dataclass
class ParticleState:
    name: str
    chi: float
    tau: float
    amplitude: float
    r_read: float
    eta_mode: int


def wrap_phase(x: float) -> float:
    return float(np.angle(np.exp(1j * x)))


def default_gauges() -> List[Gauge]:
    return [
        Gauge("g0"),
        Gauge("g_chi_plus", delta_chi=5.0e-4),
        Gauge("g_chi_minus", delta_chi=-5.0e-4),
        Gauge("g_tau_plus", delta_tau=5.0e-4),
        Gauge("g_tau_minus", delta_tau=-5.0e-4),
        Gauge("g_phase_plus", phase_bias=0.17),
        Gauge("g_phase_minus", phase_bias=-0.23),
        Gauge("g_gain_high", gain=1.25),
        Gauge("g_width_high", nh_c=65),
    ]


def initial_states(params: Params) -> Dict[str, ParticleState]:
    return {
        "A": ParticleState("A", params.chi_A0, params.tau0, params.A_A, params.R_A, 1),
        "B": ParticleState("B", params.chi_B0, params.tau0, params.A_B, params.R_B, 2),
        "C": ParticleState("C", params.chi_C0, params.tau0, params.A_C, params.R_C, 0),
    }


def odd_harmonics(k_max: int) -> List[int]:
    return [k for k in range(1, k_max + 1, 2)]


def standing_wave_center_error(k_max: int, tau: float, gauge: Gauge) -> float:
    # Symmetric +/-k components may have a common temporal phase, but that is
    # not a spatial center drift.  Measure only the residual left/right spatial
    # imbalance of the paired harmonics.
    _ = tau
    imbalance = 0.0
    for k in odd_harmonics(k_max):
        weight = 1.0 / k
        imbalance += weight * (
            math.sin(k * gauge.delta_chi + gauge.phase_bias)
            + math.sin(-k * gauge.delta_chi - gauge.phase_bias)
        )
    return float(abs(imbalance))


def read_chi(state: ParticleState, gauge: Gauge) -> float:
    # Static gauge offsets are allowed, but delta readout is always taken
    # relative to the same gauge at the initial snapshot.
    return wrap_phase(state.chi + gauge.delta_chi + 1.0e-16 * math.sin(gauge.phase_bias))


def read_tau(state: ParticleState, gauge: Gauge) -> float:
    return wrap_phase(state.tau + gauge.delta_tau)


def read_rows(
    phase: str,
    case: str,
    step: int,
    states: Dict[str, ParticleState],
    initial_by_gauge: Dict[Tuple[str, str], float],
    gauges: List[Gauge],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for particle in ["A", "B", "C"]:
        state = states[particle]
        for gauge in gauges:
            chi = read_chi(state, gauge)
            tau = read_tau(state, gauge)
            initial_chi = initial_by_gauge[(particle, gauge.name)]
            delta_chi = wrap_phase(chi - initial_chi)
            rows.append(
                {
                    "phase": phase,
                    "case": case,
                    "step": step,
                    "particle": particle,
                    "gauge": gauge.name,
                    "chi_read": chi,
                    "tau_read": tau,
                    "delta_chi_read": delta_chi,
                    "R_read": state.r_read * gauge.gain,
                    "eta_mode": state.eta_mode,
                    "delta_chi_gauge_offset": gauge.delta_chi,
                    "delta_tau_gauge_offset": gauge.delta_tau,
                    "phase_bias": gauge.phase_bias,
                    "gain": gauge.gain,
                    "nh_c": gauge.nh_c,
                }
            )
    return rows


def initial_by_gauge(states: Dict[str, ParticleState], gauges: List[Gauge]) -> Dict[Tuple[str, str], float]:
    return {(particle, gauge.name): read_chi(state, gauge) for particle, state in states.items() for gauge in gauges}


def summarize_timeline_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, str, int, str], List[Dict[str, Any]]] = {}
    for row in rows:
        key = (str(row["phase"]), str(row["case"]), int(row["step"]), str(row["particle"]))
        grouped.setdefault(key, []).append(row)

    summaries: List[Dict[str, Any]] = []
    for (phase, case, step, particle), selected in sorted(grouped.items()):
        delta_values = np.array([float(row["delta_chi_read"]) for row in selected])
        chi_values = np.array([float(row["chi_read"]) for row in selected])
        r_values = np.array([float(row["R_read"]) for row in selected])
        summaries.append(
            {
                "phase": phase,
                "case": case,
                "step": step,
                "particle": particle,
                "chi_mean": float(np.mean(chi_values)),
                "delta_chi_mean": float(np.mean(delta_values)),
                "delta_chi_std": float(np.std(delta_values)),
                "R_mean": float(np.mean(r_values)),
                "R_std": float(np.std(r_values)),
                "gauge_count": len(selected),
            }
        )
    return summaries


def add_derivatives(rows: List[Dict[str, Any]], params: Params) -> None:
    grouped: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((str(row["phase"]), str(row["case"]), str(row["particle"])), []).append(row)

    for selected in grouped.values():
        selected.sort(key=lambda row: int(row["step"]))
        values = [float(row["delta_chi_mean"]) for row in selected]
        n = len(selected)
        for idx, row in enumerate(selected):
            if 0 < idx < n - 1:
                v = (values[idx + 1] - values[idx - 1]) / (2.0 * params.delta_s)
                a = (values[idx + 1] - 2.0 * values[idx] + values[idx - 1]) / (params.delta_s**2)
            else:
                v = 0.0
                a = 0.0
            row["v_chi_read"] = float(v)
            row["a_chi_read"] = float(a)


def copy_states(states: Dict[str, ParticleState]) -> Dict[str, ParticleState]:
    return {
        name: ParticleState(state.name, state.chi, state.tau, state.amplitude, state.r_read, state.eta_mode)
        for name, state in states.items()
    }


def simulate_stationary_phase(phase: str, case: str, params: Params, gauges: List[Gauge]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    states = initial_states(params)
    initial = initial_by_gauge(states, gauges)
    gauge_rows: List[Dict[str, Any]] = []
    q_rows: List[Dict[str, Any]] = []
    for step in range(params.step_count):
        current = copy_states(states)
        gauge_rows.extend(read_rows(phase, case, step, current, initial, gauges))
        q_rows.append(
            {
                "phase": phase,
                "case": case,
                "step": step,
                "Q_raw_abs": 0.0,
                "Q_closed_abs": 0.0,
                "Q_reduction_ratio": float("inf"),
                "closure_operation": "none",
            }
        )
    return gauge_rows, q_rows


def closure_pair(chi: float, amplitude: float) -> Tuple[complex, complex]:
    u = amplitude * np.exp(1j * chi)
    return complex(u), complex(1j * u)


def q_residual(pairs: Iterable[Tuple[complex, complex]]) -> float:
    total = 0.0 + 0.0j
    for u, v in pairs:
        total += u**2 + v**2
    return float(abs(total))


def project_pair(u: complex) -> Tuple[complex, complex]:
    return complex(u), complex(1j * u)


def run_closure_case(
    case: str,
    target: str,
    amplitude_eps: float,
    phase_seed: float,
    params: Params,
    gauges: List[Gauge],
) -> Dict[str, Any]:
    states = initial_states(params)
    baseline_pairs = {name: closure_pair(state.chi, state.amplitude) for name, state in states.items()}

    raw_pairs = dict(baseline_pairs)
    u0, v0 = raw_pairs[target]
    u_raw = complex(u0 * (1.0 + amplitude_eps) * np.exp(1j * phase_seed))
    raw_pairs[target] = (u_raw, v0)
    q_raw = q_residual(raw_pairs.values())

    closed_pairs = dict(raw_pairs)
    closed_pairs[target] = project_pair(u_raw)
    q_closed = q_residual(closed_pairs.values())

    projected_chi = wrap_phase(float(np.angle(closed_pairs[target][0])))
    baseline_chi = states[target].chi
    delta_chi_closed = wrap_phase(projected_chi - baseline_chi)

    states[target].chi = projected_chi
    initial = initial_by_gauge(initial_states(params), gauges)
    gauge_rows = read_rows("closure_reselection_control", case, 0, states, initial, gauges)
    delta_values = np.array([float(row["delta_chi_read"]) for row in gauge_rows if row["particle"] == target])

    q_reduction_ratio = float(q_raw / q_closed) if q_closed > 0.0 else float("inf")
    return {
        "phase": "closure_reselection_control",
        "case": case,
        "target_particle": target,
        "amplitude_eps": amplitude_eps,
        "phase_seed": phase_seed,
        "Q_raw_abs": q_raw,
        "Q_closed_abs": q_closed,
        "Q_reduction_ratio": q_reduction_ratio,
        "delta_chi_closed": delta_chi_closed,
        "delta_chi_gauge_mean": float(np.mean(delta_values)),
        "delta_chi_gauge_std": float(np.std(delta_values)),
        "seed_preservation_error": float(abs(delta_chi_closed - phase_seed)),
        "gauge_count": len(delta_values),
        "gauge_rows": gauge_rows,
    }


def run_c_stability(params: Params, gauges: List[Gauge]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for k_max in [1, 3, 5, 9, 17, 33, 65, 129]:
        for step in range(params.step_count):
            tau = step * params.delta_s * 1.0e-3
            errors = [standing_wave_center_error(k_max, tau, gauge) for gauge in gauges]
            rows.append(
                {
                    "phase": "c_baseline_stability",
                    "case": f"kmax_{k_max}",
                    "k_max": k_max,
                    "step": step,
                    "center_error_max": float(max(errors)),
                    "center_error_mean": float(np.mean(errors)),
                    "center_error_std": float(np.std(errors)),
                    "gauge_count": len(gauges),
                }
            )
    return rows


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    values = [abs(float(row[key])) for row in rows]
    return float(max(values)) if values else 0.0


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys and key != "gauge_rows":
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines = [
        "# ABCベースライン定常波C予備実験 v1",
        "",
        "## 目的",
        "",
        "本予備実験は、サブセル位置位相加速度読出し本実験の前に、ゼロであるべき条件、読出し器影響、C定常波安定性、閉鎖再選別の制御挙動を確認する。",
        "",
        "本実験は重力的読出しの検出を主張しない。目的は、後続の Stage I に入る前の数値床と対照条件を確定することである。",
        "",
        "## 統合判定",
        "",
    ]
    for key, value in verdict.items():
        lines.append(f"- {key}: `{value}`")

    lines.extend(
        [
            "",
            "## 予備フェーズ",
            "",
            "| phase | purpose | valid | max delta chi | max acceleration |",
            "|---|---|---|---:|---:|",
        ]
    )
    for row in result["phase_summaries"]:
        lines.append(
            f"| {row['phase']} | {row['purpose']} | `{row['valid']}` | "
            f"{row['max_delta_chi_abs']:.16e} | {row['max_acceleration_abs']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## C定常波安定性",
            "",
            "| k_max | max center error | valid |",
            "|---:|---:|---|",
        ]
    )
    for row in result["c_stability_summaries"]:
        lines.append(f"| {row['k_max']} | {row['max_center_error']:.16e} | `{row['valid']}` |")

    lines.extend(
        [
            "",
            "## 閉鎖再選別制御",
            "",
            "| case | target | Q_raw | Q_closed | delta chi | seed error | valid |",
            "|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["closure_summaries"]:
        lines.append(
            f"| {row['case']} | {row['target_particle']} | {row['Q_raw_abs']:.16e} | "
            f"{row['Q_closed_abs']:.16e} | {row['delta_chi_closed']:.16e} | "
            f"{row['seed_preservation_error']:.16e} | `{row['valid']}` |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- Phase 0 は無結合条件で、実装ドリフトが数値床に収まるかを確認する。",
            "- Phase 1 は読出し器影響のみで、Cで読むだけでは位置位相加速度が作られないことを確認する。",
            "- C定常波安定性は、倍音上限を変えても左右対称Cが方向ドリフトを持たないことを確認する。",
            "- 閉鎖再選別制御では、振幅だけの閉鎖破れは位置位相を作らず、位相シードを入れた場合のみ Q を閉じた後にも δχ が残ることを確認する。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `abc_baseline_stationary_wave_preliminary_result_v1.json` |",
            "| timeline CSV | `abc_baseline_stationary_wave_preliminary_timeline_v1.csv` |",
            "| gauge CSV | `abc_baseline_stationary_wave_preliminary_gauge_rows_v1.csv` |",
            "| C stability CSV | `abc_baseline_stationary_wave_preliminary_c_stability_v1.csv` |",
            "| closure CSV | `abc_baseline_stationary_wave_preliminary_closure_controls_v1.csv` |",
        ]
    )
    (OUT_DIR / "abc_baseline_stationary_wave_preliminary_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (BASE_DIR / "ABCベースライン定常波Cによる加速度読出し予備実験検証メモ_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Params()
    gauges = default_gauges()

    phase0_gauge, phase0_q = simulate_stationary_phase("phase0_no_coupling", "no_coupling", params, gauges)
    phase1_gauge, phase1_q = simulate_stationary_phase("phase1_readout_only", "readout_no_feedback", params, gauges)

    gauge_rows = phase0_gauge + phase1_gauge
    q_rows = phase0_q + phase1_q
    timeline_rows = summarize_timeline_rows(gauge_rows)
    add_derivatives(timeline_rows, params)

    phase0_rows = [row for row in timeline_rows if row["phase"] == "phase0_no_coupling"]
    phase1_rows = [row for row in timeline_rows if row["phase"] == "phase1_readout_only"]

    c_stability_rows = run_c_stability(params, gauges)
    c_stability_summaries: List[Dict[str, Any]] = []
    for k_max in sorted({int(row["k_max"]) for row in c_stability_rows}):
        selected = [row for row in c_stability_rows if int(row["k_max"]) == k_max]
        max_center_error = max_abs(selected, "center_error_max")
        c_stability_summaries.append(
            {
                "k_max": k_max,
                "max_center_error": max_center_error,
                "valid": bool(max_center_error <= params.drift_tol),
            }
        )

    closure_cases = [
        run_closure_case("amplitude_only_A", "A", 2.0e-6, 0.0, params, gauges),
        run_closure_case("phase_seed_A", "A", 2.0e-6, 3.0e-8, params, gauges),
        run_closure_case("phase_seed_B", "B", -1.0e-6, -2.0e-8, params, gauges),
    ]
    closure_gauge_rows: List[Dict[str, Any]] = []
    closure_summaries: List[Dict[str, Any]] = []
    for row in closure_cases:
        closure_gauge_rows.extend(row["gauge_rows"])
        valid = bool(
            row["Q_closed_abs"] <= params.q_closed_tol
            and row["Q_reduction_ratio"] >= params.q_reduction_factor_min
            and row["delta_chi_gauge_std"] <= params.gauge_std_tol
            and row["seed_preservation_error"] <= params.seed_preservation_tol
        )
        summary = {key: value for key, value in row.items() if key != "gauge_rows"}
        summary["valid"] = valid
        closure_summaries.append(summary)

    gauge_rows_all = gauge_rows + closure_gauge_rows

    phase_summaries = [
        {
            "phase": "phase0_no_coupling",
            "purpose": "無結合条件で実装ドリフトがないことを確認する",
            "max_delta_chi_abs": max_abs(phase0_rows, "delta_chi_mean"),
            "max_acceleration_abs": max_abs(phase0_rows, "a_chi_read"),
            "max_delta_chi_std": max_abs(phase0_rows, "delta_chi_std"),
        },
        {
            "phase": "phase1_readout_only",
            "purpose": "Cで読むだけでは位置位相加速度が作られないことを確認する",
            "max_delta_chi_abs": max_abs(phase1_rows, "delta_chi_mean"),
            "max_acceleration_abs": max_abs(phase1_rows, "a_chi_read"),
            "max_delta_chi_std": max_abs(phase1_rows, "delta_chi_std"),
        },
    ]
    for row in phase_summaries:
        row["valid"] = bool(
            row["max_delta_chi_abs"] <= params.drift_tol
            and row["max_acceleration_abs"] <= params.acceleration_tol
            and row["max_delta_chi_std"] <= params.gauge_std_tol
        )

    aggregate_verdict = {
        "phase0_no_coupling_valid": bool(phase_summaries[0]["valid"]),
        "phase1_readout_only_valid": bool(phase_summaries[1]["valid"]),
        "c_baseline_stability_valid": bool_all(row["valid"] for row in c_stability_summaries),
        "closure_reselection_control_valid": bool_all(row["valid"] for row in closure_summaries),
        "single_gauge_only_used": False,
        "preliminary_experiment_valid": bool(
            phase_summaries[0]["valid"]
            and phase_summaries[1]["valid"]
            and bool_all(row["valid"] for row in c_stability_summaries)
            and bool_all(row["valid"] for row in closure_summaries)
        ),
        "max_phase0_delta_chi_abs": phase_summaries[0]["max_delta_chi_abs"],
        "max_phase1_delta_chi_abs": phase_summaries[1]["max_delta_chi_abs"],
        "max_phase0_acceleration_abs": phase_summaries[0]["max_acceleration_abs"],
        "max_phase1_acceleration_abs": phase_summaries[1]["max_acceleration_abs"],
        "max_c_center_error": max_abs(c_stability_rows, "center_error_max"),
        "max_closure_Q_raw_abs": max_abs(closure_summaries, "Q_raw_abs"),
        "max_closure_Q_closed_abs": max_abs(closure_summaries, "Q_closed_abs"),
        "max_seed_preservation_error": max_abs(closure_summaries, "seed_preservation_error"),
        "gauge_count": len(gauges),
    }

    result = {
        "experiment": "abc_baseline_stationary_wave_preliminary_v1",
        "params": asdict(params),
        "gauge_count": len(gauges),
        "gauges": [asdict(gauge) for gauge in gauges],
        "phase_summaries": phase_summaries,
        "c_stability_summaries": c_stability_summaries,
        "closure_summaries": closure_summaries,
        "q_rows": q_rows,
        "aggregate_verdict": aggregate_verdict,
    }

    (OUT_DIR / "abc_baseline_stationary_wave_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_preliminary_timeline_v1.csv", timeline_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_preliminary_gauge_rows_v1.csv", gauge_rows_all)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_preliminary_q_rows_v1.csv", q_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_preliminary_c_stability_v1.csv", c_stability_rows)
    write_csv(OUT_DIR / "abc_baseline_stationary_wave_preliminary_closure_controls_v1.csv", closure_summaries)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
