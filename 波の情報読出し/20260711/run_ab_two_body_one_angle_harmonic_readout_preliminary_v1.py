from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "ab_two_body_one_angle_harmonic_readout_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TAU = 2.0 * math.pi


@dataclass
class Params:
    step_count: int = 720
    omega_step: float = TAU / 96.0
    rho_AB: float = 1.0
    closure_tol: float = 1.0e-12
    protocol_tol: float = 1.0e-12
    f_consistency_tol: float = 3.0e-6
    decay_order_tol: float = 2.0e-7


@dataclass(frozen=True)
class InitialCase:
    case_id: str
    deviation_deg: float


@dataclass(frozen=True)
class ReadoutMode:
    name: str
    per_step_leak: float
    active_readout: bool


INITIAL_CASES = [
    InitialCase("near_pi_02deg", 2.0),
    InitialCase("near_pi_05deg", 5.0),
    InitialCase("near_pi_10deg", 10.0),
    InitialCase("near_pi_20deg", 20.0),
]

PROTOCOLS = ["Protocol_F", "Protocol_B"]

READOUT_MODES = [
    ReadoutMode("readout_off", 0.0, False),
    ReadoutMode("readout_weak", 1.0e-5, True),
    ReadoutMode("readout_normal", 5.0e-5, True),
    ReadoutMode("readout_strong", 2.0e-4, True),
]


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def min_max(values: List[float]) -> Tuple[float, float]:
    if not values:
        return 0.0, 0.0
    return float(min(values)), float(max(values))


def closure_rotation_series(initial_deviation_rad: float, mode: ReadoutMode, params: Params) -> List[complex]:
    z = complex(initial_deviation_rad, 0.0)
    lam = 1.0 - mode.per_step_leak
    rot = complex(math.cos(params.omega_step), math.sin(params.omega_step))
    series: List[complex] = []
    for _ in range(params.step_count + 1):
        series.append(z)
        z = lam * rot * z
    return series


def label_free_readout(signed_deviation: float) -> Tuple[float, float, float]:
    deviation = min(abs(signed_deviation), math.pi)
    d_near = math.pi - deviation
    d_far = math.pi + deviation
    v_ab = deviation * deviation
    return d_near, d_far, v_ab


def protocol_display_deviation(protocol: str, signed_deviation: float) -> float:
    if protocol == "Protocol_F":
        return abs(signed_deviation)
    if protocol == "Protocol_B":
        return signed_deviation
    raise ValueError(f"unknown protocol: {protocol}")


def estimate_decay_rate(step_values: List[float], envelope_values: List[float]) -> float:
    x = np.array(step_values, dtype=float)
    y = np.array(envelope_values, dtype=float)
    mask = y > 1.0e-30
    x = x[mask]
    y = y[mask]
    if len(x) < 3:
        return 0.0
    slope, _ = np.polyfit(x, np.log(y), 1)
    return float(slope)


def count_internal_sign_changes(series: List[complex]) -> int:
    signs: List[int] = []
    for z in series:
        value = float(z.real)
        if abs(value) <= 1.0e-14:
            continue
        signs.append(1 if value > 0.0 else -1)
    return sum(1 for prev, cur in zip(signs, signs[1:]) if prev != cur)


def rows_for_case(case: InitialCase, protocol: str, mode: ReadoutMode, params: Params) -> List[Dict[str, Any]]:
    initial_deviation_rad = math.radians(case.deviation_deg)
    series = closure_rotation_series(initial_deviation_rad, mode, params)
    omega_discrete_sq = 4.0 * math.sin(params.omega_step / 2.0) ** 2
    rows: List[Dict[str, Any]] = []
    theta_values = np.unwrap([math.atan2(z.imag, z.real) for z in series])
    for step, z in enumerate(series):
        signed_deviation = float(z.real)
        display_deviation = protocol_display_deviation(protocol, signed_deviation)
        d_near, d_far, v_ab = label_free_readout(signed_deviation)
        envelope_abs = abs(z)
        envelope_v = envelope_abs * envelope_abs
        f_center = omega_discrete_sq * abs(signed_deviation)
        if 0 < step < len(series) - 1:
            second_diff = float(series[step + 1].real - 2.0 * series[step].real + series[step - 1].real)
            f_circle = abs(second_diff)
        else:
            second_diff = 0.0
            f_circle = f_center
        q_raw = mode.per_step_leak * envelope_v
        q_closed = 0.0
        rows.append(
            {
                "case_id": case.case_id,
                "initial_deviation_deg": case.deviation_deg,
                "protocol": protocol,
                "readout_mode": mode.name,
                "active_readout": mode.active_readout,
                "step": step,
                "internal_signed_deviation_rad": signed_deviation,
                "protocol_display_deviation_rad": display_deviation,
                "D_AB_near_rad": d_near,
                "D_AB_far_rad": d_far,
                "D_AB_near_deg": math.degrees(d_near),
                "D_AB_far_deg": math.degrees(d_far),
                "V_AB": v_ab,
                "rho_AB": params.rho_AB,
                "theta_AB_unwrapped": float(theta_values[step]),
                "closure_complement": float(z.imag),
                "envelope_AB_abs": envelope_abs,
                "envelope_V_AB": envelope_v,
                "f_AB_center": f_center,
                "f_AB_circle": f_circle,
                "f_AB_projection_consistency_error": abs(f_center - f_circle),
                "second_difference_internal_deviation": second_diff,
                "Q_raw": q_raw,
                "Q_closed": q_closed,
                "closure_relaxation": q_raw - q_closed,
                "absolute_background_axis_used": False,
                "f_A_or_f_B_used": False,
            }
        )
    return rows


def summarize_case(rows: List[Dict[str, Any]], params: Params) -> Dict[str, Any]:
    step_values = [float(row["step"]) for row in rows]
    envelope_values = [float(row["envelope_V_AB"]) for row in rows]
    v_values = [float(row["V_AB"]) for row in rows]
    near_values = [float(row["D_AB_near_rad"]) for row in rows]
    f_errors = [float(row["f_AB_projection_consistency_error"]) for row in rows[1:-1]]
    q_closed_values = [abs(float(row["Q_closed"])) for row in rows]
    q_raw_values = [abs(float(row["Q_raw"])) for row in rows]
    internal_series = [complex(float(row["internal_signed_deviation_rad"]), float(row["closure_complement"])) for row in rows]
    v_min, v_max = min_max(v_values)
    near_min, near_max = min_max(near_values)
    return {
        "case_id": rows[0]["case_id"],
        "initial_deviation_deg": rows[0]["initial_deviation_deg"],
        "protocol": rows[0]["protocol"],
        "readout_mode": rows[0]["readout_mode"],
        "active_readout": rows[0]["active_readout"],
        "step_count": len(rows) - 1,
        "D_AB_near_min_deg": math.degrees(near_min),
        "D_AB_near_max_deg": math.degrees(near_max),
        "V_AB_min": v_min,
        "V_AB_max": v_max,
        "envelope_V_AB_initial": envelope_values[0],
        "envelope_V_AB_final": envelope_values[-1],
        "envelope_ratio_final_over_initial": envelope_values[-1] / envelope_values[0] if envelope_values[0] else 0.0,
        "decay_rate_V_AB": estimate_decay_rate(step_values, envelope_values),
        "internal_sign_change_count": count_internal_sign_changes(internal_series),
        "oscillation_detected": count_internal_sign_changes(internal_series) >= 4,
        "max_Q_raw_abs": max(q_raw_values) if q_raw_values else 0.0,
        "max_Q_closed_abs": max(q_closed_values) if q_closed_values else 0.0,
        "max_f_AB_projection_consistency_error": max(f_errors) if f_errors else 0.0,
        "absolute_background_axis_used": False,
        "f_A_or_f_B_used": False,
    }


def compare_protocols(summaries: List[Dict[str, Any]], series_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows_by_key: Dict[Tuple[str, str, str], List[Dict[str, Any]]] = {}
    for row in series_rows:
        key = (str(row["case_id"]), str(row["readout_mode"]), str(row["protocol"]))
        rows_by_key.setdefault(key, []).append(row)

    comparison_rows: List[Dict[str, Any]] = []
    for case in INITIAL_CASES:
        for mode in READOUT_MODES:
            f_rows = rows_by_key[(case.case_id, mode.name, "Protocol_F")]
            b_rows = rows_by_key[(case.case_id, mode.name, "Protocol_B")]
            max_d_near_diff = max(
                abs(float(f_row["D_AB_near_rad"]) - float(b_row["D_AB_near_rad"]))
                for f_row, b_row in zip(f_rows, b_rows)
            )
            max_v_diff = max(abs(float(f_row["V_AB"]) - float(b_row["V_AB"])) for f_row, b_row in zip(f_rows, b_rows))
            max_display_diff = max(
                abs(float(f_row["protocol_display_deviation_rad"]) - float(b_row["protocol_display_deviation_rad"]))
                for f_row, b_row in zip(f_rows, b_rows)
            )
            comparison_rows.append(
                {
                    "case_id": case.case_id,
                    "readout_mode": mode.name,
                    "max_D_AB_near_protocol_diff": max_d_near_diff,
                    "max_V_AB_protocol_diff": max_v_diff,
                    "max_protocol_display_deviation_diff": max_display_diff,
                    "label_free_protocol_degenerate": bool(max_d_near_diff <= 1.0e-15 and max_v_diff <= 1.0e-15),
                }
            )
    return comparison_rows


def readout_decay_rows(summaries: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_key = {(row["case_id"], row["protocol"], row["readout_mode"]): row for row in summaries}
    rows: List[Dict[str, Any]] = []
    for case in INITIAL_CASES:
        for protocol in PROTOCOLS:
            off = by_key[(case.case_id, protocol, "readout_off")]
            weak = by_key[(case.case_id, protocol, "readout_weak")]
            normal = by_key[(case.case_id, protocol, "readout_normal")]
            strong = by_key[(case.case_id, protocol, "readout_strong")]
            rates = [
                abs(float(off["decay_rate_V_AB"])),
                abs(float(weak["decay_rate_V_AB"])),
                abs(float(normal["decay_rate_V_AB"])),
                abs(float(strong["decay_rate_V_AB"])),
            ]
            rows.append(
                {
                    "case_id": case.case_id,
                    "protocol": protocol,
                    "decay_rate_off": off["decay_rate_V_AB"],
                    "decay_rate_weak": weak["decay_rate_V_AB"],
                    "decay_rate_normal": normal["decay_rate_V_AB"],
                    "decay_rate_strong": strong["decay_rate_V_AB"],
                    "abs_decay_monotonic_with_readout_strength": bool(
                        rates[0] <= rates[1] + 1.0e-12
                        and rates[1] <= rates[2] + 1.0e-12
                        and rates[2] <= rates[3] + 1.0e-12
                    ),
                    "strong_over_off_abs_decay_ratio": rates[3] / max(rates[0], 1.0e-30),
                    "strong_over_normal_abs_decay_ratio": rates[3] / max(rates[2], 1.0e-30),
                }
            )
    return rows


def make_plots(series_rows: List[Dict[str, Any]], summaries: List[Dict[str, Any]], decay_rows: List[Dict[str, Any]]) -> None:
    plot_case = "near_pi_05deg"
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for mode in READOUT_MODES:
        rows = [
            row
            for row in series_rows
            if row["case_id"] == plot_case and row["protocol"] == "Protocol_B" and row["readout_mode"] == mode.name
        ]
        steps = [int(row["step"]) for row in rows]
        v_values = [float(row["V_AB"]) for row in rows]
        envelope = [float(row["envelope_V_AB"]) for row in rows]
        axes[0].plot(steps, v_values, label=mode.name, linewidth=1.2)
        axes[1].plot(steps, envelope, label=mode.name, linewidth=1.2)
    axes[0].set_ylabel("V_AB")
    axes[1].set_ylabel("envelope_V_AB")
    axes[1].set_xlabel("step")
    axes[0].set_title("AB one-angle harmonic readout: readout-mode comparison")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_one_angle_harmonic_readout_readout_mode_comparison_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for protocol in PROTOCOLS:
        rows = [
            row
            for row in series_rows
            if row["case_id"] == plot_case and row["protocol"] == protocol and row["readout_mode"] == "readout_off"
        ]
        steps = [int(row["step"]) for row in rows]
        v_values = [float(row["V_AB"]) for row in rows]
        ax.plot(steps, v_values, label=protocol, linewidth=1.3)
    ax.set_title("Protocol F/B label-free V_AB comparison")
    ax.set_xlabel("step")
    ax.set_ylabel("V_AB")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_one_angle_harmonic_readout_protocol_comparison_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    labels = [f"{row['case_id']}\n{row['protocol']}" for row in decay_rows]
    x = np.arange(len(labels))
    ax.plot(x, [abs(float(row["decay_rate_off"])) for row in decay_rows], label="off", marker="o")
    ax.plot(x, [abs(float(row["decay_rate_weak"])) for row in decay_rows], label="weak", marker="o")
    ax.plot(x, [abs(float(row["decay_rate_normal"])) for row in decay_rows], label="normal", marker="o")
    ax.plot(x, [abs(float(row["decay_rate_strong"])) for row in decay_rows], label="strong", marker="o")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("|decay_rate_V_AB|")
    ax.set_title("Envelope decay by readout mode")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_one_angle_harmonic_readout_envelope_decay_v1.png", dpi=180)
    plt.close(fig)


def aggregate_verdict(
    params: Params,
    summaries: List[Dict[str, Any]],
    protocol_rows: List[Dict[str, Any]],
    decay_rows: List[Dict[str, Any]],
) -> Dict[str, Any]:
    max_q_closed = max(float(row["max_Q_closed_abs"]) for row in summaries)
    max_f_error = max(float(row["max_f_AB_projection_consistency_error"]) for row in summaries)
    nonstrong_summaries = [row for row in summaries if row["readout_mode"] != "readout_strong"]
    strong_summaries = [row for row in summaries if row["readout_mode"] == "readout_strong"]
    max_f_error_nonstrong = max(float(row["max_f_AB_projection_consistency_error"]) for row in nonstrong_summaries)
    max_f_error_strong = max(float(row["max_f_AB_projection_consistency_error"]) for row in strong_summaries)
    max_protocol_d = max(float(row["max_D_AB_near_protocol_diff"]) for row in protocol_rows)
    max_protocol_v = max(float(row["max_V_AB_protocol_diff"]) for row in protocol_rows)
    f_consistent_nonstrong = bool(max_f_error_nonstrong <= params.f_consistency_tol)
    strong_readout_perturbs_f_projection = bool(max_f_error_strong > params.f_consistency_tol)
    return {
        "case_count": len(summaries),
        "initial_case_count": len(INITIAL_CASES),
        "protocol_count": len(PROTOCOLS),
        "readout_mode_count": len(READOUT_MODES),
        "observer_C_used": False,
        "single_gauge_only_used": False,
        "absolute_background_axis_used": False,
        "f_A_or_f_B_used": False,
        "standard_force_law_used": False,
        "max_Q_closed_abs": max_q_closed,
        "max_Q_raw_abs": max(float(row["max_Q_raw_abs"]) for row in summaries),
        "max_f_AB_projection_consistency_error": max_f_error,
        "max_f_AB_projection_consistency_error_nonstrong": max_f_error_nonstrong,
        "max_f_AB_projection_consistency_error_strong": max_f_error_strong,
        "f_AB_projection_consistent_nonstrong_modes": f_consistent_nonstrong,
        "strong_readout_perturbs_f_AB_projection": strong_readout_perturbs_f_projection,
        "max_D_AB_near_protocol_diff": max_protocol_d,
        "max_V_AB_protocol_diff": max_protocol_v,
        "oscillation_detected_all_cases": bool_all(bool(row["oscillation_detected"]) for row in summaries),
        "label_free_protocol_degenerate_all_cases": bool_all(
            bool(row["label_free_protocol_degenerate"]) for row in protocol_rows
        ),
        "readout_decay_monotonic_all_cases": bool_all(
            bool(row["abs_decay_monotonic_with_readout_strength"]) for row in decay_rows
        ),
        "readout_off_decay_max_abs": max(
            abs(float(row["decay_rate_off"])) for row in decay_rows
        ),
        "readout_strong_decay_min_abs": min(
            abs(float(row["decay_rate_strong"])) for row in decay_rows
        ),
        "ab_one_angle_harmonic_readout_preliminary_valid": bool(
            max_q_closed <= params.closure_tol
            and max_protocol_d <= params.protocol_tol
            and max_protocol_v <= params.protocol_tol
            and f_consistent_nonstrong
            and bool_all(bool(row["oscillation_detected"]) for row in summaries)
            and bool_all(bool(row["label_free_protocol_degenerate"]) for row in protocol_rows)
            and bool_all(bool(row["abs_decay_monotonic_with_readout_strength"]) for row in decay_rows)
        ),
    }


def write_report(result: Dict[str, Any]) -> None:
    verdict = result["aggregate_verdict"]
    lines: List[str] = [
        "# AB二体閉鎖位相系における一角度円周位相調和読出し予備実験検証メモ v1",
        "",
        "## 目的",
        "",
        "観測機 C を置かない AB 二体系で、ラベルなし二弧相対位相 `D_AB`、対称偏差 `V_AB`、AB 合成補償 `f_AB` を読む予備実験を行った。",
        "",
        "本実験では、標準重力式、標準クーロン式、標準ばね式を使わない。",
        "",
        "一角度の閉鎖補助平面における複素回転",
        "",
        "```text",
        "z(s+1) = lambda * exp(i Omega) * z(s)",
        "```",
        "",
        "を用い、`Protocol F/B` と読出し波条件の違いがラベルなし読出しへどう現れるかを検査した。",
        "",
        "## 統合判定",
        "",
        "| 量 | 値 |",
        "|---|---:|",
    ]
    for key, value in verdict.items():
        lines.append(f"| `{key}` | `{value}` |")

    lines.extend(
        [
            "",
            "## ケース別サマリー",
            "",
            "| case | protocol | readout | decay_rate | envelope final/initial | sign changes | max f error |",
            "|---|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in result["case_summaries"]:
        lines.append(
            f"| {row['case_id']} | {row['protocol']} | {row['readout_mode']} | "
            f"{row['decay_rate_V_AB']:.16e} | {row['envelope_ratio_final_over_initial']:.16e} | "
            f"{row['internal_sign_change_count']} | {row['max_f_AB_projection_consistency_error']:.16e} |"
        )

    lines.extend(
        [
            "",
            "## Protocol F/B 比較",
            "",
            "| case | readout | max D diff | max V diff | display diff | label-free degenerate |",
            "|---|---|---:|---:|---:|---|",
        ]
    )
    for row in result["protocol_comparison"]:
        lines.append(
            f"| {row['case_id']} | {row['readout_mode']} | "
            f"{row['max_D_AB_near_protocol_diff']:.16e} | {row['max_V_AB_protocol_diff']:.16e} | "
            f"{row['max_protocol_display_deviation_diff']:.16e} | {row['label_free_protocol_degenerate']} |"
        )

    lines.extend(
        [
            "",
            "## 読出し波停止反証テスト",
            "",
            "| case | protocol | off | weak | normal | strong | monotonic |",
            "|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["readout_decay_comparison"]:
        lines.append(
            f"| {row['case_id']} | {row['protocol']} | "
            f"{row['decay_rate_off']:.16e} | {row['decay_rate_weak']:.16e} | "
            f"{row['decay_rate_normal']:.16e} | {row['decay_rate_strong']:.16e} | "
            f"{row['abs_decay_monotonic_with_readout_strength']} |"
        )

    lines.extend(
        [
            "",
            "## 解釈",
            "",
            "- `Protocol F/B` は内部表示としては異なるが、`D_AB` と `V_AB` では縮退した。",
            "- `readout_off` では包絡減衰が数値丸め範囲に留まり、読出し波を強くするほど減衰率が大きくなった。",
            "- これは、読出し波が長期振幅へ影響しうるという反証テストの検出系として機能する。",
            "- ただし、この予備実験は複素回転写像の検査であり、調和読出しが第一原理から自発的に出現したことの証明ではない。",
            "- 逆二乗型は本実験の対象外であり、二角度以上の位置位相自由度拡張で検査する。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            "| JSON | `ab_two_body_one_angle_harmonic_readout_preliminary_result_v1.json` |",
            "| series CSV | `ab_two_body_one_angle_harmonic_readout_series_v1.csv` |",
            "| case summary CSV | `ab_two_body_one_angle_harmonic_readout_case_summary_v1.csv` |",
            "| protocol comparison CSV | `ab_two_body_one_angle_harmonic_readout_protocol_comparison_v1.csv` |",
            "| readout decay CSV | `ab_two_body_one_angle_harmonic_readout_readout_decay_v1.csv` |",
            "| readout mode plot | `ab_two_body_one_angle_harmonic_readout_readout_mode_comparison_v1.png` |",
            "| protocol plot | `ab_two_body_one_angle_harmonic_readout_protocol_comparison_v1.png` |",
            "| decay plot | `ab_two_body_one_angle_harmonic_readout_envelope_decay_v1.png` |",
        ]
    )
    report = "\n".join(lines) + "\n"
    (OUT_DIR / "ab_two_body_one_angle_harmonic_readout_preliminary_report_v1.md").write_text(
        report, encoding="utf-8"
    )
    (BASE_DIR / "AB二体閉鎖位相系における一角度円周位相調和読出し予備実験検証メモ_v1.md").write_text(
        report, encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    params = Params()
    all_rows: List[Dict[str, Any]] = []
    summaries: List[Dict[str, Any]] = []
    for case in INITIAL_CASES:
        for protocol in PROTOCOLS:
            for mode in READOUT_MODES:
                rows = rows_for_case(case, protocol, mode, params)
                all_rows.extend(rows)
                summaries.append(summarize_case(rows, params))

    protocol_rows = compare_protocols(summaries, all_rows)
    decay_rows = readout_decay_rows(summaries)
    verdict = aggregate_verdict(params, summaries, protocol_rows, decay_rows)
    result = {
        "experiment": "ab_two_body_one_angle_harmonic_readout_preliminary_v1",
        "params": asdict(params),
        "initial_cases": [asdict(case) for case in INITIAL_CASES],
        "readout_modes": [asdict(mode) for mode in READOUT_MODES],
        "protocols": PROTOCOLS,
        "case_summaries": summaries,
        "protocol_comparison": protocol_rows,
        "readout_decay_comparison": decay_rows,
        "aggregate_verdict": verdict,
        "note": (
            "This is a preliminary closure-rotation readout test. It does not use F=Gmm/r^2, F=kx, "
            "or signed observer-target labels. Protocol F/B differences are tested only through label-free readouts."
        ),
    }
    (OUT_DIR / "ab_two_body_one_angle_harmonic_readout_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "ab_two_body_one_angle_harmonic_readout_series_v1.csv", all_rows)
    write_csv(OUT_DIR / "ab_two_body_one_angle_harmonic_readout_case_summary_v1.csv", summaries)
    write_csv(OUT_DIR / "ab_two_body_one_angle_harmonic_readout_protocol_comparison_v1.csv", protocol_rows)
    write_csv(OUT_DIR / "ab_two_body_one_angle_harmonic_readout_readout_decay_v1.csv", decay_rows)
    make_plots(all_rows, summaries, decay_rows)
    write_report(result)
    return result


def main() -> None:
    result = run()
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
