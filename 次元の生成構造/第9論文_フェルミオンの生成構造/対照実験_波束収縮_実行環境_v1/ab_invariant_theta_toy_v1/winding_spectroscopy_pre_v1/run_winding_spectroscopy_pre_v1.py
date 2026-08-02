#!/usr/bin/env python3
"""巻き数分光・予備テスト v1

目的:
    公開論文の複素交換作用素 U_R の固有値構造
        lambda_a = e^{-2 pi i m / n},  n = ord(U), m = 一周期の巻き数
    に対し、「数えた整数を角度 theta/pi に運ぶ」生成則の候補として、
    B63 の関係構造（非共通倍音 31 個 x Z4 位相 = 124 スロット）の巡回
    レジスタを構成し、実際の AB 状態の反対称成分がどの巻き数 k を占有
    するかを測る（巻き数分光）。

判定:
    - 状態依存スペクトルが k = 23（または補 101 = 124-23）に集中すれば
      生成則の候補発見
    - 別の k なら反証情報
    - 順序（スロット並べ方）依存なら、写像の正準化が先決と判明

設計境界:
    - theta_from_ab・rotate_ab・状態構築は無変更で呼ぶだけ
    - 振幅探索は使用しない（等重みパケット構成のみ）
    - スロット順序の恣意性は複数の正準候補を全部走らせて順序依存性ごと報告
    - 写像自身が作るスペクトル（等重み対照）と状態寄与を分離して報告
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"

HIGH_N = 63
NONCOMMON = tuple(range(3, HIGH_N + 1, 2))  # 31 harmonics
N_SLOTS = 4 * len(NONCOMMON)  # 124
TARGET_M = 23
ORTHOGONALITY_TOLERANCE = 1.0e-10


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ab_theta_toy_for_winding_spectroscopy_v1",
        TOY_RUNNER_PATH,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load toy runner: {TOY_RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


toy = load_toy_module()
base = toy.base
plt = base.plt


def unit_norm(vector: np.ndarray) -> np.ndarray:
    norm = math.sqrt(float(np.vdot(vector, vector).real))
    if norm <= 0.0:
        raise ValueError("zero-norm state")
    return vector / norm


def make_state(packet: tuple[int, ...], which: str, source_params: Any) -> np.ndarray:
    case = base.explicit_packet_case(
        mode=f"winding_{which}_" + f"len{len(packet)}",
        packet_a=(1,),
        packet_b=packet,
    )
    return unit_norm(
        base.make_case_state(source_params, case, which, hair_enabled=True)
    )


def harmonic_amplitudes(
    field: np.ndarray,
    singles: dict[int, np.ndarray],
) -> dict[int, complex]:
    """反対称場を単一倍音基底へ射影して複素振幅 c_h を得る。"""
    return {h: complex(np.vdot(singles[h], field)) for h in singles}


ORDERINGS = {
    # slot value = c_h * (i ** q_effective)
    "O1_harm_major_phase_up": lambda j, q: (4 * j + q, q),
    "O2_phase_major_harm_up": lambda j, q: (31 * q + j, q),
    "O3_harm_desc_phase_up": lambda j, q: (4 * (30 - j) + q, q),
    "O4_harm_major_phase_down": lambda j, q: (4 * j + q, (-q) % 4),
}


def register_vector(
    coefficients: dict[int, complex],
    ordering: str,
) -> np.ndarray:
    v = np.zeros(N_SLOTS, dtype=complex)
    place = ORDERINGS[ordering]
    for j, h in enumerate(NONCOMMON):
        c = coefficients.get(h, 0.0)
        for q in range(4):
            slot, q_eff = place(j, q)
            v[slot] += c * (1j ** q_eff)
    return v


def winding_spectrum(v: np.ndarray) -> np.ndarray:
    spectrum = np.abs(np.fft.fft(v)) ** 2
    total = float(np.sum(spectrum))
    return spectrum / total if total > 0.0 else spectrum


def top_modes(spectrum: np.ndarray, count: int = 5) -> list[dict[str, float]]:
    order = np.argsort(spectrum)[::-1][:count]
    return [
        {"k": int(k), "k_complement": int((N_SLOTS - k) % N_SLOTS), "weight": float(spectrum[k])}
        for k in order
    ]


def main() -> None:
    params = base.Params(high_n=HIGH_N, recursive_collision_count=8)
    source_params = base.build_source_params(params)

    # 単一倍音基底（非共通 31 個 + 共通 1）と直交性検証
    singles = {
        h: make_state((h,), "B", source_params) for h in (1, *NONCOMMON)
    }
    keys = list(singles)
    max_cross = 0.0
    for idx, h1 in enumerate(keys):
        for h2 in keys[idx + 1 :]:
            max_cross = max(
                max_cross, abs(complex(np.vdot(singles[h1], singles[h2])))
            )
    orthogonality_ok = max_cross <= ORTHOGONALITY_TOLERANCE

    # 測定対象の状態
    a_state = make_state(tuple(range(1, HIGH_N + 1, 2)), "A", source_params)
    b_equal = make_state(tuple(range(1, HIGH_N + 1, 2)), "B", source_params)

    rng = np.random.default_rng(7)
    profiles: dict[str, dict[int, complex]] = {}

    def coeffs_from_b(b_field: np.ndarray, label: str) -> None:
        y = (a_state - b_field) / math.sqrt(2.0)
        profiles[label] = harmonic_amplitudes(y, {h: singles[h] for h in NONCOMMON})

    coeffs_from_b(b_equal, "S_B63_equal")

    # 振幅変形版（構造読出しなら巻き数は不変のはず）
    weights_inverse = {h: 1.0 / h for h in (1, *NONCOMMON)}
    b_inverse = unit_norm(
        sum(weights_inverse[h] * singles[h] for h in (1, *NONCOMMON))
    )
    coeffs_from_b(b_inverse, "S_B63_inverse_k")

    weights_random = {
        h: float(w) for h, w in zip((1, *NONCOMMON), rng.uniform(0.2, 1.0, 32))
    }
    b_random = unit_norm(
        sum(weights_random[h] * singles[h] for h in (1, *NONCOMMON))
    )
    coeffs_from_b(b_random, "S_B63_random_seed7")

    # 対照: 写像だけの寄与（c_h = 1 全等、状態情報なし）
    profiles["C_mapping_only_flat"] = {h: 1.0 + 0.0j for h in NONCOMMON}
    # 対照: 乱数複素係数（構造なし）
    profiles["C_random_complex"] = {
        h: complex(x, y)
        for h, x, y in zip(
            NONCOMMON, rng.normal(size=31), rng.normal(size=31)
        )
    }

    rows: list[dict[str, Any]] = []
    spectra_store: dict[tuple[str, str], np.ndarray] = {}
    for state_label, coeffs in profiles.items():
        for ordering in ORDERINGS:
            v = register_vector(coeffs, ordering)
            spectrum = winding_spectrum(v)
            spectra_store[(state_label, ordering)] = spectrum
            tops = top_modes(spectrum)
            rows.append(
                {
                    "state": state_label,
                    "ordering": ordering,
                    "dominant_k": tops[0]["k"],
                    "dominant_weight": tops[0]["weight"],
                    "weight_at_23": float(spectrum[TARGET_M]),
                    "weight_at_101": float(spectrum[(N_SLOTS - TARGET_M) % N_SLOTS]),
                    "top5": json.dumps(tops),
                }
            )

    # 図: 各順序の等重み状態スペクトル + 対照
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for ax, ordering in zip(axes.flatten(), ORDERINGS):
        for state_label, style in (
            ("S_B63_equal", {"color": "tab:blue", "label": "B63 equal"}),
            ("S_B63_inverse_k", {"color": "tab:green", "label": "B63 1/k"}),
            ("C_mapping_only_flat", {"color": "0.6", "label": "mapping only"}),
        ):
            spectrum = spectra_store[(state_label, ordering)]
            ax.plot(range(N_SLOTS), spectrum, linewidth=0.9, **style)
        ax.axvline(TARGET_M, color="tab:red", linestyle=":", label="k=23")
        ax.axvline(N_SLOTS - TARGET_M, color="tab:orange", linestyle=":", label="k=101")
        ax.set_title(ordering, fontsize=9)
        ax.set_xlabel("winding k")
        ax.set_ylabel("normalized |V_k|^2")
        ax.legend(fontsize=6)
        ax.grid(alpha=0.3)
    fig.suptitle(
        "Winding spectroscopy of B63 antisymmetric component on the 31x4 register",
        fontsize=12,
    )
    figure_names = []
    for ext in ("png", "svg"):
        path = HERE / f"winding_spectroscopy_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)

    csv_path = HERE / "winding_spectroscopy_pre_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "experiment": "winding_spectroscopy_pre_v1",
        "register": {
            "noncommon_harmonics": list(NONCOMMON),
            "quarter_phases": ["1", "i", "-1", "-i"],
            "n_slots": N_SLOTS,
            "target_m": TARGET_M,
        },
        "design_boundary": {
            "theta_readout_modified": False,
            "amplitude_search_used": False,
            "orderings_tested": list(ORDERINGS),
            "mapping_only_control_included": True,
        },
        "single_harmonic_basis_orthogonality": {
            "max_cross_overlap": max_cross,
            "verdict": "PASS" if orthogonality_ok else "CHECK",
        },
        "core_runner": {
            "path": TOY_RUNNER_PATH.name,
            "sha256": toy.sha256(TOY_RUNNER_PATH),
        },
        "rows": rows,
        "figures": figure_names,
    }
    (HERE / "winding_spectroscopy_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"basis orthogonality: {'PASS' if orthogonality_ok else 'CHECK'} (max {max_cross:.2e})")
    for row in rows:
        print(
            f"{row['state']:>22} | {row['ordering']:<24}"
            f" dominant k={row['dominant_k']:>3} (w={row['dominant_weight']:.3f})"
            f" w@23={row['weight_at_23']:.4f} w@101={row['weight_at_101']:.4f}"
        )


if __name__ == "__main__":
    main()
