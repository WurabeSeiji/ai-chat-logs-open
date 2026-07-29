#!/usr/bin/env python3
"""元論文の R=0.70 波形発展図を隔離ミラー環境で再生成する。"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_DIR = HERE.parent
ENV_DIR = TOY_DIR.parent
REPO_ROOT = ENV_DIR.parent.parent.parent
ENGINE_PATH = ENV_DIR / "20260713" / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py"
ORIGINAL_PNG = (
    REPO_ROOT
    / "波の情報読出し"
    / "20260713"
    / "exchange_scattering_matrix_fermionic_localization_transfer_preliminary_result_v1"
    / "exchange_scattering_matrix_R070_waveform_evolution_v1.png"
)
OUTPUT_PNG = HERE / "exchange_scattering_matrix_R070_waveform_evolution_reproduced_v1.png"
OUTPUT_SVG = HERE / "exchange_scattering_matrix_R070_waveform_evolution_reproduced_v1.svg"
OUTPUT_JSON = HERE / "reference_R070_reproduction_result_v1.json"
OUTPUT_REPORT = HERE / "reference_R070_reproduction_report_v1.md"

R_VALUE = 0.70
EVOLUTION_COLLISIONS = (0, 1, 2, 3, 5, 10, 20, 42)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_engine() -> Any:
    spec = importlib.util.spec_from_file_location("reference_R070_local_engine_v1", ENGINE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local engine: {ENGINE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    eng = load_engine()
    params = eng.Params()
    chi, _ = eng.make_grids(params)
    x = chi / np.pi
    snapshots = eng.recursive_snapshot_states(
        params,
        R_VALUE,
        EVOLUTION_COLLISIONS,
        n_a=1,
        n_b=63,
    )

    # 元エンジンの R070_waveform_evolution 描画ブロックと同一。
    fig, axes = eng.plt.subplots(
        4,
        2,
        figsize=(12, 12),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    metrics: list[dict[str, float | int]] = []
    for ax, collision_index in zip(axes.flatten(), EVOLUTION_COLLISIONS):
        snap = snapshots[collision_index]
        a_metrics = snap["A_channel"]
        b_metrics = snap["B_channel"]
        rho_a = snap["rho_A"] / np.max(snap["rho_A"])
        rho_b = snap["rho_B"] / np.max(snap["rho_B"])
        ax.plot(
            x,
            rho_a,
            label=f"A L={a_metrics['L']:.3g}, N={a_metrics['N_eff']:.3g}",
        )
        ax.plot(
            x,
            rho_b,
            label=f"B L={b_metrics['L']:.3g}, N={b_metrics['N_eff']:.3g}",
        )
        ax.set_title(f"R=0.70, collision={collision_index}")
        ax.set_ylabel("rho_chi / max")
        ax.legend(fontsize=7)
        metrics.append(
            {
                "collision": collision_index,
                "L_A": float(a_metrics["L"]),
                "N_eff_A": float(a_metrics["N_eff"]),
                "L_B": float(b_metrics["L"]),
                "N_eff_B": float(b_metrics["N_eff"]),
            }
        )
    for ax in axes[-1]:
        ax.set_xlabel("chi / pi")
    fig.savefig(OUTPUT_PNG, dpi=160)
    fig.savefig(OUTPUT_SVG, dpi=160)
    eng.plt.close(fig)

    reproduced_hash = sha256(OUTPUT_PNG)
    original_hash = sha256(ORIGINAL_PNG)
    byte_identical = OUTPUT_PNG.read_bytes() == ORIGINAL_PNG.read_bytes()
    reproduced_pixels = eng.plt.imread(OUTPUT_PNG)
    original_pixels = eng.plt.imread(ORIGINAL_PNG)
    same_shape = reproduced_pixels.shape == original_pixels.shape
    pixel_identical = bool(same_shape and np.array_equal(reproduced_pixels, original_pixels))
    max_pixel_abs_diff = (
        float(np.max(np.abs(reproduced_pixels.astype(float) - original_pixels.astype(float))))
        if same_shape
        else None
    )

    result = {
        "experiment": "reference_R070_reproduction_v1",
        "purpose": "reproduce the previous-paper R=0.70 waveform evolution using the isolated local mirror",
        "conditions": {
            "R": R_VALUE,
            "N_A": 1,
            "N_B": 63,
            "B_kernel": "equal-amplitude odd harmonics 1,3,...,63",
            "collisions": list(EVOLUTION_COLLISIONS),
            "display_normalization": "each rho_chi divided by its own maximum at each collision",
            "grid": [4, 2],
            "figure_size_inches": [12, 12],
            "dpi": 160,
        },
        "engine": {
            "path": str(ENGINE_PATH.relative_to(ENV_DIR)),
            "sha256": sha256(ENGINE_PATH),
        },
        "comparison_to_previous_paper_png": {
            "original_path": str(ORIGINAL_PNG.relative_to(REPO_ROOT)),
            "original_sha256": original_hash,
            "reproduced_sha256": reproduced_hash,
            "byte_identical": byte_identical,
            "same_pixel_shape": same_shape,
            "pixel_identical": pixel_identical,
            "max_pixel_abs_diff": max_pixel_abs_diff,
        },
        "metrics": metrics,
    }
    OUTPUT_JSON.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    comparison = result["comparison_to_previous_paper_png"]
    lines = [
        "# 元論文 R=0.70 波形発展図の再現結果 v1",
        "",
        "## 条件",
        "",
        "- 更新則: 隔離ミラー内の元20260713エンジン",
        "- R: 0.70",
        "- A: 基本波 N=1",
        "- B: 等振幅奇数倍音 1,3,...,63",
        "- 衝突点: 0, 1, 2, 3, 5, 10, 20, 42",
        "- 表示: 各衝突・各チャネルの rho_chi を自身の最大値で除算",
        "- 配置: 4 x 2",
        "",
        "## 元論文PNGとの比較",
        "",
        f"- バイト一致: {comparison['byte_identical']}",
        f"- 画素配列一致: {comparison['pixel_identical']}",
        f"- 最大画素絶対差: {comparison['max_pixel_abs_diff']}",
        f"- 元PNG SHA256: `{comparison['original_sha256']}`",
        f"- 再生成PNG SHA256: `{comparison['reproduced_sha256']}`",
        "",
        "## 判定",
        "",
        (
            "隔離ミラー環境で元論文の図を完全再現した。"
            if comparison["pixel_identical"]
            else "図の条件は同一だが、描画環境差を含む画素差が残った。"
        ),
        "",
    ]
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"saved: {OUTPUT_PNG}")
    print(
        "comparison:",
        f"byte_identical={byte_identical}",
        f"pixel_identical={pixel_identical}",
        f"max_pixel_abs_diff={max_pixel_abs_diff}",
    )


if __name__ == "__main__":
    main()
