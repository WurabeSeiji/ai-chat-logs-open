#!/usr/bin/env python3
"""系統A 局在性交換Rスイープ 計測版ランナー v1（instrumented）

目的（波束収縮と交差項の関係性調査.md §3.2/§5 への対応）:
  既存実験では保存されていなかった「倍音別複素係数」を、物理・既存出力を一切変えずに
  追加取得する。交差相関行列 C_mn = a_m a_n^* はこの係数から後段でオフライン再構成できる。

設計:
  - 基底ランナー run_system_A_localization_exchange_R_sweep_preliminary_v1.py を
    無修正のまま importlib でロードし、run_case のみ計測版に差し替える。
  - 計測版 run_case は基底版と同一の状態初期化・同一の更新式・同一の行生成を実行し、
    加えて各衝突時点の chi 方向 FFT 複素係数（符号付き倍音 n ∈ [-M, M]、eta 次元込み）を
    記録して npz に保存する。既存 CSV 出力の数値経路には一切触れない。
  - 保存先は出力ディレクトリ配下の harmonic_dump_v1/。1 (case, R) につき 1 npz。

npz 内容:
  coeffs        complex128 (n_records, 2, n_harmonics, eta_grid_n)
                チャネル軸は [A_channel, B_channel]。FFT は np.fft.fft(axis=0, norm="ortho")。
                状態ベクトルは基底実装どおり各衝突後に normalize 済み（ノルム 1）。
  collisions    int64 (n_records,)   記録した衝突番号（--dump-stride 適用後）
  harmonics     int64 (n_harmonics,) 符号付き倍音番号 [-M..M]
  coverage      float64 (n_records, 2) 捕捉パワー / 全パワー（切り捨て監査用、通常 ≈1）
  meta          JSON 文字列: case_id, packet 仕様, R_input, R, T, t, r, 格子サイズ, 規約

使い方: 基底ランナーと同じ引数に加えて
  --dump-max-n N      符号付き倍音の記録範囲 M（既定: metrics.max_n = min(chi_n/2, high_n+2)）
  --dump-stride K     衝突記録間引き（既定 1 = 全衝突）
  --no-dump           ダンプ無効（純粋な基底互換動作）
  ※ --refine-minima は本計測版では未対応（基底ランナーを使うこと）
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py"


def load_base_module() -> Any:
    spec = importlib.util.spec_from_file_location("system_a_r_sweep_base_v1", BASE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load base runner: {BASE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = load_base_module()

# --- 計測設定（main で CLI から設定） ---
DUMP_ENABLED = True
DUMP_MAX_N: int | None = None
DUMP_STRIDE = 1
DUMP_MANIFEST: List[Dict[str, Any]] = []


def signed_harmonic_indices(metrics: Any, max_abs_n: int) -> tuple[np.ndarray, np.ndarray]:
    """符号付き倍音 n ∈ [-M, M] に対応する chi-FFT 行インデックスを返す。"""
    rounded = np.rint(metrics.freqs).astype(int)
    harmonics = np.arange(-max_abs_n, max_abs_n + 1, dtype=np.int64)
    index_of = {int(f): i for i, f in enumerate(rounded)}
    rows = np.array([index_of[int(n)] for n in harmonics], dtype=np.int64)
    return harmonics, rows


def spectral_coeffs(source_params: Any, vector: np.ndarray, rows: np.ndarray) -> tuple[np.ndarray, float]:
    """chi 方向 FFT（norm='ortho'）の複素係数行（eta 次元込み）と捕捉率を返す。"""
    arr = vector.reshape(source_params.chi_grid_n, source_params.eta_grid_n)
    transformed = np.fft.fft(arr, axis=0, norm="ortho")
    coeffs = transformed[rows, :]
    total = float(np.vdot(transformed, transformed).real)
    captured = float(np.vdot(coeffs, coeffs).real)
    coverage = captured / total if total > 0.0 else float("nan")
    return coeffs, coverage


def instrumented_run_case(
    source_params: Any,
    metrics: Any,
    case: Any,
    r_value: float,
    max_collision: int,
) -> List[Dict[str, Any]]:
    """基底 run_case と同一の物理経路＋倍音別複素係数の追加記録。

    状態初期化・更新式・行生成は基底実装の呼び出しをそのまま用いる。
    （基底 run_case: delta→(t,r,T,R)→make_case_state→衝突ループ
      a' = normalize(r*a + t*b), b' = normalize(t*a + r*b)）
    """
    delta_f = base.src.delta_from_reflection_rate(r_value)
    t, r, T, R = base.src.scattering_coefficients(delta_f)
    hair_enabled = True
    a = base.make_case_state(source_params, case, "A", hair_enabled)
    b = base.make_case_state(source_params, case, "B", hair_enabled)
    initial_a = a.copy()
    initial_b = b.copy()
    h_a0 = metrics.harmonic_distribution(initial_a)
    h_b0 = metrics.harmonic_distribution(initial_b)

    max_abs_n = DUMP_MAX_N if DUMP_MAX_N is not None else int(metrics.max_n)
    harmonics, rows_idx = signed_harmonic_indices(metrics, max_abs_n)
    dump_collisions: List[int] = []
    dump_coeffs: List[np.ndarray] = []
    dump_coverage: List[List[float]] = []

    rows: List[Dict[str, Any]] = []
    for collision in range(max_collision + 1):
        rows.append(base.row_for_state(source_params, metrics, case, r_value, t, r, T, R, collision, "A_channel", a, h_a0, h_b0, initial_a, initial_b))
        rows.append(base.row_for_state(source_params, metrics, case, r_value, t, r, T, R, collision, "B_channel", b, h_a0, h_b0, initial_a, initial_b))
        if DUMP_ENABLED and (collision % DUMP_STRIDE == 0 or collision == max_collision):
            ca, cov_a = spectral_coeffs(source_params, a, rows_idx)
            cb, cov_b = spectral_coeffs(source_params, b, rows_idx)
            dump_collisions.append(collision)
            dump_coeffs.append(np.stack([ca, cb]))
            dump_coverage.append([cov_a, cov_b])
        if collision >= max_collision:
            break
        a_next = base.src.normalize(r * a + t * b)
        b_next = base.src.normalize(t * a + r * b)
        a, b = a_next, b_next

    if DUMP_ENABLED:
        dump_dir = base.OUT_DIR / "harmonic_dump_v1"
        dump_dir.mkdir(parents=True, exist_ok=True)
        stem = f"harmonic_coeffs_{base.safe_slug(case.case_id, max_len=60)}_R{base.compact_float(r_value)}_v1"
        path = dump_dir / f"{stem}.npz"
        meta = {
            "case_id": case.case_id,
            "mode": case.mode,
            "packet_a": list(case.packet_a),
            "packet_b": list(case.packet_b),
            "packet_a_text": case.packet_a_text,
            "packet_b_text": case.packet_b_text,
            "R_input": float(r_value),
            "R": float(R),
            "T": float(T),
            "t": [float(t.real), float(t.imag)],
            "r": [float(r.real), float(r.imag)],
            "delta_f": float(delta_f),
            "max_collision": int(max_collision),
            "dump_stride": int(DUMP_STRIDE),
            "dump_max_n": int(max_abs_n),
            "chi_grid_n": int(source_params.chi_grid_n),
            "eta_grid_n": int(source_params.eta_grid_n),
            "fft_convention": "np.fft.fft(axis=0, norm='ortho') on reshape(chi, eta); states normalized to unit norm",
            "channel_order": ["A_channel", "B_channel"],
            "base_runner": BASE_PATH.name,
            "source_engine": str(base.SOURCE_PATH.name),
        }
        np.savez_compressed(
            path,
            coeffs=np.stack(dump_coeffs),
            collisions=np.array(dump_collisions, dtype=np.int64),
            harmonics=harmonics,
            coverage=np.array(dump_coverage, dtype=np.float64),
            meta=np.array(json.dumps(meta, ensure_ascii=False)),
        )
        DUMP_MANIFEST.append({"file": str(path.relative_to(base.OUT_DIR)), "case_id": case.case_id, "R_input": float(r_value), "records": len(dump_collisions), "coverage_min": float(np.min(dump_coverage))})
    return rows


def main() -> None:
    global DUMP_ENABLED, DUMP_MAX_N, DUMP_STRIDE

    # 計測専用フラグを argv から分離してから、残りを基底パーサへ渡す
    argv = sys.argv[1:]
    filtered: List[str] = []
    i = 0
    while i < len(argv):
        token = argv[i]
        if token == "--no-dump":
            DUMP_ENABLED = False
        elif token == "--dump-max-n":
            DUMP_MAX_N = int(argv[i + 1]); i += 1
        elif token == "--dump-stride":
            DUMP_STRIDE = max(1, int(argv[i + 1])); i += 1
        else:
            filtered.append(token)
        i += 1
    sys.argv = [sys.argv[0], *filtered]

    cli_args = base.parse_args()
    if cli_args.refine_minima:
        raise SystemExit("--refine-minima is not supported by the instrumented runner; use the base runner")
    cli_params = base.Params()
    selected_cases = base.selected_cases_from_args(cli_params, cli_args)
    selected_r = base.selected_r_values_from_args(cli_params, cli_args)
    selected_output_dir = base.output_dir_from_args(cli_args)

    # run_case を計測版に差し替え（基底モジュールのグローバル参照経由で有効化）
    base.run_case = instrumented_run_case

    data = base.run(
        selected_cases=selected_cases,
        selected_r_values=selected_r,
        output_dir=selected_output_dir,
        max_collision=cli_args.max_collision,
        make_plots=not cli_args.no_plots,
        run_id=cli_args.run_id,
        file_stem=cli_args.file_stem,
        fixed_l_norm=cli_args.fixed_l_norm,
        fixed_n_norm=cli_args.fixed_n_norm,
    )

    if DUMP_ENABLED:
        manifest_path = base.OUT_DIR / "harmonic_dump_v1" / "harmonic_dump_manifest_v1.json"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(
            json.dumps({"runner": Path(__file__).name, "dump_max_n": DUMP_MAX_N, "dump_stride": DUMP_STRIDE, "entries": DUMP_MANIFEST}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        print(json.dumps({"harmonic_dump_files": len(DUMP_MANIFEST), "manifest": str(manifest_path)}, ensure_ascii=False))
    print(json.dumps({"file_stem": data.get("file_stem"), "outputs": data.get("outputs")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
