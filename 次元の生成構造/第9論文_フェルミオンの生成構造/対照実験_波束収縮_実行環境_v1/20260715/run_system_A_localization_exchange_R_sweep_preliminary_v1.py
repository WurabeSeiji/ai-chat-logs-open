from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
SOURCE_PATH = BASE_DIR.parent / "20260713" / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py"
DEFAULT_OUT_DIR = BASE_DIR / "system_A_localization_exchange_R_sweep_result_v1"
OUT_DIR = DEFAULT_OUT_DIR
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_source_module() -> Any:
    spec = importlib.util.spec_from_file_location("exchange_scattering_source_v1", SOURCE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load source module: {SOURCE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


src = load_source_module()


def configure_output_dir(path: Path) -> None:
    global OUT_DIR
    OUT_DIR = path
    OUT_DIR.mkdir(parents=True, exist_ok=True)


@dataclass
class Params:
    chi_grid_n: int = 512
    eta_grid_n: int = 16
    high_n: int = 63
    recursive_collision_count: int = 256
    coarse_r_min: float = 0.600
    coarse_r_max: float = 0.900
    coarse_r_step: float = 0.010
    fine_r_min: float = 0.680
    fine_r_max: float = 0.710
    fine_r_step: float = 0.001
    r137: float = 1.0 - math.sqrt(4.0 * math.pi / 137.035999084)
    r128: float = 1.0 - math.sqrt(4.0 * math.pi / 128.0)
    complete_r_min: float = 0.600
    complete_r_max: float = 0.900
    complete_r_step: float = 0.001
    refine_half_width: float = 0.001
    refine_r_step: float = 0.00001
    refine_top_k: int = 3
    pairs: Tuple[Tuple[int, int], ...] = (
        (1, 1),
        (1, 2),
        (1, 3),
        (1, 5),
        (1, 10),
        (1, 15),
        (1, 31),
        (1, 63),
        (31, 31),
        (63, 63),
    )
    tail_window: int = 24


@dataclass(frozen=True)
class HarmonicCase:
    mode: str
    state_family: str
    n_a: int
    n_b: int
    packet_a: Tuple[int, ...]
    packet_b: Tuple[int, ...]
    packet_a_weights: Tuple[float, ...]
    packet_b_weights: Tuple[float, ...]
    packet_a_phases: Tuple[float, ...]
    packet_b_phases: Tuple[float, ...]
    packet_a_wavelength_scales: Tuple[float, ...]
    packet_b_wavelength_scales: Tuple[float, ...]
    packet_a_shift: float
    packet_b_shift: float

    @property
    def packet_a_text(self) -> str:
        return packet_text(self.packet_a)

    @property
    def packet_b_text(self) -> str:
        return packet_text(self.packet_b)

    @property
    def packet_a_weight_text(self) -> str:
        return float_packet_text(self.packet_a_weights)

    @property
    def packet_b_weight_text(self) -> str:
        return float_packet_text(self.packet_b_weights)

    @property
    def packet_a_phase_text(self) -> str:
        return float_packet_text(self.packet_a_phases)

    @property
    def packet_b_phase_text(self) -> str:
        return float_packet_text(self.packet_b_phases)

    @property
    def packet_a_wavelength_scale_text(self) -> str:
        return float_packet_text(self.packet_a_wavelength_scales)

    @property
    def packet_b_wavelength_scale_text(self) -> str:
        return float_packet_text(self.packet_b_wavelength_scales)

    @property
    def packet_a_spec_text(self) -> str:
        return packet_spec_text(
            self.packet_a,
            self.packet_a_weights,
            self.packet_a_phases,
            self.packet_a_wavelength_scales,
            self.packet_a_shift,
        )

    @property
    def packet_b_spec_text(self) -> str:
        return packet_spec_text(
            self.packet_b,
            self.packet_b_weights,
            self.packet_b_phases,
            self.packet_b_wavelength_scales,
            self.packet_b_shift,
        )

    @property
    def case_id(self) -> str:
        return f"{self.mode}|A={self.packet_a_spec_text}|B={self.packet_b_spec_text}"


def packet_text(packet: Tuple[int, ...]) -> str:
    return ",".join(str(value) for value in packet)


def float_packet_text(values: Tuple[float, ...]) -> str:
    return ",".join(f"{float(value):.12g}" for value in values)


def default_weights(packet: Tuple[int, ...]) -> Tuple[float, ...]:
    return tuple(1.0 for _ in packet)


def default_phases(packet: Tuple[int, ...]) -> Tuple[float, ...]:
    return tuple(0.0 for _ in packet)


def default_wavelength_scales(packet: Tuple[int, ...]) -> Tuple[float, ...]:
    return tuple(1.0 for _ in packet)


def packet_spec_text(
    packet: Tuple[int, ...],
    weights: Tuple[float, ...],
    phases: Tuple[float, ...],
    wavelength_scales: Tuple[float, ...],
    shift: float,
) -> str:
    parts = [packet_text(packet)]
    if any(abs(float(value) - 1.0) > 1.0e-14 for value in weights):
        parts.append(f"w={float_packet_text(weights)}")
    if any(abs(float(value)) > 1.0e-14 for value in phases):
        parts.append(f"p={float_packet_text(phases)}")
    if any(abs(float(value) - 1.0) > 1.0e-14 for value in wavelength_scales):
        parts.append(f"lambda={float_packet_text(wavelength_scales)}")
    if abs(float(shift)) > 1.0e-14:
        parts.append(f"s={float(shift):.12g}")
    return ";".join(parts)


def validate_packet_parameters(
    packet: Tuple[int, ...],
    weights: Tuple[float, ...],
    phases: Tuple[float, ...],
    wavelength_scales: Tuple[float, ...],
    label: str,
) -> None:
    if len(weights) != len(packet):
        raise ValueError(f"{label} weights length must match packet length")
    if len(phases) != len(packet):
        raise ValueError(f"{label} phases length must match packet length")
    if len(wavelength_scales) != len(packet):
        raise ValueError(f"{label} wavelength scales length must match packet length")
    if not any(abs(float(weight)) > 0.0 for weight in weights):
        raise ValueError(f"{label} weights must contain at least one nonzero value")
    if any(float(scale) <= 0.0 for scale in wavelength_scales):
        raise ValueError(f"{label} wavelength scales must be positive")


def odd_kernel_case(pair: Tuple[int, int]) -> HarmonicCase:
    packet_a = (int(pair[0]),)
    packet_b = (int(pair[1]),)
    return HarmonicCase(
        mode="odd_kernel",
        state_family="odd_kernel",
        n_a=int(pair[0]),
        n_b=int(pair[1]),
        packet_a=packet_a,
        packet_b=packet_b,
        packet_a_weights=default_weights(packet_a),
        packet_b_weights=default_weights(packet_b),
        packet_a_phases=default_phases(packet_a),
        packet_b_phases=default_phases(packet_b),
        packet_a_wavelength_scales=default_wavelength_scales(packet_a),
        packet_b_wavelength_scales=default_wavelength_scales(packet_b),
        packet_a_shift=0.0,
        packet_b_shift=0.0,
    )


def explicit_packet_case(
    mode: str,
    packet_a: Tuple[int, ...],
    packet_b: Tuple[int, ...],
    packet_a_weights: Optional[Tuple[float, ...]] = None,
    packet_b_weights: Optional[Tuple[float, ...]] = None,
    packet_a_phases: Optional[Tuple[float, ...]] = None,
    packet_b_phases: Optional[Tuple[float, ...]] = None,
    packet_a_wavelength_scales: Optional[Tuple[float, ...]] = None,
    packet_b_wavelength_scales: Optional[Tuple[float, ...]] = None,
    packet_a_shift: float = 0.0,
    packet_b_shift: float = 0.0,
) -> HarmonicCase:
    if not packet_a or not packet_b:
        raise ValueError("harmonic packet must not be empty")
    weights_a = packet_a_weights if packet_a_weights is not None else default_weights(packet_a)
    weights_b = packet_b_weights if packet_b_weights is not None else default_weights(packet_b)
    phases_a = packet_a_phases if packet_a_phases is not None else default_phases(packet_a)
    phases_b = packet_b_phases if packet_b_phases is not None else default_phases(packet_b)
    wavelengths_a = packet_a_wavelength_scales if packet_a_wavelength_scales is not None else default_wavelength_scales(packet_a)
    wavelengths_b = packet_b_wavelength_scales if packet_b_wavelength_scales is not None else default_wavelength_scales(packet_b)
    validate_packet_parameters(packet_a, weights_a, phases_a, wavelengths_a, "A packet")
    validate_packet_parameters(packet_b, weights_b, phases_b, wavelengths_b, "B packet")
    return HarmonicCase(
        mode=mode,
        state_family="explicit_packet",
        n_a=max(packet_a),
        n_b=max(packet_b),
        packet_a=tuple(int(value) for value in packet_a),
        packet_b=tuple(int(value) for value in packet_b),
        packet_a_weights=tuple(float(value) for value in weights_a),
        packet_b_weights=tuple(float(value) for value in weights_b),
        packet_a_phases=tuple(float(value) for value in phases_a),
        packet_b_phases=tuple(float(value) for value in phases_b),
        packet_a_wavelength_scales=tuple(float(value) for value in wavelengths_a),
        packet_b_wavelength_scales=tuple(float(value) for value in wavelengths_b),
        packet_a_shift=float(packet_a_shift),
        packet_b_shift=float(packet_b_shift),
    )


def default_cases(params: Params) -> Tuple[HarmonicCase, ...]:
    return tuple(odd_kernel_case(pair) for pair in params.pairs)


def built_in_packet_cases() -> Tuple[HarmonicCase, ...]:
    return (
        explicit_packet_case("even_packet", (1,), (1, 2, 4, 6)),
        explicit_packet_case("alternating_packet_1", (1,), (1, 2, 3, 4, 5)),
        explicit_packet_case("alternating_packet_2", (1,), (1, 3, 4, 5, 6)),
    )


def build_source_params(params: Params) -> Any:
    return src.Params(
        chi_grid_n=params.chi_grid_n,
        eta_grid_n=params.eta_grid_n,
        high_n=params.high_n,
        recursive_collision_count=params.recursive_collision_count,
    )


class MetricContext:
    def __init__(self, source_params: Any):
        self.params = source_params
        self.freqs = np.fft.fftfreq(source_params.chi_grid_n, d=1.0 / source_params.chi_grid_n)
        self.max_n = min(source_params.chi_grid_n // 2, source_params.high_n + 2)
        self.indices_by_abs_n: Dict[int, List[int]] = {}
        rounded = np.rint(self.freqs).astype(int)
        for n_abs in range(self.max_n + 1):
            indices = [int(i) for i, freq in enumerate(rounded) if abs(int(freq)) == n_abs]
            self.indices_by_abs_n[n_abs] = indices

    def norm2(self, vector: np.ndarray) -> float:
        return float(np.vdot(vector, vector).real)

    def fft_power(self, vector: np.ndarray) -> Tuple[np.ndarray, float]:
        denom = self.norm2(vector)
        arr = vector.reshape(self.params.chi_grid_n, self.params.eta_grid_n)
        transformed = np.fft.fft(arr, axis=0, norm="ortho")
        power = np.sum(np.abs(transformed) ** 2, axis=1)
        return power, denom

    def harmonic_distribution(self, vector: np.ndarray) -> Dict[int, float]:
        power, denom = self.fft_power(vector)
        if denom <= 0.0:
            return {0: 1.0}
        raw: Dict[int, float] = {}
        total = 0.0
        for n_abs, indices in self.indices_by_abs_n.items():
            amount = float(np.sum(power[indices]).real) / denom if indices else 0.0
            if amount > 1.0e-14:
                raw[n_abs] = max(amount, 0.0)
                total += raw[n_abs]
        if total <= 0.0:
            return {0: 1.0}
        return {k: v / total for k, v in raw.items()}

    def p_chi(self, vector: np.ndarray) -> float:
        power, denom = self.fft_power(vector)
        if denom <= 0.0:
            return float("nan")
        return float(np.sum(self.freqs * power).real / denom)


def r_values(params: Params) -> List[float]:
    values = [0.0, 0.5, 1.0, params.r137, params.r128]
    coarse_steps = int(round((params.coarse_r_max - params.coarse_r_min) / params.coarse_r_step))
    values.extend(params.coarse_r_min + params.coarse_r_step * i for i in range(coarse_steps + 1))
    fine_steps = int(round((params.fine_r_max - params.fine_r_min) / params.fine_r_step))
    values.extend(params.fine_r_min + params.fine_r_step * i for i in range(fine_steps + 1))
    return sorted({round(float(v), 12) for v in values})


def uniform_r_values(start: float, stop: float, step: float) -> List[float]:
    if step <= 0.0:
        raise ValueError("R step must be positive")
    if stop < start:
        raise ValueError("R stop must be greater than or equal to R start")
    count = int(math.floor((stop - start) / step + 1.0e-12))
    values = [round(start + step * i, 12) for i in range(count + 1)]
    if not values or abs(values[-1] - stop) > step * 0.5:
        values.append(round(stop, 12))
    return sorted({float(value) for value in values})


def parse_pair_token(token: str) -> Tuple[int, int]:
    normalized = token.strip().replace("x", ":").replace("-", ":")
    parts = normalized.split(":")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(f"pair must be like 1:2, got {token!r}")
    return int(parts[0]), int(parts[1])


def parse_pairs_text(text: str) -> Tuple[Tuple[int, int], ...]:
    if not text.strip():
        return tuple()
    return tuple(parse_pair_token(token) for token in text.split(",") if token.strip())


def parse_float_list(text: str) -> List[float]:
    if not text.strip():
        return []
    return sorted({round(float(token.strip()), 12) for token in text.split(",") if token.strip()})


def parse_packet_text(text: str) -> Tuple[int, ...]:
    if not text.strip():
        return tuple()
    return tuple(int(token.strip()) for token in text.split(",") if token.strip())


def parse_weight_text(text: str) -> Tuple[float, ...]:
    if not text.strip():
        return tuple()
    return tuple(float(token.strip()) for token in text.split(",") if token.strip())


def packet_parameter_or_default(values: Optional[str], packet: Tuple[int, ...], default: float) -> Tuple[float, ...]:
    if values is None:
        return tuple(float(default) for _ in packet)
    parsed = parse_weight_text(values)
    if len(parsed) == 1 and len(packet) > 1:
        return tuple(float(parsed[0]) for _ in packet)
    return parsed


def selected_pairs_from_args(params: Params, args: argparse.Namespace) -> Tuple[Tuple[int, int], ...]:
    if args.pairs:
        return parse_pairs_text(args.pairs)
    start = args.pair_start if args.pair_start is not None else 0
    stop = args.pair_stop if args.pair_stop is not None else len(params.pairs)
    return tuple(params.pairs[start:stop])


def selected_cases_from_args(params: Params, args: argparse.Namespace) -> Tuple[HarmonicCase, ...]:
    if args.packet_a or args.packet_b:
        packet_a = parse_packet_text(args.packet_a or "1")
        packet_b = parse_packet_text(args.packet_b or "")
        packet_a_weights = packet_parameter_or_default(args.packet_a_weights, packet_a, 1.0)
        packet_b_weights = packet_parameter_or_default(args.packet_b_weights, packet_b, 1.0)
        packet_a_phases = packet_parameter_or_default(args.packet_a_phases, packet_a, 0.0)
        packet_b_phases = packet_parameter_or_default(args.packet_b_phases, packet_b, 0.0)
        packet_a_wavelength_scales = packet_parameter_or_default(args.packet_a_wavelength_scales, packet_a, 1.0)
        packet_b_wavelength_scales = packet_parameter_or_default(args.packet_b_wavelength_scales, packet_b, 1.0)
        return (
            explicit_packet_case(
                "custom_packet",
                packet_a,
                packet_b,
                packet_a_weights=packet_a_weights,
                packet_b_weights=packet_b_weights,
                packet_a_phases=packet_a_phases,
                packet_b_phases=packet_b_phases,
                packet_a_wavelength_scales=packet_a_wavelength_scales,
                packet_b_wavelength_scales=packet_b_wavelength_scales,
                packet_a_shift=float(args.packet_a_shift),
                packet_b_shift=float(args.packet_b_shift),
            ),
        )

    pairs = selected_pairs_from_args(params, args)
    single_cases = tuple(odd_kernel_case(pair) for pair in pairs)
    packet_cases = built_in_packet_cases()

    if args.harmonic_mode == "single":
        return single_cases
    if args.harmonic_mode == "even-packet":
        return (packet_cases[0],)
    if args.harmonic_mode == "alternating-packet-1":
        return (packet_cases[1],)
    if args.harmonic_mode == "alternating-packet-2":
        return (packet_cases[2],)
    if args.harmonic_mode == "all":
        return single_cases + packet_cases
    raise ValueError(f"unknown harmonic mode: {args.harmonic_mode}")


def selected_r_values_from_args(params: Params, args: argparse.Namespace) -> List[float]:
    if args.r_values:
        values = parse_float_list(args.r_values)
    elif args.complete_sweep:
        start = args.complete_r_min if args.complete_r_min is not None else params.complete_r_min
        stop = args.complete_r_max if args.complete_r_max is not None else params.complete_r_max
        step = args.complete_r_step if args.complete_r_step is not None else params.complete_r_step
        values = uniform_r_values(float(start), float(stop), float(step))
    else:
        values = r_values(params)
    if args.r_min is not None:
        values = [value for value in values if value >= args.r_min]
    if args.r_max is not None:
        values = [value for value in values if value <= args.r_max]
    if not values:
        raise ValueError("no R values selected")
    return values


def output_dir_from_args(args: argparse.Namespace) -> Path:
    if args.output_dir:
        return Path(args.output_dir).expanduser().resolve()
    if args.run_id:
        return DEFAULT_OUT_DIR / args.run_id
    return DEFAULT_OUT_DIR


def safe_slug(text: str, max_len: int = 80) -> str:
    allowed = []
    for ch in text:
        if ch.isalnum() or ch in {"-", "_"}:
            allowed.append(ch)
        else:
            allowed.append("-")
    slug = "".join(allowed).strip("-_")
    while "--" in slug:
        slug = slug.replace("--", "-")
    if len(slug) <= max_len:
        return slug
    digest = hashlib.sha1(slug.encode("utf-8")).hexdigest()[:8]
    return f"{slug[: max_len - 9]}-{digest}"


def compact_float(value: float) -> str:
    return f"{value:.12g}".replace(".", "p").replace("-", "m")


def cases_slug(cases: Tuple[HarmonicCase, ...]) -> str:
    if all(case.mode == "odd_kernel" and case.packet_a == (1,) for case in cases):
        b_values = "-".join(case.packet_b_text.replace(",", "p") for case in cases)
        return safe_slug(f"odd_B{b_values}", max_len=90)
    if len(cases) == 1:
        case = cases[0]
        return safe_slug(f"{case.mode}_A{case.packet_a_spec_text}_B{case.packet_b_spec_text}", max_len=90)
    digest_src = "|".join(case.case_id for case in cases)
    digest = hashlib.sha1(digest_src.encode("utf-8")).hexdigest()[:8]
    modes = "-".join(sorted({case.mode for case in cases}))
    return safe_slug(f"cases{len(cases)}_{modes}_{digest}", max_len=90)


def r_values_slug(values: List[float], params: Params) -> str:
    default_values = r_values(params)
    if values == default_values:
        return "Rdefault"
    if len(values) <= 4:
        return safe_slug("R" + "-".join(compact_float(value) for value in values), max_len=90)
    digest_src = ",".join(f"{float(value):.12g}" for value in values)
    digest = hashlib.sha1(digest_src.encode("utf-8")).hexdigest()[:8]
    return safe_slug(f"R{compact_float(min(values))}-{compact_float(max(values))}_n{len(values)}_{digest}", max_len=90)


def build_file_stem(cases: Tuple[HarmonicCase, ...], values: List[float], params: Params, run_id: Optional[str], explicit_stem: Optional[str]) -> str:
    if explicit_stem:
        return safe_slug(explicit_stem, max_len=120)
    parts = ["system_A", cases_slug(cases), r_values_slug(values, params), f"C{params.recursive_collision_count}"]
    if run_id:
        parts.insert(1, safe_slug(run_id, max_len=40))
    return safe_slug("_".join(parts), max_len=150)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="System A localization exchange R sweep")
    parser.add_argument("--pairs", help="comma-separated pair list, e.g. 1:2,1:3,1:63")
    parser.add_argument("--pair-start", type=int, help="zero-based start index in default pair list")
    parser.add_argument("--pair-stop", type=int, help="zero-based stop index in default pair list")
    parser.add_argument(
        "--harmonic-mode",
        choices=("single", "even-packet", "alternating-packet-1", "alternating-packet-2", "all"),
        default="single",
        help="harmonic construction mode",
    )
    parser.add_argument("--packet-a", help="custom explicit harmonic packet for A, e.g. 1 or 1,3,5")
    parser.add_argument("--packet-b", help="custom explicit harmonic packet for B, e.g. 1,2,4,6")
    parser.add_argument("--packet-a-weights", help="comma-separated A harmonic weights; one value is broadcast")
    parser.add_argument("--packet-b-weights", help="comma-separated B harmonic weights; one value is broadcast")
    parser.add_argument("--packet-a-phases", help="comma-separated A harmonic phase offsets in radians; one value is broadcast")
    parser.add_argument("--packet-b-phases", help="comma-separated B harmonic phase offsets in radians; one value is broadcast")
    parser.add_argument("--packet-a-wavelength-scales", help="comma-separated A wavelength scales; frequency uses n / scale")
    parser.add_argument("--packet-b-wavelength-scales", help="comma-separated B wavelength scales; frequency uses n / scale")
    parser.add_argument("--packet-a-shift", type=float, default=0.0, help="A packet translation in chi radians")
    parser.add_argument("--packet-b-shift", type=float, default=0.0, help="B packet translation in chi radians")
    parser.add_argument("--r-min", type=float, help="minimum R to include")
    parser.add_argument("--r-max", type=float, help="maximum R to include")
    parser.add_argument("--r-values", help="comma-separated explicit R values")
    parser.add_argument("--complete-sweep", action="store_true", help="use a uniform R sweep without injecting R_137/R_128 probe points")
    parser.add_argument("--complete-r-min", type=float, help="uniform complete-sweep minimum R")
    parser.add_argument("--complete-r-max", type=float, help="uniform complete-sweep maximum R")
    parser.add_argument("--complete-r-step", type=float, help="uniform complete-sweep R step")
    parser.add_argument("--refine-minima", action="store_true", help="rerun around local minima found in the first sweep")
    parser.add_argument("--refine-top-k", type=int, help="number of local minima to refine per case")
    parser.add_argument("--refine-half-width", type=float, help="half width around each local minimum for refined sweep")
    parser.add_argument("--refine-r-step", type=float, help="R step for refined sweep")
    parser.add_argument("--fixed-l-norm", type=float, help="fixed normalization denominator for L gap scores")
    parser.add_argument("--fixed-n-norm", type=float, help="fixed normalization denominator for N_eff gap scores")
    parser.add_argument("--max-collision", type=int, help="override recursive collision count")
    parser.add_argument("--run-id", help="write outputs under the default result directory with this subdirectory name")
    parser.add_argument("--file-stem", help="explicit output file stem")
    parser.add_argument("--output-dir", help="write outputs to this exact directory")
    parser.add_argument("--no-plots", action="store_true", help="skip png generation")
    return parser.parse_args()


def distribution_similarity(a: Dict[int, float], b: Dict[int, float]) -> float:
    keys = set(a) | set(b)
    va = np.asarray([float(a.get(k, 0.0)) for k in keys], dtype=float)
    vb = np.asarray([float(b.get(k, 0.0)) for k in keys], dtype=float)
    na = float(np.linalg.norm(va))
    nb = float(np.linalg.norm(vb))
    if na <= 0.0 or nb <= 0.0:
        return float("nan")
    return float(max(0.0, min(1.0, np.dot(va, vb) / (na * nb))))


def explicit_packet_kernel(
    source_params: Any,
    harmonics: Tuple[int, ...],
    weights: Tuple[float, ...],
    phases: Tuple[float, ...],
    wavelength_scales: Tuple[float, ...],
    shift: float,
) -> np.ndarray:
    chi, _ = src.make_grids(source_params)
    u = chi - source_params.chi_center
    kernel = np.zeros_like(chi, dtype=float)
    for harmonic, weight, phase, wavelength_scale in zip(harmonics, weights, phases, wavelength_scales):
        effective_frequency = float(harmonic) / float(wavelength_scale)
        kernel += float(weight) * np.cos(effective_frequency * (u - shift) + float(phase))
    norm = math.sqrt(float(sum(float(weight) ** 2 for weight in weights)))
    return kernel / norm if norm > 0.0 else kernel


def make_explicit_packet_state(
    source_params: Any,
    harmonics: Tuple[int, ...],
    weights: Tuple[float, ...],
    phases: Tuple[float, ...],
    wavelength_scales: Tuple[float, ...],
    shift: float,
    q: float,
    m: int,
    hair_enabled: bool,
    amplitude: float,
) -> np.ndarray:
    chi, eta = src.make_grids(source_params)
    kernel = explicit_packet_kernel(source_params, harmonics, weights, phases, wavelength_scales, shift)
    phase_chi = np.exp(1j * q * source_params.p0 * (chi - source_params.chi_center))
    eta_phase = np.exp(1j * m * eta) if hair_enabled else np.ones_like(eta, dtype=complex)
    psi = (kernel * phase_chi)[:, None] * eta_phase[None, :]
    return amplitude * src.normalize(psi.reshape(-1))


def make_case_state(source_params: Any, case: HarmonicCase, side: str, hair_enabled: bool) -> np.ndarray:
    if side == "A":
        q = source_params.q_A
        m = source_params.m_A
        amplitude = source_params.A_A
        n_value = case.n_a
        packet = case.packet_a
        weights = case.packet_a_weights
        phases = case.packet_a_phases
        wavelength_scales = case.packet_a_wavelength_scales
        shift = case.packet_a_shift
    elif side == "B":
        q = source_params.q_B
        m = source_params.m_B
        amplitude = source_params.A_B
        n_value = case.n_b
        packet = case.packet_b
        weights = case.packet_b_weights
        phases = case.packet_b_phases
        wavelength_scales = case.packet_b_wavelength_scales
        shift = case.packet_b_shift
    else:
        raise ValueError(f"unknown side: {side}")

    if case.state_family == "odd_kernel":
        return src.make_state(source_params, n_value, q, m, hair_enabled, amplitude)
    if case.state_family == "explicit_packet":
        return make_explicit_packet_state(source_params, packet, weights, phases, wavelength_scales, shift, q, m, hair_enabled, amplitude)
    raise ValueError(f"unknown state family: {case.state_family}")


def row_for_state(
    source_params: Any,
    metrics: MetricContext,
    case: HarmonicCase,
    r_value: float,
    t: complex,
    r: complex,
    T: float,
    R: float,
    collision: int,
    channel: str,
    vector: np.ndarray,
    h_a0: Dict[int, float],
    h_b0: Dict[int, float],
    initial_a: np.ndarray,
    initial_b: np.ndarray,
) -> Dict[str, Any]:
    h = metrics.harmonic_distribution(vector)
    n_eff, n_eff_2 = src.effective_n(h)
    l_value = src.localization(vector)
    p_chi = metrics.p_chi(vector)
    sim_a = distribution_similarity(h, h_a0)
    sim_b = distribution_similarity(h, h_b0)
    return {
        "case_id": case.case_id,
        "mode": case.mode,
        "state_family": case.state_family,
        "harmonic_packet_A": case.packet_a_text,
        "harmonic_packet_B": case.packet_b_text,
        "harmonic_packet_A_weights": case.packet_a_weight_text,
        "harmonic_packet_B_weights": case.packet_b_weight_text,
        "harmonic_packet_A_phases": case.packet_a_phase_text,
        "harmonic_packet_B_phases": case.packet_b_phase_text,
        "harmonic_packet_A_wavelength_scales": case.packet_a_wavelength_scale_text,
        "harmonic_packet_B_wavelength_scales": case.packet_b_wavelength_scale_text,
        "harmonic_packet_A_shift": case.packet_a_shift,
        "harmonic_packet_B_shift": case.packet_b_shift,
        "N_A": case.n_a,
        "N_B": case.n_b,
        "R_input": r_value,
        "R": R,
        "T": T,
        "Delta_F": src.delta_from_reflection_rate(r_value),
        "collision": collision,
        "channel": channel,
        "L": l_value,
        "N_eff": n_eff,
        "N_eff_2": n_eff_2,
        "p_chi": p_chi,
        "origin_A": src.projection_weight(vector, initial_a),
        "origin_B": src.projection_weight(vector, initial_b),
        "sim_to_A0": sim_a,
        "sim_to_B0": sim_b,
    }


def run_case(source_params: Any, metrics: MetricContext, case: HarmonicCase, r_value: float, max_collision: int) -> List[Dict[str, Any]]:
    delta_f = src.delta_from_reflection_rate(r_value)
    t, r, T, R = src.scattering_coefficients(delta_f)
    hair_enabled = True
    a = make_case_state(source_params, case, "A", hair_enabled)
    b = make_case_state(source_params, case, "B", hair_enabled)
    initial_a = a.copy()
    initial_b = b.copy()
    h_a0 = metrics.harmonic_distribution(initial_a)
    h_b0 = metrics.harmonic_distribution(initial_b)
    rows: List[Dict[str, Any]] = []
    for collision in range(max_collision + 1):
        rows.append(row_for_state(source_params, metrics, case, r_value, t, r, T, R, collision, "A_channel", a, h_a0, h_b0, initial_a, initial_b))
        rows.append(row_for_state(source_params, metrics, case, r_value, t, r, T, R, collision, "B_channel", b, h_a0, h_b0, initial_a, initial_b))
        if collision >= max_collision:
            break
        a_next = src.normalize(r * a + t * b)
        b_next = src.normalize(t * a + r * b)
        a, b = a_next, b_next
    return rows


def pair_rows_by_collision(rows: Iterable[Dict[str, Any]]) -> Dict[int, Dict[str, Dict[str, Any]]]:
    grouped: Dict[int, Dict[str, Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(int(row["collision"]), {})[str(row["channel"])] = row
    return grouped


def summarize_case(rows: List[Dict[str, Any]], params: Params) -> Dict[str, Any]:
    grouped = pair_rows_by_collision(rows)
    first_row = rows[0]
    n_a = int(first_row["N_A"])
    n_b = int(first_row["N_B"])
    r_value = float(first_row["R_input"])
    pair_kind = "one_side_high_harmonic" if n_a != n_b else "same_harmonic_control"
    if params.fine_r_min <= r_value <= params.fine_r_max:
        sweep_region = "fine"
    elif params.coarse_r_min <= r_value <= params.coarse_r_max:
        sweep_region = "coarse"
    else:
        sweep_region = "control"
    records: List[Dict[str, float]] = []
    for collision, by_channel in grouped.items():
        if "A_channel" not in by_channel or "B_channel" not in by_channel:
            continue
        a = by_channel["A_channel"]
        b = by_channel["B_channel"]
        l_gap = abs(float(a["L"]) - float(b["L"]))
        n_gap = abs(float(a["N_eff"]) - float(b["N_eff"]))
        transfer_b_to_a = float(a["sim_to_B0"])
        transfer_a_to_b = float(b["sim_to_A0"])
        records.append(
            {
                "collision": float(collision),
                "L_gap": l_gap,
                "N_eff_gap": n_gap,
                "B_to_A_transfer": transfer_b_to_a,
                "A_to_B_transfer": transfer_a_to_b,
            }
        )
    max_l = max([record["L_gap"] for record in records] or [1.0])
    max_n = max([record["N_eff_gap"] for record in records] or [1.0])
    max_l = max(max_l, 1.0e-300)
    max_n = max(max_n, 1.0e-300)
    for record in records:
        record["joint_score"] = (record["L_gap"] / max_l) + (record["N_eff_gap"] / max_n) + (1.0 - record["B_to_A_transfer"])

    min_l = min(records, key=lambda item: item["L_gap"])
    min_n = min(records, key=lambda item: item["N_eff_gap"])
    max_transfer = max(records, key=lambda item: item["B_to_A_transfer"])
    min_joint = min(records, key=lambda item: item["joint_score"])
    tail_from = max(0, params.recursive_collision_count - params.tail_window)
    tail = [record for record in records if int(record["collision"]) >= tail_from]
    return {
        "case_id": str(first_row["case_id"]),
        "mode": str(first_row["mode"]),
        "state_family": str(first_row["state_family"]),
        "harmonic_packet_A": str(first_row["harmonic_packet_A"]),
        "harmonic_packet_B": str(first_row["harmonic_packet_B"]),
        "harmonic_packet_A_weights": str(first_row["harmonic_packet_A_weights"]),
        "harmonic_packet_B_weights": str(first_row["harmonic_packet_B_weights"]),
        "harmonic_packet_A_phases": str(first_row["harmonic_packet_A_phases"]),
        "harmonic_packet_B_phases": str(first_row["harmonic_packet_B_phases"]),
        "harmonic_packet_A_wavelength_scales": str(first_row["harmonic_packet_A_wavelength_scales"]),
        "harmonic_packet_B_wavelength_scales": str(first_row["harmonic_packet_B_wavelength_scales"]),
        "harmonic_packet_A_shift": float(first_row["harmonic_packet_A_shift"]),
        "harmonic_packet_B_shift": float(first_row["harmonic_packet_B_shift"]),
        "N_A": n_a,
        "N_B": n_b,
        "pair_kind": pair_kind,
        "sweep_region": sweep_region,
        "R": r_value,
        "T": float(first_row["T"]),
        "Delta_F": float(first_row["Delta_F"]),
        "L_gap_min": min_l["L_gap"],
        "L_gap_min_collision": int(min_l["collision"]),
        "N_eff_gap_at_L_gap_min": min_l["N_eff_gap"],
        "N_eff_gap_min": min_n["N_eff_gap"],
        "N_eff_gap_min_collision": int(min_n["collision"]),
        "max_B_to_A_transfer": max_transfer["B_to_A_transfer"],
        "B_to_A_transfer_collision": int(max_transfer["collision"]),
        "A_to_B_transfer_at_B_to_A_max": max_transfer["A_to_B_transfer"],
        "joint_score_min": min_joint["joint_score"],
        "joint_score_min_collision": int(min_joint["collision"]),
        "tail_from_collision": tail_from,
        "tail_L_gap_min": min(record["L_gap"] for record in tail) if tail else float("nan"),
        "tail_L_gap_max": max(record["L_gap"] for record in tail) if tail else float("nan"),
        "tail_N_eff_gap_min": min(record["N_eff_gap"] for record in tail) if tail else float("nan"),
        "tail_N_eff_gap_max": max(record["N_eff_gap"] for record in tail) if tail else float("nan"),
    }


def normalize_column(rows: List[Dict[str, Any]], key: str) -> Dict[Tuple[str, float], float]:
    values = [float(row[key]) for row in rows if math.isfinite(float(row[key]))]
    max_value = max(values) if values else 1.0
    max_value = max(max_value, 1.0e-300)
    return {(str(row["case_id"]), float(row["R"])): float(row[key]) / max_value for row in rows}


def add_pair_level_scores(
    summaries: List[Dict[str, Any]],
    fixed_l_norm: Optional[float] = None,
    fixed_n_norm: Optional[float] = None,
) -> None:
    by_pair: Dict[str, List[Dict[str, Any]]] = {}
    for row in summaries:
        by_pair.setdefault(str(row["case_id"]), []).append(row)
    for pair_rows in by_pair.values():
        if fixed_l_norm is not None:
            fixed_l = max(float(fixed_l_norm), 1.0e-300)
            l_norm = {(str(row["case_id"]), float(row["R"])): float(row["L_gap_min"]) / fixed_l for row in pair_rows}
        else:
            l_norm = normalize_column(pair_rows, "L_gap_min")
        if fixed_n_norm is not None:
            fixed_n = max(float(fixed_n_norm), 1.0e-300)
            n_norm = {(str(row["case_id"]), float(row["R"])): float(row["N_eff_gap_min"]) / fixed_n for row in pair_rows}
        else:
            n_norm = normalize_column(pair_rows, "N_eff_gap_min")
        for row in pair_rows:
            key = (str(row["case_id"]), float(row["R"]))
            row["joint_R_score"] = l_norm[key] + n_norm[key] + (1.0 - float(row["max_B_to_A_transfer"]))


def best_rows_for_pair(
    summaries: List[Dict[str, Any]],
    params: Params,
    fixed_l_norm: Optional[float] = None,
    fixed_n_norm: Optional[float] = None,
) -> List[Dict[str, Any]]:
    add_pair_level_scores(summaries, fixed_l_norm=fixed_l_norm, fixed_n_norm=fixed_n_norm)
    out: List[Dict[str, Any]] = []
    for case_id in sorted({str(row["case_id"]) for row in summaries}):
        all_pair_rows = [row for row in summaries if str(row["case_id"]) == case_id]
        pair_rows = [row for row in all_pair_rows if str(row.get("sweep_region")) in {"coarse", "fine"}]
        control_rows = [row for row in all_pair_rows if str(row.get("sweep_region")) == "control"]
        if not pair_rows:
            pair_rows = all_pair_rows
        best_l = min(pair_rows, key=lambda row: float(row["L_gap_min"]))
        best_n = min(pair_rows, key=lambda row: float(row["N_eff_gap_min"]))
        best_transfer = max(pair_rows, key=lambda row: float(row["max_B_to_A_transfer"]))
        best_joint = min(pair_rows, key=lambda row: float(row["joint_R_score"]))
        coarse_rows = [row for row in all_pair_rows if str(row.get("sweep_region")) == "coarse"]
        fine_rows = [row for row in all_pair_rows if str(row.get("sweep_region")) == "fine"]
        best_control = min(control_rows, key=lambda row: float(row["joint_R_score"])) if control_rows else None
        best_coarse = min(coarse_rows, key=lambda row: float(row["joint_R_score"])) if coarse_rows else None
        best_fine = min(fine_rows, key=lambda row: float(row["joint_R_score"])) if fine_rows else None
        min_score = float(best_joint["joint_R_score"])
        within_5 = [row for row in pair_rows if float(row["joint_R_score"]) <= min_score * 1.05 + 1.0e-15]
        within_10 = [row for row in pair_rows if float(row["joint_R_score"]) <= min_score * 1.10 + 1.0e-15]
        def width(rows: List[Dict[str, Any]]) -> float:
            if not rows:
                return float("nan")
            r_list = [float(row["R"]) for row in rows]
            return max(r_list) - min(r_list)
        out.append(
            {
                "case_id": case_id,
                "mode": str(best_joint["mode"]),
                "state_family": str(best_joint["state_family"]),
                "harmonic_packet_A": str(best_joint["harmonic_packet_A"]),
                "harmonic_packet_B": str(best_joint["harmonic_packet_B"]),
                "harmonic_packet_A_weights": str(best_joint["harmonic_packet_A_weights"]),
                "harmonic_packet_B_weights": str(best_joint["harmonic_packet_B_weights"]),
                "harmonic_packet_A_phases": str(best_joint["harmonic_packet_A_phases"]),
                "harmonic_packet_B_phases": str(best_joint["harmonic_packet_B_phases"]),
                "harmonic_packet_A_wavelength_scales": str(best_joint["harmonic_packet_A_wavelength_scales"]),
                "harmonic_packet_B_wavelength_scales": str(best_joint["harmonic_packet_B_wavelength_scales"]),
                "harmonic_packet_A_shift": float(best_joint["harmonic_packet_A_shift"]),
                "harmonic_packet_B_shift": float(best_joint["harmonic_packet_B_shift"]),
                "N_A": int(best_joint["N_A"]),
                "N_B": int(best_joint["N_B"]),
                "pair_kind": best_joint["pair_kind"],
                "R_star_L": float(best_l["R"]),
                "collision_at_R_star_L": int(best_l["L_gap_min_collision"]),
                "R_star_N": float(best_n["R"]),
                "collision_at_R_star_N": int(best_n["N_eff_gap_min_collision"]),
                "R_star_transfer": float(best_transfer["R"]),
                "collision_at_R_star_transfer": int(best_transfer["B_to_A_transfer_collision"]),
                "R_star_joint": float(best_joint["R"]),
                "collision_at_R_star_joint": int(best_joint["joint_score_min_collision"]),
                "joint_score_min": float(best_joint["joint_R_score"]),
                "R_band_width_5": width(within_5),
                "R_band_width_10": width(within_10),
                "distance_to_R_137": abs(float(best_joint["R"]) - params.r137),
                "distance_to_R_128": abs(float(best_joint["R"]) - params.r128),
                "best_control_R": float(best_control["R"]) if best_control is not None else float("nan"),
                "best_control_joint_score": float(best_control["joint_R_score"]) if best_control is not None else float("nan"),
                "best_coarse_R": float(best_coarse["R"]) if best_coarse is not None else float("nan"),
                "best_coarse_joint_score": float(best_coarse["joint_R_score"]) if best_coarse is not None else float("nan"),
                "best_fine_R": float(best_fine["R"]) if best_fine is not None else float("nan"),
                "best_fine_joint_score": float(best_fine["joint_R_score"]) if best_fine is not None else float("nan"),
            }
        )
    return out


def refined_r_values_from_summaries(
    summaries: List[Dict[str, Any]],
    base_values: List[float],
    half_width: float,
    step: float,
    top_k: int,
    fixed_l_norm: Optional[float] = None,
    fixed_n_norm: Optional[float] = None,
) -> List[float]:
    if half_width <= 0.0:
        raise ValueError("refine half width must be positive")
    if step <= 0.0:
        raise ValueError("refine R step must be positive")
    if top_k <= 0:
        raise ValueError("refine top-k must be positive")
    if not base_values:
        raise ValueError("base R values must not be empty")

    add_pair_level_scores(summaries, fixed_l_norm=fixed_l_norm, fixed_n_norm=fixed_n_norm)
    r_min = min(base_values)
    r_max = max(base_values)
    values = set(round(float(value), 12) for value in base_values)
    by_case: Dict[str, List[Dict[str, Any]]] = {}
    for row in summaries:
        if not math.isfinite(float(row.get("joint_R_score", float("nan")))):
            continue
        if str(row.get("sweep_region")) == "control":
            continue
        by_case.setdefault(str(row["case_id"]), []).append(row)

    for rows in by_case.values():
        ordered = sorted(rows, key=lambda row: float(row["R"]))
        local_minima: List[Dict[str, Any]] = []
        for i, row in enumerate(ordered):
            score = float(row["joint_R_score"])
            left = float(ordered[i - 1]["joint_R_score"]) if i > 0 else float("inf")
            right = float(ordered[i + 1]["joint_R_score"]) if i + 1 < len(ordered) else float("inf")
            if score <= left and score <= right:
                local_minima.append(row)
        candidates = local_minima or ordered
        selected = sorted(candidates, key=lambda row: float(row["joint_R_score"]))[:top_k]
        for row in selected:
            center = float(row["R"])
            start = max(r_min, center - half_width)
            stop = min(r_max, center + half_width)
            values.update(round(float(value), 12) for value in uniform_r_values(start, stop, step))

    return sorted(float(value) for value in values)


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def make_score_plot(best: List[Dict[str, Any]], summaries: List[Dict[str, Any]], file_stem: str) -> str:
    cols = 2
    rows_n = math.ceil(len(best) / cols)
    fig, axes = plt.subplots(rows_n, cols, figsize=(12, 3.2 * rows_n), sharex=True)
    axes_flat = np.atleast_1d(axes).reshape(-1)
    for ax, best_row in zip(axes_flat, best):
        case_id = str(best_row["case_id"])
        pair_rows = sorted(
            [row for row in summaries if str(row["case_id"]) == case_id],
            key=lambda row: float(row["R"]),
        )
        ax.plot([row["R"] for row in pair_rows], [row["joint_R_score"] for row in pair_rows], marker="o", linewidth=1.2)
        ax.axvline(best_row["R_star_joint"], color="black", linestyle="--", linewidth=1)
        ax.axvline(Params.r137, color="tab:red", linestyle=":", linewidth=1)
        ax.axvline(Params.r128, color="tab:purple", linestyle=":", linewidth=1)
        ax.set_title(f"{best_row['mode']} B=[{best_row['harmonic_packet_B']}] R*={best_row['R_star_joint']:.6f}")
        ax.set_ylabel("joint score")
        ax.grid(alpha=0.25)
    for ax in axes_flat[len(best):]:
        ax.axis("off")
    for ax in axes_flat[-cols:]:
        ax.set_xlabel("R")
    fig.tight_layout()
    filename = f"{file_stem}_scores_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=160)
    plt.close(fig)
    return filename


def terrain_rows(
    all_rows: List[Dict[str, Any]],
    params: Params,
    fixed_l_norm: Optional[float] = None,
    fixed_n_norm: Optional[float] = None,
) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, float, int], Dict[str, Dict[str, Any]]] = {}
    for row in all_rows:
        key = (str(row["case_id"]), float(row["R_input"]), int(row["collision"]))
        grouped.setdefault(key, {})[str(row["channel"])] = row

    out: List[Dict[str, Any]] = []
    for key in sorted(grouped):
        _case_id, r_value, collision = key
        by_channel = grouped[key]
        if "A_channel" not in by_channel or "B_channel" not in by_channel:
            continue
        a = by_channel["A_channel"]
        b = by_channel["B_channel"]
        n_a = int(a["N_A"])
        n_b = int(a["N_B"])
        if params.fine_r_min <= r_value <= params.fine_r_max:
            sweep_region = "fine"
        elif params.coarse_r_min <= r_value <= params.coarse_r_max:
            sweep_region = "coarse"
        else:
            sweep_region = "control"
        out.append(
            {
                "case_id": str(a["case_id"]),
                "mode": str(a["mode"]),
                "state_family": str(a["state_family"]),
                "harmonic_packet_A": str(a["harmonic_packet_A"]),
                "harmonic_packet_B": str(a["harmonic_packet_B"]),
                "harmonic_packet_A_weights": str(a["harmonic_packet_A_weights"]),
                "harmonic_packet_B_weights": str(a["harmonic_packet_B_weights"]),
                "harmonic_packet_A_phases": str(a["harmonic_packet_A_phases"]),
                "harmonic_packet_B_phases": str(a["harmonic_packet_B_phases"]),
                "harmonic_packet_A_wavelength_scales": str(a["harmonic_packet_A_wavelength_scales"]),
                "harmonic_packet_B_wavelength_scales": str(a["harmonic_packet_B_wavelength_scales"]),
                "harmonic_packet_A_shift": float(a["harmonic_packet_A_shift"]),
                "harmonic_packet_B_shift": float(a["harmonic_packet_B_shift"]),
                "N_A": n_a,
                "N_B": n_b,
                "pair_kind": "one_side_high_harmonic" if n_a != n_b else "same_harmonic_control",
                "R": r_value,
                "T": float(a["T"]),
                "Delta_F": float(a["Delta_F"]),
                "collision": collision,
                "sweep_region": sweep_region,
                "L_A": float(a["L"]),
                "L_B": float(b["L"]),
                "L_gap": abs(float(a["L"]) - float(b["L"])),
                "N_eff_A": float(a["N_eff"]),
                "N_eff_B": float(b["N_eff"]),
                "N_eff_gap": abs(float(a["N_eff"]) - float(b["N_eff"])),
                "B_to_A_transfer": float(a["sim_to_B0"]),
                "A_to_B_transfer": float(b["sim_to_A0"]),
            }
        )

    by_pair: Dict[str, List[Dict[str, Any]]] = {}
    for row in out:
        by_pair.setdefault(str(row["case_id"]), []).append(row)
    for rows in by_pair.values():
        max_l = float(fixed_l_norm) if fixed_l_norm is not None else max([float(row["L_gap"]) for row in rows] or [1.0])
        max_n = float(fixed_n_norm) if fixed_n_norm is not None else max([float(row["N_eff_gap"]) for row in rows] or [1.0])
        max_l = max(max_l, 1.0e-300)
        max_n = max(max_n, 1.0e-300)
        for row in rows:
            gap_score = (float(row["L_gap"]) / max_l) + (float(row["N_eff_gap"]) / max_n)
            score = (float(row["L_gap"]) / max_l) + (float(row["N_eff_gap"]) / max_n) + (1.0 - float(row["B_to_A_transfer"]))
            row["gap_terrain_score"] = gap_score
            row["log10_gap_terrain_score"] = math.log10(max(gap_score, 1.0e-300))
            row["gap_depth"] = -row["log10_gap_terrain_score"]
            row["gap_depth_plot"] = min(row["gap_depth"], 6.0)
            row["joint_terrain_score"] = score
            row["log10_joint_terrain_score"] = math.log10(max(score, 1.0e-300))
            row["joint_depth"] = -row["log10_joint_terrain_score"]
    return out


def primary_case_ids(rows: List[Dict[str, Any]]) -> List[str]:
    case_ids = {
        str(row["case_id"])
        for row in rows
        if int(row["N_A"]) == 1 and str(row["sweep_region"]) in {"coarse", "fine"}
    }
    return sorted(case_ids)


def case_label_from_rows(rows: List[Dict[str, Any]], case_id: str) -> str:
    row = next(item for item in rows if str(item["case_id"]) == case_id)
    return f"{row['mode']} A=[{row['harmonic_packet_A']}] B=[{row['harmonic_packet_B']}]"


def best_depth_rows_by_r(rows: List[Dict[str, Any]], case_id: str) -> List[Dict[str, Any]]:
    by_r: Dict[float, List[Dict[str, Any]]] = {}
    for row in rows:
        if str(row["case_id"]) != case_id:
            continue
        if str(row["sweep_region"]) not in {"coarse", "fine"}:
            continue
        by_r.setdefault(float(row["R"]), []).append(row)
    return [max(r_rows, key=lambda row: float(row["gap_depth"])) for _, r_rows in sorted(by_r.items())]


def make_gap_depth_distribution_overview_plot(rows: List[Dict[str, Any]], params: Params, file_stem: str) -> str:
    case_ids = primary_case_ids(rows)
    fig, axes = plt.subplots(len(case_ids), 1, figsize=(12, 2.6 * len(case_ids)), sharex=True, sharey=True, constrained_layout=True)
    axes_list = list(np.atleast_1d(axes))
    for ax, case_id in zip(axes_list, case_ids):
        pair_rows = [
            row
            for row in rows
            if str(row["case_id"]) == case_id and str(row["sweep_region"]) in {"coarse", "fine"}
        ]
        ax.scatter(
            [float(row["R"]) for row in pair_rows],
            [float(row["gap_depth_plot"]) for row in pair_rows],
            s=8,
            color="0.70",
            alpha=0.18,
            linewidths=0,
            label="all collisions",
        )
        envelope = best_depth_rows_by_r(rows, case_id)
        ax.plot(
            [float(row["R"]) for row in envelope],
            [float(row["gap_depth_plot"]) for row in envelope],
            color="black",
            linewidth=1.4,
            label="deepest collision at each R",
        )
        best_row = max(pair_rows, key=lambda row: float(row["gap_depth"]))
        ax.scatter([float(best_row["R"])], [float(best_row["gap_depth_plot"])], marker="x", s=90, color="black", linewidths=1.8)
        ax.axvline(params.r137, color="black", linestyle=":", linewidth=1.1)
        ax.axvline(params.r128, color="0.35", linestyle="--", linewidth=1.1)
        ax.set_title(f"{case_label_from_rows(rows, case_id)} R-depth distribution")
        ax.set_ylabel("depth = -log10(gap score), capped at 6")
        ax.grid(alpha=0.20)
    axes_list[-1].set_xlabel("R")
    axes_list[0].legend(loc="upper right", fontsize=8)
    filename = f"{file_stem}_gap_depth_distribution_overview_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=170)
    plt.close(fig)
    return filename


def make_gap_depth_distribution_deep_plot(rows: List[Dict[str, Any]], params: Params, file_stem: str) -> str:
    case_ids = primary_case_ids(rows)
    fig, axes = plt.subplots(len(case_ids), 1, figsize=(12, 2.6 * len(case_ids)), sharex=False, sharey=True, constrained_layout=True)
    axes_list = list(np.atleast_1d(axes))
    for ax, case_id in zip(axes_list, case_ids):
        pair_rows = [
            row
            for row in rows
            if str(row["case_id"]) == case_id and str(row["sweep_region"]) in {"coarse", "fine"}
        ]
        best_row = max(pair_rows, key=lambda row: float(row["gap_depth"]))
        r_center = float(best_row["R"])
        focused = [row for row in pair_rows if abs(float(row["R"]) - r_center) <= 0.020]
        ax.scatter(
            [float(row["R"]) for row in focused],
            [float(row["gap_depth_plot"]) for row in focused],
            s=12,
            color="0.68",
            alpha=0.22,
            linewidths=0,
        )
        envelope = [
            row
            for row in best_depth_rows_by_r(rows, case_id)
            if abs(float(row["R"]) - r_center) <= 0.020
        ]
        ax.plot([float(row["R"]) for row in envelope], [float(row["gap_depth_plot"]) for row in envelope], color="black", linewidth=1.6)
        top_rows = sorted(focused, key=lambda row: float(row["gap_depth"]), reverse=True)[:8]
        ax.scatter(
            [float(row["R"]) for row in top_rows],
            [float(row["gap_depth_plot"]) for row in top_rows],
            s=32,
            facecolors="white",
            edgecolors="black",
            linewidths=0.9,
        )
        ax.axvline(params.r137, color="black", linestyle=":", linewidth=1.1)
        ax.axvline(params.r128, color="0.35", linestyle="--", linewidth=1.1)
        ax.scatter([r_center], [float(best_row["gap_depth_plot"])], marker="x", s=95, color="black", linewidths=1.8)
        ax.set_title(f"{case_label_from_rows(rows, case_id)} deepest R-depth, R={r_center:.12g}, collision={int(best_row['collision'])}")
        ax.set_xlabel("R")
        ax.set_ylabel("depth")
        ax.grid(alpha=0.20)
    filename = f"{file_stem}_gap_depth_distribution_deep_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=170)
    plt.close(fig)
    return filename


def build_report(params: Params, best: List[Dict[str, Any]], outputs: Dict[str, str]) -> str:
    lines: List[str] = [
        "# 系統A 局在性交換R近傍斉一スイープ 予備実験結果 v1",
        "",
        "## 1. 実験条件",
        "",
        "```text",
        f"coarse R range = {params.coarse_r_min} .. {params.coarse_r_max}",
        f"coarse Delta R = {params.coarse_r_step}",
        f"fine R range = {params.fine_r_min} .. {params.fine_r_max}",
        f"fine Delta R = {params.fine_r_step}",
        f"R_137 = {params.r137}",
        f"R_128 = {params.r128}",
        f"max_collision = {params.recursive_collision_count}",
        f"pairs = {params.pairs}",
        f"cases = {[row['case_id'] for row in best]}",
        "```",
        "",
        "## 2. 系統A 判定サマリー",
        "",
        "| mode | A packet | B packet | A weights | B weights | A phases | B phases | A wavelengths | B wavelengths | A shift | B shift | N_A | N_B | kind | R*_L | col_L | R*_N | col_N | R*_transfer | col_transfer | R*_joint | col_joint | band5 | band10 | d137 | d128 | coarse R | fine R | control R |",
        "|---|---|---|---|---|---|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in best:
        lines.append(
            "| {mode} | {harmonic_packet_A} | {harmonic_packet_B} | {harmonic_packet_A_weights} | {harmonic_packet_B_weights} | "
            "{harmonic_packet_A_phases} | {harmonic_packet_B_phases} | {harmonic_packet_A_wavelength_scales} | {harmonic_packet_B_wavelength_scales} | "
            "{harmonic_packet_A_shift:.12g} | {harmonic_packet_B_shift:.12g} | "
            "{N_A} | {N_B} | {pair_kind} | {R_star_L:.12g} | {collision_at_R_star_L} | "
            "{R_star_N:.12g} | {collision_at_R_star_N} | {R_star_transfer:.12g} | {collision_at_R_star_transfer} | "
            "{R_star_joint:.12g} | {collision_at_R_star_joint} | {R_band_width_5:.12g} | {R_band_width_10:.12g} | "
            "{distance_to_R_137:.12g} | {distance_to_R_128:.12g} | {best_coarse_R:.12g} | {best_fine_R:.12g} | {best_control_R:.12g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## 3. 読み",
            "",
            "本予備実験では、片側高次倍音条件と同次数対照を同じ `R` 掃引で比較した。",
            "",
            "`R_star_joint` は、広域粗スイープと局所精密スイープを合わせた候補範囲の中で、局在性差、実効次数差、B側初期倍音分布のA側への類似度を合わせて読んだ値である。",
            "",
            "`R=0.0`, `R=0.5`, `R=1.0` は対照であり、`R_star` 判定からは除外した。",
            "",
            "同次数対照では `R_band_width` が主掃引幅全体に広がり、`R` に対して鋭い判定点を持たない。",
            "",
            "一方、片側高次倍音条件では、`R_star_L`, `R_star_N`, `R_star_transfer`, `R_star_joint` が同じ位置に集まった。",
            "",
            "ただし、広域粗スイープでは複数の浅い谷候補も見える。",
            "",
            "特に片側高次倍音条件では、粗スイープ上で `R=0.74` 付近が最良の粗候補となり、`R=0.63`, `R=0.67`, `R=0.78`, `R=0.86`, `R=0.89` 付近にも浅い谷が現れた。",
            "",
            "したがって、`R_137` 近傍だけを細分化するのではなく、これらの谷候補も次段階の局所精密スイープ対象に含める。",
            "",
            "`R_137` と `R_128` は判定値ではなく、全系統で同じ位置を読むための固定プローブ点である。",
            "",
            "本予備実験では `R_137` を明示的にプローブ点として含めたため、次段階では `R_137` 近傍を特別扱いせず、より細かい一様掃引で谷幅を確認する必要がある。",
            "",
            "また、点としての `R_star` だけでは、反射係数の最小点と衝突回数方向の干渉縞を分離できない。",
            "",
            "そのため、本実装では `R` ごとの全衝突回を縦方向へ並べる分布図を追加した。",
            "",
            "分布図では、`L_gap` と `N_eff_gap` から作った `gap_terrain_score` を `depth=-log10(gap_terrain_score)` として縦軸に置く。",
            "",
            "完全同一対照では数値的に深さが過大になるため、図の表示では `depth` を 6 で上限表示する。CSV には上限前の `gap_depth` も保存する。",
            "",
            "灰色点は各 `R` における全衝突回の分布、黒線は各 `R` で最も深い衝突回の包絡線である。",
            "",
            "最深部の拡大分布図は、次段階で局所スイープ幅を決めるための候補図である。",
            "",
            "## 4. 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
        ]
    )
    for label, filename in outputs.items():
        lines.append(f"| {label} | `{filename}` |")
    lines.append("")
    return "\n".join(lines)


def run(
    selected_cases: Optional[Tuple[HarmonicCase, ...]] = None,
    selected_r_values: Optional[List[float]] = None,
    output_dir: Optional[Path] = None,
    max_collision: Optional[int] = None,
    make_plots: bool = True,
    run_id: Optional[str] = None,
    file_stem: Optional[str] = None,
    fixed_l_norm: Optional[float] = None,
    fixed_n_norm: Optional[float] = None,
) -> Dict[str, Any]:
    params = Params()
    if max_collision is not None:
        params.recursive_collision_count = max_collision
    cases = selected_cases if selected_cases is not None else default_cases(params)
    if not cases:
        raise ValueError("no harmonic case selected")

    r_list = selected_r_values if selected_r_values is not None else r_values(params)
    if not r_list:
        raise ValueError("no R value selected")
    if output_dir is not None:
        configure_output_dir(output_dir)
    resolved_file_stem = build_file_stem(cases, r_list, params, run_id, file_stem)

    source_params = build_source_params(params)
    metrics = MetricContext(source_params)
    all_rows: List[Dict[str, Any]] = []
    summaries: List[Dict[str, Any]] = []
    for case in cases:
        for r_value in r_list:
            rows = run_case(source_params, metrics, case, r_value, params.recursive_collision_count)
            all_rows.extend(rows)
            summaries.append(summarize_case(rows, params))
    best = best_rows_for_pair(summaries, params, fixed_l_norm=fixed_l_norm, fixed_n_norm=fixed_n_norm)
    outputs = {
        "rows": f"{resolved_file_stem}_rows_v1.csv",
        "summary": f"{resolved_file_stem}_summary_v1.csv",
        "best": f"{resolved_file_stem}_best_v1.csv",
        "terrain": f"{resolved_file_stem}_collision_terrain_v1.csv",
        "json": f"{resolved_file_stem}_result_v1.json",
        "report": f"{resolved_file_stem}_report_v1.md",
    }
    if make_plots:
        outputs["score_plot"] = make_score_plot(best, summaries, resolved_file_stem)
    terrain = terrain_rows(all_rows, params, fixed_l_norm=fixed_l_norm, fixed_n_norm=fixed_n_norm)
    if make_plots and primary_case_ids(terrain):
        outputs["depth_distribution_overview_plot"] = make_gap_depth_distribution_overview_plot(terrain, params, resolved_file_stem)
        outputs["depth_distribution_deep_plot"] = make_gap_depth_distribution_deep_plot(terrain, params, resolved_file_stem)
    write_csv(OUT_DIR / outputs["rows"], all_rows)
    write_csv(OUT_DIR / outputs["summary"], summaries)
    write_csv(OUT_DIR / outputs["best"], best)
    write_csv(OUT_DIR / outputs["terrain"], terrain)
    result: Dict[str, Any] = {
        "experiment": "system_A_localization_exchange_R_sweep_preliminary_v1",
        "params": asdict(params),
        "cases": [asdict(case) for case in cases],
        "file_stem": resolved_file_stem,
        "run_id": run_id,
        "output_dir": str(OUT_DIR),
        "source_script": str(SOURCE_PATH.relative_to(BASE_DIR.parent)),
        "r_values": r_list,
        "best": best,
        "summary": summaries,
        "outputs": outputs,
        "fixed_l_norm": fixed_l_norm,
        "fixed_n_norm": fixed_n_norm,
    }
    (OUT_DIR / outputs["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    report = build_report(params, best, outputs)
    (OUT_DIR / outputs["report"]).write_text(report, encoding="utf-8")
    return result


if __name__ == "__main__":
    cli_args = parse_args()
    cli_params = Params()
    selected_cases = selected_cases_from_args(cli_params, cli_args)
    selected_r = selected_r_values_from_args(cli_params, cli_args)
    selected_output_dir = output_dir_from_args(cli_args)
    fixed_l_norm = cli_args.fixed_l_norm
    fixed_n_norm = cli_args.fixed_n_norm
    if cli_args.refine_minima:
        initial_file_stem = f"{cli_args.file_stem}_initial" if cli_args.file_stem else None
        initial_data = run(
            selected_cases=selected_cases,
            selected_r_values=selected_r,
            output_dir=selected_output_dir,
            max_collision=cli_args.max_collision,
            make_plots=not cli_args.no_plots,
            run_id=cli_args.run_id,
            file_stem=initial_file_stem,
            fixed_l_norm=fixed_l_norm,
            fixed_n_norm=fixed_n_norm,
        )
        refine_half_width = cli_args.refine_half_width if cli_args.refine_half_width is not None else cli_params.refine_half_width
        refine_r_step = cli_args.refine_r_step if cli_args.refine_r_step is not None else cli_params.refine_r_step
        refine_top_k = cli_args.refine_top_k if cli_args.refine_top_k is not None else cli_params.refine_top_k
        refined_r = refined_r_values_from_summaries(
            initial_data["summary"],
            selected_r,
            half_width=float(refine_half_width),
            step=float(refine_r_step),
            top_k=int(refine_top_k),
            fixed_l_norm=fixed_l_norm,
            fixed_n_norm=fixed_n_norm,
        )
        data = run(
            selected_cases=selected_cases,
            selected_r_values=refined_r,
            output_dir=selected_output_dir,
            max_collision=cli_args.max_collision,
            make_plots=not cli_args.no_plots,
            run_id=cli_args.run_id,
            file_stem=cli_args.file_stem,
            fixed_l_norm=fixed_l_norm,
            fixed_n_norm=fixed_n_norm,
        )
        print(
            json.dumps(
                {
                    "output_dir": data["output_dir"],
                    "initial_file_stem": initial_data["file_stem"],
                    "file_stem": data["file_stem"],
                    "initial_r_count": len(selected_r),
                    "refined_r_count": len(refined_r),
                    "refine_half_width": float(refine_half_width),
                    "refine_r_step": float(refine_r_step),
                    "refine_top_k": int(refine_top_k),
                    "fixed_l_norm": fixed_l_norm,
                    "fixed_n_norm": fixed_n_norm,
                    "best": data["best"],
                    "outputs": data["outputs"],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
    else:
        data = run(
            selected_cases=selected_cases,
            selected_r_values=selected_r,
            output_dir=selected_output_dir,
            max_collision=cli_args.max_collision,
            make_plots=not cli_args.no_plots,
            run_id=cli_args.run_id,
            file_stem=cli_args.file_stem,
            fixed_l_norm=fixed_l_norm,
            fixed_n_norm=fixed_n_norm,
        )
        print(
            json.dumps(
                {
                    "output_dir": data["output_dir"],
                    "file_stem": data["file_stem"],
                    "best": data["best"],
                    "outputs": data["outputs"],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
