#!/usr/bin/env python3
"""上流GENESIS結果の閉包・倍音構造頑健性監査。

事前登録は同じフォルダの
`事前登録_上流インフレーション倍音頑健性監査_v1.md` を参照。
親コードはread-only importし、既存成果物を上書きしない。
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
S3_PATH = PARENT / "run_stage3_sharedO_v2_and_hair_v1.py"
GENESIS_PATH = PARENT / "run_genesis_v1.py"
PUBLISHED_RESULT = PARENT / "genesis_result_v1.json"

spec = importlib.util.spec_from_file_location("upstream_s3", S3_PATH)
s3 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = s3
assert spec.loader is not None
spec.loader.exec_module(s3)

abl = s3.abl
gen3 = s3.gen3
V2 = s3.VertexEngineV2

N = 5
M = 10
T_LONG = 4000
DELTAS = (1e-6, 1e-2)
WIN_LAT = (50, 250)
WIN_BURST = (400, 1100)
WIN_META = (2000, 4000)
ODD63 = tuple(range(1, 64, 2))
EVEN63 = tuple(range(2, 63, 2))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bilinear(x: np.ndarray, y: np.ndarray) -> complex:
    """複素双線形積。Hermite内積ではない。"""
    return complex(x @ y)


def make_closed_seed(z: np.ndarray, old_seed: np.ndarray) -> np.ndarray:
    """S^T S=Z^T S=0 を満たす単位種を決定論的に構成する。"""
    basis: list[np.ndarray] = []

    def add_real_direction(candidate: np.ndarray) -> None:
        v = np.asarray(candidate, dtype=float).copy()
        for q in basis:
            v -= q * float(q @ v)
        nv = float(np.linalg.norm(v))
        if nv > 1e-12:
            basis.append(v / nv)

    # 最初の2本がポンプの実平面。その後の2本を種に使う。
    add_real_direction(z.real)
    add_real_direction(z.imag)
    pump_rank = len(basis)
    if pump_rank != 2:
        raise RuntimeError(f"pump real-plane rank must be 2, got {pump_rank}")
    add_real_direction(old_seed.real)
    add_real_direction(old_seed.imag)
    for e in np.eye(z.size):
        if len(basis) >= 4:
            break
        add_real_direction(e)
    if len(basis) < 4:
        raise RuntimeError("could not construct two transverse real directions")
    a, b = basis[2], basis[3]
    seed = (a + 1j * b) / np.sqrt(2.0)
    # 旧種とのHermite重なりが正の実数になる位相を選ぶ。
    overlap = np.vdot(old_seed, seed)
    if abs(overlap) > 0:
        seed *= np.exp(-1j * np.angle(overlap))
    return seed / np.linalg.norm(seed)


class FastVertexEngineV2(V2):
    """V2の共有Cayley写像を列ごとの反復から行列積へ置換した等価実装。"""

    def _linear(self) -> None:
        zsum = np.sum(self.C, axis=1)
        self.sys_shared.set_theta(np.angle(zsum))
        se, self.wp = self.sys_shared.sigma_max_power(self.wp)
        eye = np.eye(self.m)
        transform = np.column_stack(
            [self.sys_shared.cayley_step(eye[:, j], se) for j in range(self.m)]
        )
        self.C = transform @ self.C


def fval_factory(z0: np.ndarray):
    p = z0.real / np.linalg.norm(z0.real)
    q = z0.imag - (z0.imag @ p) * p
    q = q / np.linalg.norm(q)

    def fval(z: np.ndarray) -> float:
        denom = float(np.real(np.vdot(z, z)))
        if denom <= 0:
            return float("nan")
        zp = z - p * (p @ z) - q * (q @ z)
        return float(np.real(np.vdot(zp, zp))) / denom

    return fval


def window_rate(ts: np.ndarray, ln_p: np.ndarray, lo: int, hi: int) -> float:
    mask = (ts >= lo) & (ts < hi)
    a = np.vstack([ts[mask].astype(float), np.ones(np.count_nonzero(mask))]).T
    coef, _, _, _ = np.linalg.lstsq(a, ln_p[mask], rcond=None)
    return float(coef[0])


def initial_state(layout: str, delta: float, z: np.ndarray,
                  old_seed: np.ndarray, closed_seed: np.ndarray):
    if layout == "L-old":
        nreg, evens, odds, seed = 5, (2,), (1,), old_seed
    elif layout == "L-cl":
        nreg, evens, odds, seed = 5, (2,), (1,), closed_seed
    elif layout == "O63-cl":
        nreg, evens, odds, seed = 512, (2,), ODD63, closed_seed
    elif layout == "EO63-cl":
        nreg, evens, odds, seed = 512, EVEN63, ODD63, closed_seed
    else:
        raise ValueError(layout)

    c0 = np.zeros((z.size, nreg), dtype=complex)
    for k in evens:
        c0[:, k] = z / np.sqrt(len(evens))
    for k in odds:
        c0[:, k] = delta * seed / np.sqrt(len(odds))
    return c0, evens, odds


def pointwise_closure(c: np.ndarray) -> float:
    w = np.fft.ifft(c, axis=1) * c.shape[1]
    return float(np.max(np.abs(np.sum(w ** 2, axis=0))))


def validate_fast_engine(c0: np.ndarray, wp0: np.ndarray, steps: int = 3) -> dict:
    slow = V2(N, c0, wp0, vertex_on=True)
    fast = FastVertexEngineV2(N, c0, wp0, vertex_on=True)
    max_rel = 0.0
    for _ in range(steps):
        slow.step()
        fast.step()
        rel = float(np.linalg.norm(slow.C - fast.C) /
                    max(np.linalg.norm(slow.C), 1e-300))
        max_rel = max(max_rel, rel)
    return {"steps": steps, "max_relative_C_difference": max_rel,
            "pass_le_1e-12": bool(max_rel <= 1e-12)}


def run_one(layout: str, delta: float, z: np.ndarray, old_seed: np.ndarray,
            closed_seed: np.ndarray, wp0: np.ndarray, steps: int,
            full_engine: str = "fast") -> tuple[dict, dict]:
    c0, evens, odds = initial_state(layout, delta, z, old_seed, closed_seed)
    engine_cls = (V2 if c0.shape[1] == 5 or full_engine == "slow"
                  else FastVertexEngineV2)
    eng = engine_cls(N, c0, wp0, vertex_on=True)
    fval = fval_factory(z)
    scale = np.sqrt(len(evens))
    initial_support = np.zeros(c0.shape[1], dtype=bool)
    initial_support[list(evens) + list(odds)] = True

    f_k2 = np.zeros(steps)
    f_even = np.zeros(steps)
    f_odd = np.zeros(steps)
    closure = np.zeros(steps + 1)
    norm = np.zeros(steps + 1)
    crossing_k2 = None
    crossing_even = None
    closure[0] = pointwise_closure(eng.C)
    norm[0] = eng.diagnostics()["norm"]

    for t in range(steps):
        eng.step()
        z_k2 = eng.C[:, 2]
        z_even = np.sum(eng.C[:, evens], axis=1) / scale
        f_k2[t] = fval(z_k2)
        f_even[t] = fval(z_even)
        diag = eng.diagnostics()
        f_odd[t] = diag["f_seed"]
        closure[t + 1] = diag["closure_max"]
        norm[t + 1] = diag["norm"]
        if crossing_k2 is None and f_k2[t] > 0.05:
            crossing_k2 = t + 1
        if crossing_even is None and f_even[t] > 0.05:
            crossing_even = t + 1

    ts = np.arange(1, steps + 1)
    ln_odd = np.log(np.maximum(f_odd, 1e-300))
    rates = {}
    for name, (lo, hi) in {
        "latency": WIN_LAT, "burst": WIN_BURST, "metastable": WIN_META
    }.items():
        rates[name] = (window_rate(ts, ln_odd, lo, min(hi, steps + 1))
                       if steps >= hi else None)

    p_final = np.sum(np.abs(eng.C) ** 2, axis=0)
    p_total = float(p_final.sum())
    outside = float(p_final[~initial_support].sum() / p_total) if p_total > 0 else 0.0
    odd_mask = np.arange(c0.shape[1]) % 2 == 1
    even_mask = (np.arange(c0.shape[1]) % 2 == 0)
    even_mask[0] = False
    if c0.shape[1] % 2 == 0:
        even_mask[c0.shape[1] // 2] = False
    outside_odd = (~initial_support) & odd_mask
    outside_even = (~initial_support) & even_mask
    summary = {
        "layout": layout,
        "delta": delta,
        "nreg": int(c0.shape[1]),
        "even_modes": list(evens),
        "odd_modes": list(odds),
        "initial_total_power": float(np.sum(np.abs(c0) ** 2)),
        "initial_even_power": float(np.sum(np.abs(c0[:, evens]) ** 2)),
        "initial_odd_power": float(np.sum(np.abs(c0[:, odds]) ** 2)),
        "initial_closure_max": float(closure[0]),
        "closure_max_over_run": float(closure.max()),
        "closure_change_max": float(np.max(np.abs(closure - closure[0]))),
        "norm_relative_drift_max": float(np.max(np.abs(norm - norm[0])) / norm[0]),
        "crossing_k2": crossing_k2,
        "crossing_even_aggregate": crossing_even,
        "rates_dlog_fodd": rates,
        "f_odd_first": float(f_odd[0]),
        "f_odd_final": float(f_odd[-1]),
        "final_power_outside_initial_modes": outside,
        "final_odd_power_outside_initial_modes": (
            float(p_final[outside_odd].sum() / p_total) if p_total > 0 else 0.0),
        "final_even_power_outside_initial_modes": (
            float(p_final[outside_even].sum() / p_total) if p_total > 0 else 0.0),
        "final_dc_power_fraction": float(p_final[0] / p_total) if p_total > 0 else 0.0,
        "final_nyquist_power_fraction": (
            float(p_final[c0.shape[1] // 2] / p_total)
            if p_total > 0 and c0.shape[1] % 2 == 0 else None),
    }
    series = {"f_k2": f_k2, "f_even": f_even, "f_odd": f_odd,
              "closure": closure, "norm": norm, "final_mode_power": p_final}
    return summary, series


def relative_crossing_difference(a, b) -> float | None:
    if a is None or b is None:
        return None
    return abs(float(a) - float(b)) / max(abs(float(b)), 1.0)


def compare_to_published(baseline: dict, published: dict) -> dict:
    old = published["main"] if baseline["delta"] == 1e-6 else published["backreaction"]
    diffs = {
        "crossing_equal": baseline["crossing_k2"] == old["crossing"],
        "g_latency_abs_diff": abs(baseline["rates_dlog_fodd"]["latency"] - old["g_latency"]),
        "g_burst_abs_diff": abs(baseline["rates_dlog_fodd"]["burst"] - old["g_burst"]),
        "g_metastable_abs_diff": abs(baseline["rates_dlog_fodd"]["metastable"] - old["g_metastable"]),
        "f_seed0_abs_diff": abs(baseline["f_odd_first"] - old["f_seed0"]),
        "f_seed_final_abs_diff": abs(baseline["f_odd_final"] - old["f_seed_final"]),
    }
    diffs["pass"] = bool(diffs["crossing_equal"] and
                         max(v for k, v in diffs.items() if k.endswith("abs_diff")) <= 1e-15)
    return diffs


def sensitivity(reference: dict, candidate: dict) -> dict:
    c_ref = reference["crossing_even_aggregate"]
    c_new = candidate["crossing_even_aggregate"]
    cross_rel = relative_crossing_difference(c_new, c_ref)
    crossing_sensitive = ((c_ref is None) != (c_new is None) or
                          (cross_rel is not None and cross_rel > 0.05))
    rate_flags = {}
    for window in ("latency", "burst", "metastable"):
        x = reference["rates_dlog_fodd"][window]
        y = candidate["rates_dlog_fodd"][window]
        if x is None or y is None:
            rate_flags[window] = None
            continue
        sign_change = np.sign(x) != np.sign(y)
        ax, ay = abs(x), abs(y)
        decade_change = max(ax, ay) > 10.0 * max(min(ax, ay), 1e-300)
        rate_flags[window] = bool(sign_change or decade_change)
    return {
        "crossing_relative_difference": cross_rel,
        "crossing_sensitive_gt_5pct_or_missing": bool(crossing_sensitive),
        "rate_sign_or_decade_change": rate_flags,
        "representation_sensitive": bool(crossing_sensitive or
            any(v is True for v in rate_flags.values())),
    }


def make_figure(all_series: dict, delta: float, steps: int, suffix: str) -> Path:
    ts = np.arange(1, steps + 1)
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    colors = {"L-old": "#7f7f7f", "L-cl": "#4c78a8",
              "O63-cl": "#f58518", "EO63-cl": "#54a24b"}
    for layout, color in colors.items():
        s = all_series[(layout, delta)]
        axes[0].semilogy(ts, np.clip(s["f_even"], 1e-34, None),
                        label=layout, color=color, lw=1.0)
        axes[1].semilogy(ts, np.clip(s["f_odd"], 1e-300, None),
                        label=layout, color=color, lw=1.0)
        axes[2].semilogy(np.arange(steps + 1), np.clip(s["closure"], 1e-34, None),
                        label=layout, color=color, lw=1.0)
    axes[0].axhline(0.05, color="red", ls=":", lw=0.8)
    axes[0].set_ylabel("f outside pump plane")
    axes[1].set_ylabel("odd power fraction")
    axes[2].set_ylabel("max pointwise closure")
    axes[2].set_xlabel("step")
    axes[0].set_title(f"Upstream harmonic robustness, delta={delta:g}")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.18)
    fig.tight_layout()
    path = HERE / f"fig_upstream_harmonic_robustness_delta{delta:g}_{suffix}_v1.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true", help="20 steps, smoke artifacts")
    parser.add_argument("--full-engine", choices=("fast", "slow"), default="fast",
                        help="512点条件の共有線形部。slowは親V2を無変更使用")
    args = parser.parse_args()
    steps = 20 if args.smoke else T_LONG
    suffix = "smoke" if args.smoke else "T4000"
    if args.full_engine == "slow":
        suffix += "_parentV2"
    started = time.time()

    _, _, _, _, _, _, _, z, wp0 = abl.build_init(N, False)
    parent = gen3.make_parent(N, seed=2)
    csec = np.fft.fft(parent.relation_waves, axis=1) / N
    old_seed = csec[:, 1] / np.linalg.norm(csec[:, 1])
    closed_seed = make_closed_seed(z, old_seed)

    construction = {
        "ZtZ_abs": abs(bilinear(z, z)),
        "old_StS_abs": abs(bilinear(old_seed, old_seed)),
        "old_ZtS_abs": abs(bilinear(z, old_seed)),
        "closed_StS_abs": abs(bilinear(closed_seed, closed_seed)),
        "closed_ZtS_abs": abs(bilinear(z, closed_seed)),
        "closed_seed_norm": float(np.linalg.norm(closed_seed)),
        "hermitian_overlap_old_closed_abs": float(abs(np.vdot(old_seed, closed_seed))),
    }

    c_full, _, _ = initial_state("EO63-cl", 1e-2, z, old_seed, closed_seed)
    fast_validation = validate_fast_engine(c_full, wp0, steps=3)
    if args.full_engine == "fast" and not fast_validation["pass_le_1e-12"]:
        raise RuntimeError(f"fast engine validation failed: {fast_validation}")

    summaries: dict[tuple[str, float], dict] = {}
    series: dict[tuple[str, float], dict] = {}
    for delta in DELTAS:
        for layout in ("L-old", "L-cl", "O63-cl", "EO63-cl"):
            print(f"run {layout:8s} delta={delta:g} T={steps}", flush=True)
            summary, trace = run_one(layout, delta, z, old_seed, closed_seed, wp0,
                                     steps, full_engine=args.full_engine)
            summaries[(layout, delta)] = summary
            series[(layout, delta)] = trace
            print(f"  crossing={summary['crossing_even_aggregate']} "
                  f"closure0={summary['initial_closure_max']:.3e} "
                  f"closure_max={summary['closure_max_over_run']:.3e}", flush=True)

    published = json.loads(PUBLISHED_RESULT.read_text(encoding="utf-8"))
    reproduction = {
        str(delta): compare_to_published(summaries[("L-old", delta)], published)
        for delta in DELTAS
    } if not args.smoke else None
    comparisons = {}
    for delta in DELTAS:
        comparisons[str(delta)] = {
            "closure_correction_Lcl_vs_Lold": sensitivity(
                summaries[("L-old", delta)], summaries[("L-cl", delta)]),
            "odd_packet_O63_vs_Lcl": sensitivity(
                summaries[("L-cl", delta)], summaries[("O63-cl", delta)]),
            "full_packet_EO63_vs_Lcl": sensitivity(
                summaries[("L-cl", delta)], summaries[("EO63-cl", delta)]),
        }

    output = {
        "status": "smoke" if args.smoke else "science",
        "preregistration": "事前登録_上流インフレーション倍音頑健性監査_v1.md",
        "parameters": {"N": N, "M": M, "T": steps, "deltas": list(DELTAS),
                       "odd_harmonics": list(ODD63), "even_harmonics": list(EVEN63),
                       "full_engine": args.full_engine},
        "source_sha256": {"this_script": sha256(Path(__file__)),
                          "stage3_engine": sha256(S3_PATH),
                          "genesis_v1": sha256(GENESIS_PATH),
                          "published_result": sha256(PUBLISHED_RESULT)},
        "closed_seed_construction": construction,
        "fast_engine_validation": fast_validation,
        "summaries": {f"{layout}|delta={delta:g}": value
                      for (layout, delta), value in summaries.items()},
        "published_reproduction": reproduction,
        "comparisons": comparisons,
        "runtime_sec": time.time() - started,
    }

    json_path = HERE / f"upstream_harmonic_robustness_result_{suffix}_v1.json"
    json_path.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    npz_payload = {}
    for (layout, delta), trace in series.items():
        safe = layout.replace("-", "_") + "_d" + f"{delta:g}".replace("-", "m")
        for name, values in trace.items():
            npz_payload[f"{safe}_{name}"] = values
    np.savez_compressed(HERE / f"upstream_harmonic_robustness_series_{suffix}_v1.npz",
                        **npz_payload)
    for delta in DELTAS:
        make_figure(series, delta, steps, suffix)
    print(f"saved {json_path.name}; runtime={output['runtime_sec']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
