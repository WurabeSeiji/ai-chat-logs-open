#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文A v2：倍音構造の対応誤りを調べる事前登録済み再実験。

既存の力学は一切変更せず、unified_interaction_v1.UnifiedEngine を
read-only import する。論文A v1 の単一モード16帯と、先行模型の
奇数倍音 1,3,...,63／偶数倍音 2,4,...,62 を同じ総パワーで比較する。
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
SERIES = HERE.parent.parent
F_PATH = SERIES / "統一万能関数_v1" / "unified_interaction_v1.py"
PREREG = HERE / "事前登録_倍音構造訂正再実験_v1.md"

N = 12
DELTA = 1.0e-2
SEED = 2
H_VALUES = (1, 3, 7, 15, 31, 63)
FULL_NK = 512
FULL_NETA = 16
SMALL_NK = 16
SMALL_NETA = 8
ETA_INDEX = 0
RATE_THRESHOLD = 1.0e-20
LOG_LINES: list[str] = []


def emit(message: str) -> None:
    print(message, flush=True)
    LOG_LINES.append(message)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"module load failed: {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


F = load_module("paperA_v2_unified_F1", F_PATH)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def jsonable(x):
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, complex):
        return {"real": x.real, "imag": x.imag}
    raise TypeError(type(x).__name__)


def parent_data(n: int):
    built = F.abl.build_init(n, False)
    Z0c = built[7]
    wp0 = built[8]
    parent = F.gen3.make_parent(n, seed=SEED)
    Csec = np.fft.fft(parent.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    return Z0c, wp0, seed_state


def harmonic_sets(H: int):
    odds = tuple(range(1, H + 1, 2))
    evens = tuple(range(2, H, 2)) if H >= 3 else (2,)
    return evens, odds


def make_packet(
    Z0c: np.ndarray,
    seed_state: np.ndarray,
    H: int,
    nk: int,
    neta: int,
) -> tuple[np.ndarray, tuple[int, ...], tuple[int, ...]]:
    evens, odds = harmonic_sets(H)
    if max(evens + odds) >= nk // 2:
        raise ValueError(f"H={H} cannot be represented without aliasing at Nk={nk}")
    C = np.zeros((Z0c.size, nk, neta), dtype=np.complex128)
    pump_scale = 1.0 / np.sqrt(len(evens))
    seed_scale = DELTA / np.sqrt(len(odds))
    for k in evens:
        C[:, k, ETA_INDEX] = pump_scale * Z0c
    for k in odds:
        C[:, k, ETA_INDEX] = seed_scale * seed_state
    return C, evens, odds


def fold_selected_k(
    C: np.ndarray,
    selected_modes: tuple[int, ...],
    out_nk: int = SMALL_NK,
    out_neta: int = SMALL_NETA,
) -> np.ndarray:
    """事前登録 Q：選んだ倍音部分空間を k mod 16 へ正規化線形射影。

    eta=0 の実験なので eta 軸は 0 をそのまま写す。未選択の、頂点が新たに
    生成した倍音は捨てる。この捨て方も縮約の非可換性に含めて報告する。
    """
    out = np.zeros((C.shape[0], out_nk, out_neta), dtype=np.complex128)
    groups: dict[int, list[int]] = {}
    for k in selected_modes:
        groups.setdefault(k % out_nk, []).append(k)
    for residue, modes in groups.items():
        out[:, residue, 0] = C[:, modes, 0].sum(axis=1) / np.sqrt(len(modes))
    return out


def signed_modes(nk: int) -> np.ndarray:
    return np.rint(np.fft.fftfreq(nk, d=1.0 / nk)).astype(int)


def parity_powers(C: np.ndarray) -> dict[str, float]:
    Pk = np.sum(np.abs(C) ** 2, axis=(0, 2))
    idx = np.arange(C.shape[1])
    return {
        "odd": float(Pk[idx % 2 == 1].sum()),
        "even_nonzero": float(Pk[(idx % 2 == 0) & (idx != 0)].sum()),
        "zero": float(Pk[0]),
        "total": float(Pk.sum()),
    }


def state_metrics(C: np.ndarray, wp0: np.ndarray) -> tuple[dict, dict[str, np.ndarray]]:
    nk, neta = C.shape[1], C.shape[2]
    eng = F.UnifiedEngine(N, C, wp0)
    R = eng._readout()
    W = np.fft.ifft2(C, axes=(1, 2)) * (nk * neta)
    closure = np.abs(np.sum(W ** 2, axis=0))
    chi_power = np.sum(np.abs(W) ** 2, axis=(0, 2))
    chi_prob = chi_power / chi_power.sum()
    pr = 1.0 / np.sum(chi_prob ** 2)

    rateW = eng._vertex_rate(W.reshape(C.shape[0], -1), R).reshape(C.shape)
    dC = np.fft.fft2(rateW, axes=(1, 2)) / (nk * neta)
    rate_spectrum = np.sum(np.abs(dC) ** 2, axis=(0, 2))
    freqs = signed_modes(nk)
    rate_total = float(rate_spectrum.sum())
    out_mask = np.abs(freqs) > 7
    rate_out = float(rate_spectrum[out_mask].sum())
    max_spec = float(rate_spectrum.max(initial=0.0))
    sig_mask = rate_spectrum > RATE_THRESHOLD * max(max_spec, np.finfo(float).tiny)
    max_abs_sig_k = int(np.max(np.abs(freqs[sig_mask]))) if np.any(sig_mask) else 0
    amp_max = float(np.max(np.abs(W)))
    rate_amp_max = float(np.max(np.abs(rateW)))
    estimated_nsub = max(
        1,
        int(np.ceil((rate_amp_max / max(amp_max, 1.0e-300)) / F.s2.H_MAX)),
    )
    powers = parity_powers(C)
    metrics = {
        "nk": nk,
        "neta": neta,
        "odd_power": powers["odd"],
        "even_nonzero_power": powers["even_nonzero"],
        "zero_power": powers["zero"],
        "total_power": powers["total"],
        "readout_mean": float(R.mean()),
        "readout_min": float(R.min()),
        "readout_max": float(R.max()),
        "chi_PR": float(pr),
        "chi_PR_fraction": float(pr / nk),
        "chi_peak_over_mean": float(chi_power.max() / chi_power.mean()),
        "closure_max": float(closure.max()),
        "rate_norm": float(np.sqrt(rate_total)),
        "rate_relative_norm": float(np.sqrt(rate_total) / np.linalg.norm(C)),
        "rate_power_outside_abs_k_le_7": rate_out / rate_total if rate_total > 0 else 0.0,
        "rate_max_abs_k_above_threshold": max_abs_sig_k,
        "estimated_RK4_substeps": estimated_nsub,
    }
    arrays = {
        "R": R,
        "chi_prob": chi_prob,
        "rate_spectrum": rate_spectrum,
        "freqs": freqs,
        "dC": dC,
    }
    return metrics, arrays


def one_step(C: np.ndarray, wp0: np.ndarray) -> tuple[np.ndarray, dict]:
    eng = F.UnifiedEngine(N, C, wp0)
    before = parity_powers(C)
    eng.step()
    after_C = eng.C2().copy()
    after = parity_powers(after_C)
    W = np.fft.ifft2(after_C, axes=(1, 2)) * (C.shape[1] * C.shape[2])
    closure = np.abs(np.sum(W ** 2, axis=0))
    return after_C, {
        "delta_C_relative": float(np.linalg.norm(after_C - C) / np.linalg.norm(C)),
        "odd_power_before": before["odd"],
        "odd_power_after": after["odd"],
        "even_nonzero_power_before": before["even_nonzero"],
        "even_nonzero_power_after": after["even_nonzero"],
        "zero_power_after": after["zero"],
        "total_power_relative_drift": float(
            abs(after["total"] - before["total"]) / before["total"]
        ),
        "closure_max_after": float(closure.max()),
    }


def save_figures(records: dict, arrays: dict, fold_info: dict) -> list[Path]:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Hiragino Sans", "Yu Gothic", "DejaVu Sans"],
        "axes.unicode_minus": False,
    })
    outputs = []

    # 図1：論理構造（実空間波形→FFT→倍音族→縮約候補）
    fig, ax = plt.subplots(figsize=(13.2, 4.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    boxes = [
        (0.02, "局在波\n$Z_H(u)=K^{-1}\\sum e^{i(2j+1)u}$"),
        (0.27, "FFT\n奇数倍音 $1,3,\\ldots,H$\n偶数対照 $2,4,\\ldots,H-1$"),
        (0.53, "先行標準\n$H=63$, 32奇数倍音\n$N_\\chi=512, N_\\eta=16$"),
        (0.78, "論文A v1\n$16\\times8$ 巡回レジスタ\n縮約 $Q$ が未提示"),
    ]
    for x, label in boxes:
        rect = plt.Rectangle((x, 0.31), 0.20, 0.38, facecolor="#eef4ff",
                             edgecolor="#315c9b", linewidth=1.8)
        ax.add_patch(rect)
        ax.text(x + 0.10, 0.50, label, ha="center", va="center", fontsize=12)
    for x in (0.225, 0.475, 0.725):
        ax.annotate("", xy=(x + 0.035, 0.50), xytext=(x, 0.50),
                    arrowprops=dict(arrowstyle="->", lw=2, color="#333333"))
    ax.text(0.88, 0.20, "$k\\mapsto k\\;\\mathrm{mod}\\;16$ は候補であり、\n非線形力学との可換性を要検定",
            ha="center", va="top", fontsize=11, color="#9b2f2f")
    ax.set_title("図A-v2-1　欠けていた導出の位置", fontsize=15, pad=12)
    fig.tight_layout()
    p = HERE / "figA_v2_01_missing_derivation_v1.png"
    fig.savefig(p, dpi=180, bbox_inches="tight"); plt.close(fig); outputs.append(p)

    # 図2：倍音数と局在波形
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 8.0), gridspec_kw={"height_ratios": [2, 1]})
    x = np.linspace(-np.pi, np.pi, 2048, endpoint=False)
    for H in (1, 7, 15, 31, 63):
        odds = np.arange(1, H + 1, 2)
        z = np.exp(1j * np.outer(odds, x)).sum(axis=0) / len(odds)
        axes[0].plot(x / np.pi, np.abs(z) ** 2, label=f"H={H}", lw=1.6)
    axes[0].set(xlabel="$u/\\pi$", ylabel="$|Z_H(u)|^2$",
                title="奇数倍音を増やすと局在ピークが細くなる")
    axes[0].legend(ncol=5); axes[0].grid(alpha=0.25)
    odds63 = np.arange(1, 64, 2)
    axes[1].stem(odds63, np.ones_like(odds63), linefmt="#b23a48", markerfmt=" ", basefmt=" ")
    axes[1].set(xlabel="内在倍音 k", ylabel="成分", title="先行標準 H=63：奇数倍音32本")
    axes[1].set_xlim(0, 65); axes[1].grid(alpha=0.25)
    fig.tight_layout()
    p = HERE / "figA_v2_02_harmonic_localization_v1.png"
    fig.savefig(p, dpi=180); plt.close(fig); outputs.append(p)

    # 図3：実エンジンに入れた波形プロファイル
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    for H in H_VALUES:
        key = f"H{H}_full"
        prob = arrays[key]["chi_prob"]
        ax.plot(np.arange(prob.size) / prob.size, prob * prob.size, label=f"H={H}")
    ax.set(xlabel="u / 2π", ylabel="局所パワー / 平均",
           title="同じ総パワーでも倍音上限が局所波形を変える")
    ax.legend(ncol=3); ax.grid(alpha=0.25)
    fig.tight_layout()
    p = HERE / "figA_v2_03_engine_wave_profiles_v1.png"
    fig.savefig(p, dpi=180); plt.close(fig); outputs.append(p)

    # 図4：H掃引
    hs = np.array(H_VALUES)
    peak = np.array([records[f"H{h}_full"]["initial"]["chi_peak_over_mean"] for h in hs])
    prf = np.array([records[f"H{h}_full"]["initial"]["chi_PR_fraction"] for h in hs])
    rate = np.array([records[f"H{h}_full"]["initial"]["rate_relative_norm"] for h in hs])
    outside = np.array([records[f"H{h}_full"]["initial"]["rate_power_outside_abs_k_le_7"] for h in hs])
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2))
    for ax, y, ylabel in (
        (axes[0, 0], peak, "ピーク / 平均"),
        (axes[0, 1], prf, "PR / $N_\\chi$"),
        (axes[1, 0], rate, "$||\\dot C||/||C||$"),
        (axes[1, 1], outside, "$|k|>7$ の速度パワー比"),
    ):
        ax.plot(hs, y, "o-", lw=1.8)
        ax.set_xlabel("最高奇数倍音 H"); ax.set_ylabel(ylabel); ax.grid(alpha=0.3)
        ax.set_xscale("symlog", linthresh=1)
    fig.suptitle("図A-v2-4　倍音上限を変えると局在と頂点応答が変わる", fontsize=14)
    fig.tight_layout()
    p = HERE / "figA_v2_04_harmonic_sweep_metrics_v1.png"
    fig.savefig(p, dpi=180); plt.close(fig); outputs.append(p)

    # 図5：63次倍音の16剰余への折畳み
    evens, odds = harmonic_sets(63)
    fig, ax = plt.subplots(figsize=(11.5, 5.6))
    odd_counts = np.array([sum(k % 16 == r for k in odds) for r in range(16)])
    even_counts = np.array([sum(k % 16 == r for k in evens) for r in range(16)])
    ax.bar(np.arange(16) - 0.18, odd_counts, width=0.36, label="奇数倍音", color="#b23a48")
    ax.bar(np.arange(16) + 0.18, even_counts, width=0.36, label="偶数倍音", color="#315c9b")
    ax.axvspan(-0.48, 0.48, color="#f0b429", alpha=0.25)
    ax.annotate("16, 32, 48次は零モードへ\n（現読出しではボゾン帯外）",
                xy=(0, even_counts[0]), xytext=(3.2, 3.2),
                arrowprops=dict(arrowstyle="->", color="#8a5b00"), fontsize=11)
    ax.set(xlabel="$k\\;\\mathrm{mod}\\;16$", ylabel="畳み込まれる倍音数",
           title=f"H=63 の剰余縮約：初期読出し差={fold_info['readout_max_abs_difference']:.3e}")
    ax.set_xticks(range(16)); ax.legend(); ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    p = HERE / "figA_v2_05_mod16_folding_v1.png"
    fig.savefig(p, dpi=180); plt.close(fig); outputs.append(p)

    return outputs


def main() -> None:
    t0 = time.time()
    emit("論文A v2 倍音構造訂正再実験 v1")
    emit(f"N={N}, delta={DELTA}, full={FULL_NK}x{FULL_NETA}, small={SMALL_NK}x{SMALL_NETA}")
    Z0c, wp0, seed_state = parent_data(N)

    # v1 の構成を独立に再構成し、正本ビルダーと bitwise 比較する。
    single_C, _, _ = make_packet(Z0c, seed_state, 1, SMALL_NK, SMALL_NETA)
    standard_eng, _, _ = F.build_standard_universe(N, DELTA, Nn=SMALL_NK, Neta=SMALL_NETA)
    control_bitwise = bool(np.array_equal(single_C, standard_eng.C2()))
    emit(f"v1 single16 control bitwise: {control_bitwise}")
    if not control_bitwise:
        raise RuntimeError("正本 single16 対照が不一致なので中止")

    records: dict[str, dict] = {}
    arrays: dict[str, dict[str, np.ndarray]] = {}
    initial_states: dict[str, np.ndarray] = {"single16": single_C}

    m, a = state_metrics(single_C, wp0)
    records["single16"] = {"initial": m}
    arrays["single16"] = a
    emit(f"single16: R={m['readout_mean']:.12g} peak={m['chi_peak_over_mean']:.6g} "
         f"rate={m['rate_relative_norm']:.6g} nsub={m['estimated_RK4_substeps']}")

    full_mode_lists: dict[str, tuple[int, ...]] = {}
    for H in H_VALUES:
        key = f"H{H}_full"
        C, evens, odds = make_packet(Z0c, seed_state, H, FULL_NK, FULL_NETA)
        initial_states[key] = C
        full_mode_lists[key] = tuple(sorted(set(evens + odds)))
        m, a = state_metrics(C, wp0)
        records[key] = {"initial": m, "even_modes": list(evens), "odd_modes": list(odds)}
        arrays[key] = a
        emit(f"{key}: R={m['readout_mean']:.12g} peak={m['chi_peak_over_mean']:.6g} "
             f"PR/N={m['chi_PR_fraction']:.6g} rate={m['rate_relative_norm']:.6g} "
             f"out7={m['rate_power_outside_abs_k_le_7']:.6g} "
             f"nsub={m['estimated_RK4_substeps']}")

    selected63 = full_mode_lists["H63_full"]
    fold_C = fold_selected_k(initial_states["H63_full"], selected63)
    initial_states["H63_fold16"] = fold_C
    m, a = state_metrics(fold_C, wp0)
    records["H63_fold16"] = {"initial": m, "projection": "Q=sum/sqrt(n_r), selected k=1..63"}
    arrays["H63_fold16"] = a
    emit(f"H63_fold16: R={m['readout_mean']:.12g} peak={m['chi_peak_over_mean']:.6g} "
         f"rate={m['rate_relative_norm']:.6g} nsub={m['estimated_RK4_substeps']}")

    # 読出し不変予測（full系列）とfold零モード効果。
    full_R = np.vstack([arrays[f"H{H}_full"]["R"] for H in H_VALUES])
    full_R_spread = float(np.max(np.ptp(full_R, axis=0)))
    fold_R_diff = float(np.max(np.abs(arrays["H63_full"]["R"] - arrays["H63_fold16"]["R"])))

    # 事前登録した1更新。全系列を実エンジンで走らせる。
    emit("one-step evolution begins")
    for key in ("single16",) + tuple(f"H{H}_full" for H in H_VALUES) + ("H63_fold16",):
        s0 = time.time()
        C1, step_metrics = one_step(initial_states[key], wp0)
        records[key]["one_step"] = step_metrics
        arrays[key]["C1"] = C1
        emit(f"  {key}: {time.time()-s0:.2f}s drift={step_metrics['total_power_relative_drift']:.3e} "
             f"closure={step_metrics['closure_max_after']:.3e}")

    # Q∘step_full と step_16∘Q の可換性。
    Q_after_full = fold_selected_k(arrays["H63_full"]["C1"], selected63)
    step_after_Q = arrays["H63_fold16"]["C1"]
    comm_abs = float(np.linalg.norm(Q_after_full - step_after_Q))
    comm_rel = comm_abs / max(float(np.linalg.norm(Q_after_full)), np.finfo(float).tiny)
    fold_info = {
        "definition": "Q(C)_r=sum_{k in 1..63, k mod 16=r} C_k/sqrt(n_r); eta0 copied",
        "full_readout_spread_max_abs": full_R_spread,
        "readout_max_abs_difference": fold_R_diff,
        "commutator_abs_L2": comm_abs,
        "commutator_relative_L2": comm_rel,
        "commutes_at_1e-10": bool(comm_rel <= 1.0e-10),
    }
    emit(f"full R spread={full_R_spread:.3e}; fold R difference={fold_R_diff:.3e}")
    emit(f"Q(step(full)) vs step(Q(full)): relative L2={comm_rel:.6g}")

    source_paths = {
        "experiment": Path(__file__).resolve(),
        "preregistration": PREREG,
        "unified_F1": Path(F.__file__).resolve(),
        "stage2_vertex": Path(F.s2.__file__).resolve(),
        "parent_generator": Path(F.gen3.__file__).resolve(),
        "lowrank_parent": Path(F.abl.__file__).resolve(),
    }
    source_hashes = {name: {"path": str(path), "sha256": sha256(path)}
                     for name, path in source_paths.items()}

    report = {
        "title": "論文A v2 倍音構造訂正再実験 v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed",
        "conditions": {
            "N": N, "M": N * (N - 1) // 2, "delta": DELTA, "parent_seed": SEED,
            "full_Nk": FULL_NK, "full_Neta": FULL_NETA,
            "small_Nk": SMALL_NK, "small_Neta": SMALL_NETA,
            "eta_index": ETA_INDEX,
            "normalization": "fixed total pump-family and seed-family powers",
            "complex_wave": "positive-k analytic harmonic packet; carrier stripped",
        },
        "control": {"single16_vs_build_standard_universe_bitwise": control_bitwise},
        "records": records,
        "fold_test": fold_info,
        "predictions": {
            "P1_full_readout_invariant": bool(full_R_spread <= 1.0e-13),
            "P2_peak_monotone_non_decreasing": bool(all(
                records[f"H{b}_full"]["initial"]["chi_peak_over_mean"] >=
                records[f"H{a}_full"]["initial"]["chi_peak_over_mean"]
                for a, b in zip(H_VALUES[:-1], H_VALUES[1:])
            )),
            "P2_PR_fraction_monotone_non_increasing": bool(all(
                records[f"H{b}_full"]["initial"]["chi_PR_fraction"] <=
                records[f"H{a}_full"]["initial"]["chi_PR_fraction"]
                for a, b in zip(H_VALUES[:-1], H_VALUES[1:])
            )),
            "P3_single_vs_H63_rate_relative_difference": float(
                abs(records["H63_full"]["initial"]["rate_relative_norm"] -
                    records["H1_full"]["initial"]["rate_relative_norm"]) /
                max(records["H1_full"]["initial"]["rate_relative_norm"], np.finfo(float).tiny)
            ),
            "P4_H63_rate_has_outside_abs_k7": bool(
                records["H63_full"]["initial"]["rate_power_outside_abs_k_le_7"] > RATE_THRESHOLD
            ),
            "P5_fold_dynamics_noncommuting": bool(comm_rel > 1.0e-10),
        },
        "scope_limit": (
            "eta=0 only; the 16-to-8 eta reduction is not validated here. "
            "This is a structural response and one-step correction experiment, not a replacement "
            "for all v1 long-time sweeps."
        ),
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
        },
        "source_hashes": source_hashes,
        "elapsed_seconds_before_figures": time.time() - t0,
    }

    figures = save_figures(records, arrays, fold_info)
    report["figures"] = [p.name for p in figures]
    report["elapsed_seconds_total"] = time.time() - t0

    json_path = HERE / "harmonic_structure_correction_results_v1.json"
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2, default=jsonable) + "\n",
                         encoding="utf-8")

    npz_data: dict[str, np.ndarray] = {
        "H_values": np.array(H_VALUES),
        "full_R": full_R,
        "Q_after_full_C1": Q_after_full,
        "step_after_Q_C1": step_after_Q,
    }
    for key, vals in arrays.items():
        npz_data[f"{key}__R"] = vals["R"]
        npz_data[f"{key}__chi_prob"] = vals["chi_prob"]
        npz_data[f"{key}__rate_spectrum"] = vals["rate_spectrum"]
        npz_data[f"{key}__freqs"] = vals["freqs"]
    npz_path = HERE / "harmonic_structure_correction_arrays_v1.npz"
    np.savez_compressed(npz_path, **npz_data)

    log_path = HERE / "run_harmonic_structure_correction_v1.log"
    emit(f"saved {json_path.name}, {npz_path.name}, {len(figures)} figures")
    emit(f"elapsed={time.time()-t0:.2f}s")
    log_path.write_text("\n".join(LOG_LINES) + "\n", encoding="utf-8")

    # JSON/NPZ/図/ログを含む実行後マニフェスト。
    artifact_paths = [json_path, npz_path, log_path, *figures]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "artifacts": [
            {"file": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in artifact_paths
        ],
        "sources": source_hashes,
    }
    manifest_path = HERE / "manifest_harmonic_structure_correction_v1.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                             encoding="utf-8")

if __name__ == "__main__":
    main()
