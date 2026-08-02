#!/usr/bin/env python3
"""交差項電荷文法 v2/v3 の正式ランナー化

README の「v2（無効）」「v3（成立）」を生成したインライン実行の再現コード。

v2: 搬送波なし＋パリティ反転マスク構成の無効性の診断
    （搬送波なしでは倍音が想定 chi ビンに乗らず、分率が全てゼロになる）
v3: 搬送波なし・プローブ回転での電荷文法検証
    （kappa = |<A5|B7>| = 0.5 厳密、中性性、符号反転。
     既コミットの cross_term_charge_hairoff_M4_fits_v1.json の生成元。
     上位互換の完全版は paperA_two_grammar_evidence_v1 を参照）
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
V1_PATH = HERE / "run_cross_term_charge_pre_v1.py"

spec = importlib.util.spec_from_file_location("ct_for_v2v3_v1", V1_PATH)
ct = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ct
ct.__name__ = "ct_for_v2v3_v1"
spec.loader.exec_module(ct)

toy, base = ct.toy, ct.base


def mk_hairoff(k: int, which: str, sp) -> np.ndarray:
    case = base.explicit_packet_case(mode=f"v2v3_{which}_{k}", packet_a=(k,), packet_b=(k,))
    return ct.unit_norm(base.make_case_state(sp, case, which, hair_enabled=False))


def masked_fraction_m4(state: np.ndarray, sp) -> float:
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    sf = np.fft.fft(state.reshape(shape), axis=0, norm="ortho")
    freqs = np.rint(np.fft.fftfreq(sp.chi_grid_n, d=1.0 / sp.chi_grid_n)).astype(int)
    f = np.abs(freqs)
    mask = (f >= 3) & (f % 2 == 1)
    p = np.sum(np.abs(sf) ** 2, axis=1)
    return float(np.sum(p[mask])) / float(np.sum(p))


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=1)
    sp = base.build_source_params(params)

    a1, a5 = mk_hairoff(1, "A", sp), mk_hairoff(5, "A", sp)
    b1, b7 = mk_hairoff(1, "B", sp), mk_hairoff(7, "B", sp)

    # ---- v2 無効性の診断 ----
    v2 = {
        "A1_fraction_M4": masked_fraction_m4(a1, sp),
        "A5_fraction_M4": masked_fraction_m4(a5, sp),
        "B7_fraction_M4": masked_fraction_m4(b7, sp),
    }
    v2_invalid = v2["A5_fraction_M4"] < 0.5 and v2["B7_fraction_M4"] < 0.5
    print("v2 diagnosis (hair-off fractions under parity-flip mask):", v2)
    print("v2 invalid (fractions collapse, as recorded):", v2_invalid)

    # ---- v3 電荷文法（インライン版と同一プロトコル）----
    c6 = complex(np.vdot(a5, b7))
    F = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
    PH = [2 * math.pi * k / 8 for k in range(8)]
    fits = []
    for th_p in (0.1, 0.2):
        s2t = math.sin(2 * th_p)
        for fa in F:
            a0 = ct.unit_norm(math.sqrt(1 - fa) * a1 + math.sqrt(fa) * a5)
            for fb in F:
                flows = []
                for phi in PH:
                    b0 = ct.unit_norm(
                        math.sqrt(1 - fb) * b1 + np.exp(1j * phi) * math.sqrt(fb) * b7
                    )
                    a2, b2 = toy.rotate_ab(a0.copy(), b0.copy(), th_p)
                    flows.append(float(np.vdot(b2, b2).real) - 1.0)
                off, m_mod, delta = ct.fit_cosine(np.asarray(PH), np.asarray(flows))
                amp = math.sqrt(fa * fb)
                kappa = m_mod / (s2t * amp) if s2t * amp > 1e-15 else None
                fits.append(dict(theta_p=th_p, fa=fa, fb=fb, offset=off, M=m_mod,
                                 delta=delta, kappa=kappa))
    charged = [f for f in fits if f["kappa"] is not None]
    ks = np.asarray([f["kappa"] for f in charged])
    neutral_max = max(f["M"] for f in fits if f["kappa"] is None)
    print(
        f"v3: kappa mean={ks.mean():.15f} spread={ks.max()-ks.min():.3e}"
        f" (prediction |c6|={abs(c6):.15f})"
    )
    print(f"v3: neutral max modulation = {neutral_max:.3e}")

    payload = {
        "experiment": "cross_term_charge_v2_v3_formalized_v1",
        "v2_diagnosis": {**v2, "invalid_confirmed": v2_invalid},
        "v3": {
            "kappa_mean": float(ks.mean()),
            "kappa_spread": float(ks.max() - ks.min()),
            "abs_c6": abs(c6),
            "neutral_max_modulation": neutral_max,
            "charged_cells": len(charged),
        },
        "note": "完全版（4プローブ角・三項分解・異ノルム）は paperA_two_grammar_evidence_v1 が正",
    }
    (HERE / "cross_term_charge_v2_v3_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
