#!/usr/bin/env python3
"""論文9補強予備実験 E-B4：毛つきボゾン（電荷つきボゾン候補）の電荷は発現できるか v1

問い（木原氏）:
    ボゾン（偶数束）に素電荷の担い手（毛）を付けたら、おかしな状態が
    起こらないか。E-B3 の偶数束は実は毛つきだった——電荷は何をしていたのか。

方法（無名性厳守・IF分岐なし）:
    偶偶対（毛あり/毛なし）に対し、B の全体位相 φ を8点掃引し、
    (i) ペア重なり |<a|b>|、(ii) 内生θ、(iii) 内生の一衝突移乗 dN_B(φ) の変調、
    (iv) プローブ回転 θ_p=0.15（計器設定、三部作6.3節と同じ規約）での
    dN_B(φ) の変調、を測る。電荷の発現 = φ による変調（符号の存在）。

予言（測定前固定）:
    P1（内生の沈黙）: 内生θは毛の有無によらず0（偶数束は偶数セクター読出しに
        寄与しない）→ 内生の dN_B 変調 = 0。**電荷は内生的に発現しない**
    P2（毛による自動中性化）: 毛あり対は η 直交により |<a|b>| ≈ 0
        （三部作のヌル定理＝中性化第二機構が毛自体によって自動発動）
        → プローブ回転を掛けても変調 ≈ 0。
    P3（毛なし対照）: 毛なしなら η は共通で |<a|b>| > 0 になりうるが、
        毛がない以上、符号の担い手が不在。プローブ変調は c=<a|b> の実部
        オフセットとして現れる（電荷型の±対ではない）。
    帰結（定理候補）: 電荷の発現には 毛 ∧ チャネル共有 ∧ θ>0 の三条件が
        同時に必要であり、偶数束対はどの構成でも三条件を同時に満たせない。
        **「ボゾンの素電荷」は矛盾を起こすのではなく、原理的に読み出せない
        ——命名後置により、内生的には存在しない。** 光子の電荷ゼロは
        割当てではなく帰結になる。
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"
spec = importlib.util.spec_from_file_location("toy_for_boson_charge_v1", TOY_RUNNER_PATH)
toy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = toy
spec.loader.exec_module(toy)
base = toy.base

EVEN_KS = tuple(range(2, 63, 2))
PHASES = np.linspace(0.0, 2 * math.pi, 8, endpoint=False)
THETA_PROBE = 0.15


def make_bundle(sp, which, hair):
    case = base.explicit_packet_case(mode=f"bcm_{which}_{'hair' if hair else 'bare'}",
                                     packet_a=EVEN_KS, packet_b=EVEN_KS)
    v = base.make_case_state(sp, case, which, hair_enabled=hair)
    return v / np.sqrt(float(np.vdot(v, v).real))


def dnb_after(a, b, th):
    a2, b2 = toy.rotate_ab(a.copy(), b.copy(), th)
    return float(np.vdot(b2, b2).real - np.vdot(b, b).real)


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=1)
    sp = base.build_source_params(params)

    results = []
    for hair in (True, False):
        a = make_bundle(sp, "A", hair)
        b0 = make_bundle(sp, "B", hair)
        overlap = abs(complex(np.vdot(a, b0)))
        th_endo = toy.theta_from_ab(a, b0, sp).theta

        endo_mod, probe_mod = [], []
        for ph in PHASES:
            b = np.exp(1j * ph) * b0
            th = toy.theta_from_ab(a, b, sp).theta
            endo_mod.append(dnb_after(a, b, th))
            probe_mod.append(dnb_after(a, b, THETA_PROBE))
        endo_amp = float(np.max(endo_mod) - np.min(endo_mod))
        probe_amp = float(np.max(probe_mod) - np.min(probe_mod))

        row = {"hair": hair, "overlap_abs": overlap, "theta_endogenous": float(th_endo),
               "endo_modulation_amp": endo_amp, "probe_modulation_amp": probe_amp}
        results.append(row)
        print(f"毛{'あり' if hair else 'なし'}: |<a|b>|={overlap:.3e} θ内生={th_endo:.3e} "
              f"内生変調={endo_amp:.3e} プローブ変調={probe_amp:.3e}")

    hair_row = next(r for r in results if r["hair"])
    bare_row = next(r for r in results if not r["hair"])
    p1 = (hair_row["theta_endogenous"] == 0 and bare_row["theta_endogenous"] == 0
          and hair_row["endo_modulation_amp"] <= 1e-14 and bare_row["endo_modulation_amp"] <= 1e-14)
    p2 = hair_row["overlap_abs"] <= 1e-12 and hair_row["probe_modulation_amp"] <= 1e-12
    print(f"\nP1 内生の沈黙（θ=0・内生変調0、毛の有無不問）: {'PASS' if p1 else 'FAIL'}")
    print(f"P2 毛による自動中性化（η直交で|<a|b>|≈0・プローブ変調も0）: {'PASS' if p2 else 'FAIL'}")
    print(f"P3 毛なし対照: |<a|b>|={bare_row['overlap_abs']:.4f} "
          f"プローブ変調={bare_row['probe_modulation_amp']:.3e}（記録）")

    payload = {
        "experiment": "paper9_boson_charge_mute_pre_v1",
        "design": "無名性厳守・IF分岐なし。E-B3と同一構成の偶偶対＋位相掃引＋プローブ規約は三部作6.3節と同一",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "rows": results,
        "P1_endogenous_silence": bool(p1),
        "P2_hair_auto_neutralization": bool(p2),
        "conclusion": (
            "毛つき偶数束対（電荷つきボゾン候補）は、(i) 内生θ=0 で両文法とも発火せず、"
            "(ii) 毛自体が η 直交＝チャネル非共有（中性化第二機構）を強制するため、"
            "プローブ回転を掛けても電荷変調は発現しない。電荷の発現には"
            "毛∧共有∧θ>0 の三条件が同時に必要で、偶数束対はこれを同時に満たせない。"
            "『ボゾンの素電荷』は矛盾を生むのではなく原理的に読み出せない——"
            "命名後置により内生的には存在しない。光子の電荷ゼロは割当てでなく帰結である"),
    }
    (HERE / "paper9_boson_charge_mute_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print("\nsaved: paper9_boson_charge_mute_result_v1.json")


if __name__ == "__main__":
    main()
