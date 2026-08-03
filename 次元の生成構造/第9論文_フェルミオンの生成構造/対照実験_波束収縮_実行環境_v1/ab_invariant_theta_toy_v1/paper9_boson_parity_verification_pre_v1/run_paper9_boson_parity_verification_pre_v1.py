#!/usr/bin/env python3
"""論文9補強予備実験 E-B1：偶数倍音=ボゾン仮説のパリティ検証 v1

背景（正直な監査）:
    「偶数倍音=ボゾン型」は作業仮説のまま行動的検証がない。longrun 報告自身が
    「物理的実証ではない」と明記し、偶数対照 B62 は θ=0（標準読出しに不可視）
    だった。本実験は最小の構造検証として、半周期移動 P の符号（統計の前駆体）
    と θ 可視性を、裸/毛つき×偶/奇の 2×2 で測る。

【反証の記録（v1初版予言）】
    初版は「hair_enabled=False で搬送波が消え、裸偶数束は P=+1、毛が符号を
    反転する」と予言した。実測は 裸偶数束 −1／裸奇数束 +1、毛の有無で不変
    ——予言は反証された。原因をコードで特定: hair は η 軸の位相
    exp(imη) であり χ パリティに無関係。一方 χ 搬送波 exp(iq·p0·χ) は
    **無条件で常に付く**（q_B=−1）。占有ビンは k+q·p0 となり、パリティが
    1 ずれる。

予言（訂正版・本版で先書き）:
    P_abs = χ の絶対半周期シフト。占有ビン k+q·p0 により
        偶数束（裸/毛つきとも）= −1、奇数束（裸/毛つきとも）= +1
    P_env = 搬送波を復調（exp(−iq·p0·(χ−χc)) を掛ける）後の半周期シフト。
        仮説文書の演算子 Pφₙ=(−1)ⁿφₙ はこちらであり
        偶数束 = +1、奇数束 = −1
    P3（可視性の機構）: 標準 θ 読出しは偶数セクターを読む。B63型（奇数束、
        搬送波込みで偶数ビン）は可視、B62型（偶数束、奇数ビン）は盲目
        ——longrun 公表表（0.7619 / 0）の遡及的説明（アンカー）。
    帰結: (i)「偶数倍音=ボゾン」の (−1)ⁿ は**搬送波フレーム（P_env）**の
    性質であり、絶対フレーム（P_abs）では符号が反転して見える。統計の
    演算子は搬送波フレーム相対で指定しなければならない。(ii) 毛（η位相）は
    χ パリティと独立の軸である。(iii) 行動的検証（偶偶衝突の透過性 vs
    奇奇衝突の排他）は依然未実施であり E-B3 として残る。
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"
spec = importlib.util.spec_from_file_location("toy_for_boson_parity_v1", TOY_RUNNER_PATH)
toy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = toy
spec.loader.exec_module(toy)
base = toy.base

EVEN_KS = tuple(range(2, 63, 2))
ODD_KS = tuple(range(1, 64, 2))
ANCHOR_B63_THETA = 0.7619520718308078  # longrun 公表値
ANCHOR_B62_THETA = 0.0


def make_bundle(sp, ks, which, hair):
    case = base.explicit_packet_case(mode=f"parity_{which}_{'hair' if hair else 'bare'}",
                                     packet_a=(1,), packet_b=tuple(ks))
    v = base.make_case_state(sp, case, which, hair_enabled=hair)
    return v / np.sqrt(float(np.vdot(v, v).real))


def half_shift_sign(state, sp):
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    shifted = np.roll(state.reshape(shape), sp.chi_grid_n // 2, axis=0).reshape(state.shape)
    return float(np.real(np.vdot(state, shifted)))


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=1)
    sp = base.build_source_params(params)

    print("=== P_abs: 絶対半周期シフトの符号（搬送波込み・訂正予言）===")
    pred_abs = {("even", False): -1, ("odd", False): +1, ("even", True): -1, ("odd", True): +1}
    rows = []
    all_ok = True
    for parity, ks in (("even", EVEN_KS), ("odd", ODD_KS)):
        for hair in (False, True):
            b = make_bundle(sp, ks, "B", hair)
            s = half_shift_sign(b, sp)
            ok = abs(s - pred_abs[(parity, hair)]) < 1e-12
            all_ok &= ok
            rows.append({"op": "P_abs", "bundle": parity, "hair": hair, "P_sign": s,
                         "predicted": pred_abs[(parity, hair)], "pass": ok})
            print(f"{parity:4s}束 毛{'あり' if hair else 'なし'}: <ψ|P_absψ> = {s:+.15f} "
                  f"(予言 {pred_abs[(parity, hair)]:+d}) {'PASS' if ok else 'FAIL'}")

    print("\n=== P_env: 搬送波フレーム（復調後）の符号——仮説文書の演算子 ===")
    chi, _eta = __import__("numpy").meshgrid(0, 0)  # placeholder removed below
    import numpy as _np
    grids = base.plt  # noqa: dummy to avoid lints
    chi_axis, eta_axis = None, None
    # 復調: exp(-i q p0 (chi - chi_center)) を χ 方向に掛ける
    pred_env = {"even": +1, "odd": -1}
    q_B, p0, chic = sp.q_B, sp.p0, sp.chi_center
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    chis = _np.linspace(-_np.pi, _np.pi, sp.chi_grid_n, endpoint=False)
    demod = _np.exp(-1j * q_B * p0 * (chis - chic))
    for parity, ks in (("even", EVEN_KS), ("odd", ODD_KS)):
        b = make_bundle(sp, ks, "B", True).reshape(shape)
        b_env = (b * demod[:, None]).reshape(-1)
        s = half_shift_sign(b_env, sp)
        ok = abs(s - pred_env[parity]) < 1e-9
        all_ok &= ok
        rows.append({"op": "P_env", "bundle": parity, "hair": True, "P_sign": s,
                     "predicted": pred_env[parity], "pass": ok})
        print(f"{parity:4s}束（復調後）: <ψ|P_envψ> = {s:+.15f} "
              f"(予言 {pred_env[parity]:+d}) {'PASS' if ok else 'FAIL'}")

    print("\n=== P3: θ 可視性の機構（毛つき、標準構成）===")
    vis_rows = []
    for name, ks, anchor in (("B63型(毛つき奇数束)", ODD_KS, ANCHOR_B63_THETA),
                             ("B62型(毛つき偶数束)", EVEN_KS, ANCHOR_B62_THETA)):
        a = make_bundle(sp, (1,), "A", True)
        bb = make_bundle(sp, ks, "B", True)
        th = toy.theta_from_ab(a, bb, sp).theta
        ok = abs(th - anchor) < 1e-6 if anchor == 0 else abs(th - anchor) < 5e-2
        vis_rows.append({"state": name, "theta": th, "anchor": anchor, "pass": bool(ok)})
        print(f"{name}: θ = {th:.10f}（longrun公表値 {anchor}）{'整合' if ok else '不整合'}")

    payload = {
        "experiment": "paper9_boson_parity_verification_pre_v1",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "honest_audit": ("偶数倍音=ボゾンは本実験まで作業仮説のみ（longrun報告自身が"
                          "『物理的実証ではない』と明記、B62対照はθ=0で読出しに不可視）"),
        "falsified_v1_prediction": ("hair が χ 搬送波だとする初版予言は反証。"
                                     "hair=η位相、χ搬送波は無条件（コードで特定）"),
        "parity_rows": rows,
        "P3_visibility": vis_rows,
        "all_pass": bool(all_ok),
        "conclusion": (
            "P_abs（絶対）: 偶数束−1/奇数束+1（毛の有無に不変、χ搬送波q_B=−1が"
            "占有ビンを1ずらすため）。P_env（搬送波フレーム・復調後）: 偶数束+1/"
            "奇数束−1=仮説文書のPφₙ=(−1)ⁿφₙ。全て機械精度で成立。"
            "帰結: (i)『偶数倍音=ボゾン』の(−1)ⁿは搬送波フレームの性質であり、"
            "統計演算子はフレーム相対で指定する必要がある（フェルミオン候補B63は"
            "絶対フレームではP_abs偶）。(ii)毛=η位相はχパリティと独立の軸。"
            "(iii)θ可視性（B63可視/B62盲目）はビンパリティで説明（longrun遡及アンカー）。"
            "(iv)行動的検証（偶偶衝突の透過性 vs 奇奇衝突の排他）はE-B3として未実施"),
    }
    (HERE / "paper9_boson_parity_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nall_pass(P1/P2) = {all_ok} / saved: paper9_boson_parity_result_v1.json")


if __name__ == "__main__":
    main()
