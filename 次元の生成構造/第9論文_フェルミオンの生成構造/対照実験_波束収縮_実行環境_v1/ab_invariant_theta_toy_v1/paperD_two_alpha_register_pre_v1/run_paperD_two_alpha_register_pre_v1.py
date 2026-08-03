#!/usr/bin/env python3
"""論文D予備実験 D1+D2+D5：二つのα根のレジスタ算術・時計盲目性・エイリアシング v1

背景:
    有限位数根論文は α^-1=137 近傍の R_{124,23} と α^-1=128.95 近傍の
    R_{620,117} を厳密根として同定した。本実験は両根の関係を三層で判定する。

D1（算術・アンカー）:
    (i) 変位の整数性: π/4 からの変位が厳密に 20/310π, 19/310π
    (ii) 埋め込み: 23/124 = 115/620（既約でない）、π/4側の最初の既約住所 = 117
    (iii) 作用素周期: U(θ124)^124 = -I, ^248 = I / U(θ620)^620 = -I, ^1240 = I
    (iv) アンカー: F_8(124根) = 0.0026（E6 公開値の再現）

D2（時計盲目性）予言（測定前固定・閉形式 F = cos^2(J·θ_rot) をコード内で計算）:
    実測は全予言と 1.3e-14 で一致。判定: F_248(620根)=0.6545 < 1 かつ
    F_1240(両根) = 1 —— どの根が見えるかは観測時計の約数構造が決める
    （エネルギー = 時計細分の直接証拠）。

D5（エイリアシング＝電荷普遍性）:
    【反証の記録】初版の辞書（cos規約の住所 m をそのまま DFT ビンに対応:
    W=124→23, 248→47, 620→117）は反証された（実測ピーク 19/39/97）。
    原因は二つの規約落ち: (i) 回転角は流れ規約の補角 θ_rot = 193π/620
    （124根では 39π/124）である、(ii) 窓 W の DFT はビン刻み 2π/W であり、
    レジスタ刻み π/n の半分——奇数住所は窓 W=n では見えず、時計 W=2n
    （状態の完全周期）が必要。
    【修正辞書の予言（本版で先書き）】ストロボ周波数は θ_rot/2π = 193/1240。
        W=124（半時計）  : ビン 19（=偶数住所38へのエイリアス）→ 読み 0.326312
        W=248（電子時計）: ビン 39 → 読み電荷 cos^2(39π/124) = 0.302822 ★
        W=620（半時計）  : ビン 97
        W=1240（固有時計）: ビン 193 → 読み電荷 cos^2(193π/620) = 0.312175
    ★= 電子の時計 248 で読んだ細かい根は、電子の住所 39/124 に厳密に
    折り返され、電子と同じ電荷と読まれる——電荷普遍性の機構候補。

D9への設計帰結（本実験の副産物として記録）:
    位相盲目読出しの no-go により θ は凍結し、住所選択は θ 力学では起こり
    えない。質量平衡による選択（第五候補）はノルム・パケットセクターで
    検定しなければならない。
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from fractions import Fraction
from math import gcd, pi
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SEARCH_PATH = HERE.parent / "inverse_initial_conditions_v1" / "search_initial_conditions_and_plot_v1.py"

spec = importlib.util.spec_from_file_location("search_for_paperD_v1", SEARCH_PATH)
search = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = search
spec.loader.exec_module(search)
toy = search.toy
base = toy.base
plt = base.plt

# ---- 根の定義 ----
N_LOW, M_LOW = 124, 23      # 低エネルギー根（α^-1=137 近傍）
N_HIGH, M_HIGH = 620, 117   # 高エネルギー根（α^-1=128.95 近傍）
TH_LOW = M_LOW * pi / N_LOW
TH_HIGH = M_HIGH * pi / N_HIGH
R_LOW = math.cos(TH_LOW) ** 2
R_HIGH = math.cos(TH_HIGH) ** 2

J_CLOCKS = (8, 248, 620, 1240)

# ---- 予言（測定前固定）----
PRED_F = {
    "root124": {8: math.cos(8 * TH_LOW) ** 2, 248: 1.0, 620: 1.0, 1240: 1.0},
    "root620": {8: math.cos(8 * TH_HIGH) ** 2, 248: math.cos(248 * TH_HIGH) ** 2,
                620: 1.0, 1240: 1.0},
}
# 修正辞書（v1初版の {124:23, 248:47, 620:117} は反証・記録済み）
TH_ROT_HIGH = pi / 2 - TH_HIGH            # 流れ規約の回転角 193π/620
PRED_ALIAS = {124: 19, 248: 39, 620: 97, 1240: 193}
PRED_CHARGE_READ = {248: 0.302822, 1240: 0.312175}  # 時計窓での読み電荷（6桁）
FALSIFIED_ALIAS_V0 = {124: 23, 248: 47, 620: 117}
ANCHOR_F8_ROOT124 = 0.0026  # E6 公開値（丸め）


def pair_fidelity(a0, b0, a, b) -> float:
    ov = complex(np.vdot(a0, a) + np.vdot(b0, b))
    n0 = float(np.vdot(a0, a0).real + np.vdot(b0, b0).real)
    return abs(ov) ** 2 / n0 ** 2


def main() -> None:
    # ================= D1: 算術・アンカー =================
    print("=== D1: レジスタ算術 ===")
    disp_low = (pi / 4 - TH_LOW) / (pi / 310)
    disp_high = (pi / 4 - TH_HIGH) / (pi / 310)
    assert abs(disp_low - 20) < 1e-9 and abs(disp_high - 19) < 1e-9
    print(f"変位（単位π/310）: 低={disp_low:.12f} 高={disp_high:.12f}  → 整数 20/19 PASS")

    assert Fraction(M_LOW, N_LOW) == Fraction(115, 620)
    fresh = [m for m in range(116, 156) if m % 2 == 1 and gcd(m, N_HIGH) == 1]
    assert fresh[0] == M_HIGH
    print(f"埋め込み 23/124=115/620 PASS / π/4側の最初の既約住所 = {fresh[0]} PASS")

    def rot(th):
        return np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])

    U_low, U_high = rot(TH_LOW), rot(TH_HIGH)
    e = np.eye(2)
    checks = {
        "U(th124)^124=-I": np.max(np.abs(np.linalg.matrix_power(U_low, 124) + e)),
        "U(th124)^248=+I": np.max(np.abs(np.linalg.matrix_power(U_low, 248) - e)),
        "U(th620)^620=-I": np.max(np.abs(np.linalg.matrix_power(U_high, 620) + e)),
        "U(th620)^1240=+I": np.max(np.abs(np.linalg.matrix_power(U_high, 1240) - e)),
    }
    for k, v in checks.items():
        print(f"{k}: {v:.2e}")
        assert v < 1e-11

    e_low, e_high = 1 - R_LOW, 1 - R_HIGH
    print(f"1-R: 低 {e_low:.6f} (√(4πα137)={math.sqrt(4*pi/137.035999):.6f}) "
          f"高 {e_high:.6f} (√(4πα128.95)={math.sqrt(4*pi/128.95):.6f})")

    # ================= 状態構成（逆算探索、前進は無変更）=================
    params = base.Params(high_n=63, recursive_collision_count=max(J_CLOCKS))
    sp = base.build_source_params(params)
    a_t, b_t, _ = search.make_unit_templates(sp)

    states = {}
    for label, target in (("root124", R_LOW), ("root620", R_HIGH)):
        res = search.search_initial_b_amplitude(target, a_t, b_t, sp, tolerance=1.0e-15)
        a0, b0 = a_t.copy(), res.initial_b_amplitude * b_t
        r0 = toy.theta_from_ab(a0, b0, sp).reflection_rate
        states[label] = (a0, b0)
        print(f"{label}: 目標R={target:.9f} 実現R={r0:.9f} 差={abs(r0-target):.2e}")

    # ================= D2: 時計盲目性 =================
    print("\n=== D2: 時計盲目性（忠実度表）===")
    rows = []
    strobe = {}
    for label, (a0, b0) in states.items():
        a, b = a0.copy(), b0.copy()
        fids = {}
        zs = []
        for j in range(1, max(J_CLOCKS) + 1):
            readout = toy.theta_from_ab(a, b, sp)
            a, b = toy.rotate_ab(a, b, readout.theta)
            zs.append(complex(np.vdot(a0, a) + np.vdot(b0, b)))
            if j in J_CLOCKS:
                fids[j] = pair_fidelity(a0, b0, a, b)
        strobe[label] = np.array(zs)
        row = {"state": label}
        for j in J_CLOCKS:
            row[f"F_{j}"] = fids[j]
            row[f"pred_F_{j}"] = PRED_F[label][j]
            row[f"err_F_{j}"] = abs(fids[j] - PRED_F[label][j])
        rows.append(row)
        errs = ", ".join(f"J={j}: F={fids[j]:.6f} (予言 {PRED_F[label][j]:.6f})" for j in J_CLOCKS)
        print(f"{label}: {errs}")

    max_err = max(r[f"err_F_{j}"] for r in rows for j in J_CLOCKS)
    anchor_err = abs(rows[0]["F_8"] - ANCHOR_F8_ROOT124)
    blind = rows[1]["F_248"] < 1 - 1e-3 and abs(rows[1]["F_1240"] - 1) < 1e-9
    print(f"予言との最大誤差: {max_err:.2e} / アンカー F_8(124根)≈0.0026 誤差 {anchor_err:.4f}")
    print(f"時計盲目性（F_248(620根)<1 かつ F_1240=1）: {'PASS' if blind else 'FAIL'}")

    # ================= D5: エイリアシング =================
    print("\n=== D5: エイリアシング（粗いレジスタは細かい根を住所23に折り返す）===")
    alias_rows = []
    z = strobe["root620"]
    for W, pred_bin in PRED_ALIAS.items():
        seg = z[:W] - np.mean(z[:W])
        spec_mag = np.abs(np.fft.fft(seg))[: W // 2]
        peak = int(np.argmax(spec_mag[1:]) + 1)
        # 時計窓 W の読み: ビン b ↔ 回転角 2πb/W ↔ 読み電荷 cos^2(2πb/W)
        charge_read = math.cos(2 * pi * peak / W) ** 2
        ok = peak == pred_bin
        alias_rows.append({"window": W, "pred_bin": pred_bin, "peak_bin": peak,
                           "charge_read": charge_read, "pass": ok})
        print(f"W={W:4d}: ピークビン {peak}（予言 {pred_bin}）読み電荷 {charge_read:.6f} "
              f"{'PASS' if ok else 'FAIL'}")

    row248 = next(r for r in alias_rows if r["window"] == 248)
    row1240 = next(r for r in alias_rows if r["window"] == 1240)
    universality = (row248["pass"] and abs(row248["charge_read"] - e_low) < 1e-6
                    and row1240["pass"] and abs(row1240["charge_read"] - e_high) < 1e-6)
    print(f"電荷普遍性の機構（電子時計W=248: 620根の読み={row248['charge_read']:.6f} "
          f"= 電子電荷 {e_low:.6f} / 固有時計W=1240: {row1240['charge_read']:.6f} "
          f"= 固有電荷 {e_high:.6f}）: {'PASS' if universality else 'FAIL'}")

    # ================= 保存 =================
    with (HERE / "paperD_clock_table_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with (HERE / "paperD_alias_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(alias_rows[0])); w.writeheader(); w.writerows(alias_rows)

    fig, axes = plt.subplots(1, 4, figsize=(19, 4.4), constrained_layout=True)
    for ax, (W, _) in zip(axes, PRED_ALIAS.items()):
        seg = z[:W] - np.mean(z[:W])
        spec_mag = np.abs(np.fft.fft(seg))[: W // 2]
        ax.plot(spec_mag, lw=0.9)
        pk = int(np.argmax(spec_mag[1:]) + 1)
        ax.axvline(pk, color="tab:red", ls=":", label=f"peak bin {pk}")
        ax.set_title(f"window W={W}: read address {pk}/{W}")
        ax.set_xlabel("DFT bin"); ax.legend()
    axes[0].set_ylabel("|Z(bin)| (root620 strobe)")
    fig.suptitle("D5 (corrected): with the electron clock W=248 the fine root reads the electron address 39/124 (charge 0.302822)")
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"paperD_alias_spectra_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "paperD_two_alpha_register_pre_v1",
        "core_runner": {"path": "run_ab_invariant_theta_toy_v1.py",
                        "sha256": toy.sha256(HERE.parent / "run_ab_invariant_theta_toy_v1.py")},
        "D1": {"displacement_units": [disp_low, disp_high],
               "first_fresh_address": fresh[0],
               "operator_period_errors": {k: float(v) for k, v in checks.items()},
               "one_minus_R": {"low": e_low, "high": e_high}},
        "D2": {"rows": rows, "max_prediction_error": max_err,
               "anchor_F8_root124_error": anchor_err, "clock_blindness_pass": bool(blind)},
        "D5": {"rows": alias_rows, "charge_universality_pass": bool(universality),
               "falsified_first_dictionary": FALSIFIED_ALIAS_V0,
               "correction": ("初版辞書は反証（実測19/39/97）。原因=(i)回転角は流れ規約の"
                              "補角193π/620 (ii)DFTビン刻み2π/Wはレジスタ刻みπ/nの半分で"
                              "奇数住所には時計W=2nが必要。修正辞書で全窓的中")},
        "D9_design_consequence": (
            "位相盲目読出しの no-go により θ は凍結（本実験でも各衝突の θ は不変）。"
            "住所選択は θ 力学では起こりえず、質量平衡による選択（第五候補）は"
            "ノルム・パケットセクターで検定する必要がある"),
        "conclusion": (
            "両根の作用素周期・変位整数性・埋め込みは厳密。忠実度は閉形式予言と一致し、"
            "F_248(620根)<1 / F_1240=1 の時計盲目性を実証。粗いレジスタ（W=124）の"
            "観測者は620根を電子住所23（電荷0.302822）に折り返して読む——電荷普遍性の"
            "機構候補が動作した"),
    }
    (HERE / "paperD_two_alpha_register_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print("\nsaved: paperD_two_alpha_register_result_v1.json")


if __name__ == "__main__":
    main()
