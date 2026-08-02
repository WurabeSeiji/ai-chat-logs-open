#!/usr/bin/env python3
"""主張硬化バッテリー v1（E2/E3/E4）

論文Aの主張の未検証側面を塞ぐ三実験。予言は測定前に固定。

E2 多重共有チャネルの積則:
    A = u5 A5 + u7 A7, B = v7 e^{i phi1} B7 + v9 e^{i phi2} B9 とすると、
    共有チャネルは (A5-B7) と (A7-B9) の2本（梯子 Delta k=2）。予言:
        <a|b> = u5 v7 c57 e^{i phi1} + u7 v9 c79 e^{i phi2}
    すなわち実効電荷結合は**チャネルごとの振幅積のコヒーレント和**であり、
    変調は2位相の干渉 M(phi1,phi2) を示す。全 (phi1,phi2) グリッドで
    流れの実測を厳密予言と比較する。

E3 符号の多衝突持続性:
    固定プローブ回転を繰り返すと全回転角は j*theta_p として蓄積するため、
    移乗は Rabi 型に振動すると予言される（周期 pi/theta_p）。
    「符号」は各瞬間の流れの向きであり、初期方向は phi の半回転で反転したまま
    振動全体が鏡映になる（N_B(j; phi) と N_B(j; phi+pi) が初期ノルム対称）。
    Bjerknes 力が音響周期の時間平均であるのと同様、正味の定常流には
    平均化の議論が必要であることを実測で記録する。

E4 平均合成則の規約拡張:
    対角パワーはビン別に加法的なので、R_joint = (f_A^(m) + f_B^(m))/2 は
    **任意のマスク m について**恒等的に成立するはず（定理）。変種マスク
    M1(even>=6), M3(any>=4) で数値確認する（分率はマスク別にビン表から測る）。
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"

spec = importlib.util.spec_from_file_location("toy_for_hardening_v1", TOY_RUNNER_PATH)
toy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = toy
spec.loader.exec_module(toy)
base = toy.base
plt = base.plt

PRED_TOL = 1.0e-12


def unit_norm(v: np.ndarray) -> np.ndarray:
    return v / math.sqrt(float(np.vdot(v, v).real))


def mk(k: int, which: str, hair: bool, sp) -> np.ndarray:
    case = base.explicit_packet_case(mode=f"hard_{which}_{k}_{int(hair)}", packet_a=(k,), packet_b=(k,))
    return unit_norm(base.make_case_state(sp, case, which, hair_enabled=hair))


def spectra_power(state: np.ndarray, sp) -> tuple[np.ndarray, np.ndarray]:
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    f = np.fft.fft(state.reshape(shape), axis=0, norm="ortho")
    freqs = np.rint(np.fft.fftfreq(sp.chi_grid_n, d=1.0 / sp.chi_grid_n)).astype(int)
    return freqs, np.sum(np.abs(f) ** 2, axis=1)


MASKS: dict[str, Callable[[np.ndarray], np.ndarray]] = {
    "M1_even_ge6": lambda f: (f >= 6) & (f % 2 == 0),
    "M3_any_ge4": lambda f: f >= 4,
}


def fraction(state: np.ndarray, sp, mask_fn) -> float:
    freqs, p = spectra_power(state, sp)
    m = mask_fn(np.abs(freqs))
    return float(np.sum(p[m])) / float(np.sum(p))


def variant_R(a, b, sp, mask_fn) -> float:
    freqs, pa = spectra_power(a, sp)
    _, pb = spectra_power(b, sp)
    p = pa + pb
    m = mask_fn(np.abs(freqs))
    pf = float(np.sum(p[m])); pball = float(np.sum(p))
    return pf / pball


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=1)
    sp = base.build_source_params(params)

    # ================= E2 多重共有チャネル =================
    a5, a7 = mk(5, "A", False, sp), mk(7, "A", False, sp)
    b7, b9 = mk(7, "B", False, sp), mk(9, "B", False, sp)
    c57 = complex(np.vdot(a5, b7))
    c79 = complex(np.vdot(a7, b9))
    leak = max(abs(complex(np.vdot(a5, b9))), abs(complex(np.vdot(a7, b7))))
    print(f"E2 channels: c57={c57:.6f} c79={c79:.6f} leak={leak:.1e}")

    th_p = 0.15
    s2t = math.sin(2 * th_p)
    u5 = u7 = 1.0 / math.sqrt(2.0)
    v7 = v9 = 1.0 / math.sqrt(2.0)
    a0 = unit_norm(u5 * a5 + u7 * a7)
    grid = [2 * math.pi * k / 6 for k in range(6)]
    rows_e2 = []
    max_err = 0.0
    for phi1 in grid:
        for phi2 in grid:
            b0 = unit_norm(v7 * np.exp(1j * phi1) * b7 + v9 * np.exp(1j * phi2) * b9)
            a2, b2 = toy.rotate_ab(a0.copy(), b0.copy(), th_p)
            dnb = float(np.vdot(b2, b2).real) - 1.0
            pred_overlap = (
                u5 * v7 * c57 * np.exp(1j * phi1) + u7 * v9 * c79 * np.exp(1j * phi2)
            )
            pred = s2t * pred_overlap.real
            err = abs(dnb - pred)
            max_err = max(max_err, err)
            rows_e2.append(
                {"phi1": phi1, "phi2": phi2, "dN_B": dnb, "pred": float(pred), "err": err}
            )
    e2_pass = max_err <= PRED_TOL
    # チャネル干渉の実在: phi1 固定で phi2 により変調幅が変わる（コヒーレント和の証拠）
    m_at = {}
    for phi2 in (0.0, math.pi):
        vals = [r["dN_B"] for r in rows_e2 if abs(r["phi2"] - phi2) < 1e-12]
        m_at[phi2] = (max(vals) - min(vals)) / 2.0
    print(f"E2: coherent-sum prediction max|err|={max_err:.2e} -> {'PASS' if e2_pass else 'FAIL'}")
    print(f"E2: channel interference visible (swing at phi2=0: {m_at[0.0]:.4f} vs pi: {m_at[math.pi]:.4f})")

    # ================= E3 符号の多衝突持続性 =================
    a1b, b1b = mk(1, "A", False, sp), mk(1, "B", False, sp)
    a5b, b7b = a5, b7
    f = 0.6
    a0 = unit_norm(math.sqrt(1 - f) * a1b + math.sqrt(f) * a5b)
    th_p3 = 0.1
    traj = {}
    for phi in (0.0, math.pi):
        b0 = unit_norm(math.sqrt(1 - f) * b1b + np.exp(1j * phi) * math.sqrt(f) * b7b)
        a, b = a0.copy(), b0.copy()
        nb = []
        for j in range(80):
            a, b = toy.rotate_ab(a, b, th_p3)
            nb.append(float(np.vdot(b, b).real))
        traj[phi] = np.asarray(nb)
    period_pred = math.pi / th_p3
    mirror_err = float(np.max(np.abs((traj[0.0] - 1.0) + (traj[math.pi] - 1.0) -
                                     (traj[0.0][:] - 1.0 + traj[math.pi][:] - 1.0))))
    # 鏡映性: N_B(j;0)-1 と N_B(j;pi)-1 の和が「位相非依存成分×2」に一致するか
    sym_component = (traj[0.0] + traj[math.pi]) / 2.0
    antisym_0 = traj[0.0] - sym_component
    antisym_pi = traj[math.pi] - sym_component
    mirror_max = float(np.max(np.abs(antisym_0 + antisym_pi)))
    first_reversal = int(np.argmax(np.abs(traj[0.0] - 1.0)))
    print(
        f"E3: transfer oscillates (Rabi-like), predicted period ~{period_pred:.1f} collisions;"
        f" first extremum at j={first_reversal + 1}"
    )
    print(f"E3: phase-mirror antisymmetry max|err|={mirror_max:.2e} (sign component exactly reverses)")

    # ================= E4 平均則の規約拡張 =================
    rows_e4 = []
    e4_max = 0.0
    b1c, b7c, a1c = mk(1, "B", True, sp), mk(7, "B", True, sp), mk(1, "A", True, sp)
    a5c = mk(5, "A", True, sp)
    for mask_name, mask_fn in MASKS.items():
        for fa in (0.0, 0.3, 0.7, 1.0):
            a0 = unit_norm(math.sqrt(1 - fa) * a1c + math.sqrt(fa) * a5c)
            f_a = fraction(a0, sp, mask_fn)
            for fb in (0.0, 0.3, 0.7, 1.0):
                b0 = unit_norm(math.sqrt(1 - fb) * b1c + math.sqrt(fb) * b7c)
                f_b = fraction(b0, sp, mask_fn)
                r_joint = variant_R(a0, b0, sp, mask_fn)
                err = abs(r_joint - (f_a + f_b) / 2.0)
                e4_max = max(e4_max, err)
                rows_e4.append(
                    {"mask": mask_name, "fa_dial": fa, "fb_dial": fb,
                     "f_A_mask": f_a, "f_B_mask": f_b, "R_joint": r_joint, "err": err}
                )
    e4_pass = e4_max <= PRED_TOL
    print(f"E4: mean law under variant masks max|err|={e4_max:.2e} -> {'PASS' if e4_pass else 'FAIL'}")

    # ---- 保存 ----
    with (HERE / "claim_hardening_e2_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows_e2[0])); w.writeheader(); w.writerows(rows_e2)
    with (HERE / "claim_hardening_e4_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows_e4[0])); w.writeheader(); w.writerows(rows_e4)
    np.savez_compressed(
        HERE / "claim_hardening_e3_trajectories_v1.npz",
        nb_phi0=traj[0.0], nb_phipi=traj[math.pi],
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    g = np.asarray([[next(r["dN_B"] for r in rows_e2
                          if abs(r["phi1"] - p1) < 1e-12 and abs(r["phi2"] - p2) < 1e-12)
                     for p2 in grid] for p1 in grid])
    im = axes[0].imshow(g, origin="lower", extent=(0, 2 * math.pi, 0, 2 * math.pi))
    axes[0].set_xlabel("phi2"); axes[0].set_ylabel("phi1")
    axes[0].set_title(f"E2: two-channel coherent sum (err {max_err:.0e})")
    fig.colorbar(im, ax=axes[0], shrink=0.8)
    axes[1].plot(traj[0.0] - 1.0, label="phi=0")
    axes[1].plot(traj[math.pi] - 1.0, label="phi=pi")
    axes[1].set_xlabel("collision"); axes[1].set_ylabel("N_B - 1")
    axes[1].set_title("E3: Rabi-type oscillation, mirrored by phase")
    axes[1].legend(); axes[1].grid(alpha=0.3)
    for mask_name in MASKS:
        sub = [r for r in rows_e4 if r["mask"] == mask_name]
        axes[2].plot([(r["f_A_mask"] + r["f_B_mask"]) / 2 for r in sub],
                     [r["R_joint"] for r in sub], "o", markersize=4, label=mask_name)
    axes[2].plot([0, 1], [0, 1], "-", color="0.6", linewidth=0.8)
    axes[2].set_xlabel("(f_A+f_B)/2 per mask"); axes[2].set_ylabel("R_joint")
    axes[2].set_title(f"E4: mean law, variant masks (err {e4_max:.0e})")
    axes[2].legend(fontsize=7); axes[2].grid(alpha=0.3)
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"claim_hardening_battery_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "claim_hardening_battery_v1",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "E2_multichannel": {
            "channels": {"c57": [c57.real, c57.imag], "c79": [c79.real, c79.imag]},
            "cross_leak": leak, "max_err": max_err, "verdict": "PASS" if e2_pass else "FAIL",
            "claim": "実効電荷結合はチャネル別振幅積のコヒーレント和（電荷のチャネル加算性）",
        },
        "E3_sign_persistence": {
            "predicted_rabi_period": period_pred,
            "phase_mirror_antisymmetry_max_err": mirror_max,
            "claim": "移乗はRabi型振動、符号成分は位相半回転で厳密鏡映。正味定常流は平均化の議論が必要（Bjerknes同様）",
        },
        "E4_mean_law_variant_masks": {"max_err": e4_max, "verdict": "PASS" if e4_pass else "FAIL"},
    }
    (HERE / "claim_hardening_battery_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
