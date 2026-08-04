#!/usr/bin/env python3
"""点火後の長時間運命 v1：「ボゾンは一瞬で全てフェルミオンに変わるか」への実測回答

問い（2026-08-04 木原氏）: 自己触媒点火（rate=C·f²）は正帰還なので、
    点火したら全ボゾンが即フェルミオン化するのではないか。

方法: v3 写像（強さ=反射率・定数ゼロ）に部分刻み積分を追加
    （非弾性の実効位相増分 r·max|amp|² を h_max=0.01 以下に分割——
     素の中点法は f>0.1 で発散するため。発散は物理でなく刻み不足と確認済み）。
    s=8・シード f₀=4.7e-3（点火圏）・3000衝突で f(j) の全曲線を記録。

判定（実行前固定）: 暴走型（f→1）/ 往復型（極大後に下降）/ プラトー型 を
    f の極大位置と終値で機械判定。
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("v3_fate", HERE / "universal_inelastic_map_v3.py")
v3 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = v3
spec.loader.exec_module(v3)
v1, toy, base = v3.v1, v3.toy, v3.base

H_MAX = 0.01
J = 3000
S = 8.0
SEED_AMP = 0.1


def collision_step_sub(a, b, sp, h_max=H_MAX):
    """v3 衝突＋非弾性の部分刻み（実効位相増分 ≤ h_max）。"""
    ro = toy.theta_from_ab(a, b, sp)
    a, b = toy.rotate_ab(a, b, ro.theta)
    r = float(ro.reflection_rate)
    if r > 0.0:
        peak = float(max(np.max(np.abs(a)) ** 2, np.max(np.abs(b)) ** 2))
        n_sub = max(1, int(math.ceil(r * peak / h_max)))
        rs = r / n_sub
        for _ in range(n_sub):
            da1, db1 = v3.vertex(a, b, rs)
            am, bm = a + 0.5 * da1, b + 0.5 * db1
            da2, db2 = v3.vertex(am, bm, rs)
            a, b = a + da2, b + db2
    return a, b, ro


def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    a = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
    a = a + v1.make_bundle(sp, v1.ODD_KS, "A", scale=SEED_AMP * S)
    b = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)
    tot0 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    c0 = abs(complex(np.sum(a * a) + np.sum(b * b)))
    fs, norms = [], []
    for _ in range(J):
        a, b, _ = collision_step_sub(a, b, sp)
        tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
        fs.append(v1.fermionic_power_raw(a, b, sp) / tot)
        norms.append(tot / tot0)
    fs = np.array(fs)
    c1 = abs(complex(np.sum(a * a) + np.sum(b * b)))
    imax = int(fs.argmax())
    f_eq = float(fs[-500:].mean())
    if fs[-1] > 0.9:
        fate = "暴走型（ほぼ全変換）"
    elif imax < J - 100 and fs[-1] < 0.9 * fs.max():
        fate = f"往復型（逆変換あり）→ 統計的平衡 f*≈{f_eq:.3f}"
    else:
        fate = f"プラトー型 f*≈{f_eq:.3f}"
    # 参考: フェルミオンマスクの位相空間割合（等分配平衡の理論値候補）
    freqs, _ = toy.combined_chi_power(a, b, sp)
    af = np.abs(freqs)
    mask_frac = float(np.mean((af >= 4) & ((af % 2) == 0)))
    print(f"点火後の運命（s={S}, f0={fs[0]:.4f}, {J}衝突, h_max={H_MAX}）")
    for j in (0, 200, 400, 500, 800, 1000, 1500, 2000, 2500, J - 1):
        print(f"  j={j:4d}: f={fs[j]:.4f} ノルム比={norms[j]:.6f}")
    print(f"  f最大={fs.max():.4f}(j={imax}) 終値={fs[-1]:.4f} 後半平均={f_eq:.4f}")
    print(f"  判定: {fate}")
    print(f"  参考: マスク位相空間割合={mask_frac:.4f}（等分配平衡の理論値候補）")
    print(f"  ノルムドリフト最大={max(abs(n-1) for n in norms):.2e} 閉塞 {c0:.1e}→{c1:.1e}")
    out = {"S": S, "seed_amp": SEED_AMP, "J": J, "h_max": H_MAX,
           "f_series_every10": [float(x) for x in fs[::10]],
           "f_max": float(fs.max()), "j_max": imax, "f_final": float(fs[-1]),
           "f_equilibrium_late500": f_eq, "mask_phase_space_fraction": mask_frac,
           "fate": fate, "norm_drift_max": float(max(abs(n - 1) for n in norms)),
           "closure_drift": c1 - c0, "runtime_sec": time.time() - t0}
    (HERE / "ignition_fate_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
