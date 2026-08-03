#!/usr/bin/env python3
"""E-M3補遺: N=300 判別（辺800本サブサンプル、実行済み・結果はJSON参照）
判定基準・観測量は run_paper8_matter_birth_N40_pre_v1.py と同一。
実測: crossing=4849, 準安定開始=7849, 比ずれ最大=3.674e-3, 非自明整数比対=0
→ 仮説A支持（ユニゾン・ロックなし）
"""
import importlib.util, sys, time, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CODE = HERE.parent / "code" / "run_preliminary_seed_ablation_v1.py"
spec = importlib.util.spec_from_file_location("abl_n300", CODE)
abl = importlib.util.module_from_spec(spec); sys.modules[spec.name] = abl
spec.loader.exec_module(abl)

N, XMAX = 300, 12000

def main():
    t0 = time.time()
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = abl.build_init(N, initial_seed=False)
    M = sys_lr.m
    rng = np.random.default_rng(7)
    sub = np.sort(rng.choice(M, size=800, replace=False))

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    phases = np.zeros((XMAX + 1, 800), dtype=np.float32)
    fs = np.zeros(XMAX + 1)
    phases[0] = np.angle(Z[sub]); fs[0] = fval(Z)
    crossing = None
    for t in range(1, XMAX + 1):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z[sub]); fs[t] = fval(Z)
        if crossing is None and fs[t] > 0.05:
            crossing = t
    meta = crossing + abl.GUARD if crossing else None

    u = np.unwrap(phases.astype(np.float64), axis=0)
    seg = u[XMAX - 4000:]
    fr = np.polyfit(np.arange(seg.shape[0]), seg, 1)[0]
    fb = np.abs(fr); fb = fb[fb > 1e-6]
    r = fb[:, None] / np.maximum(fb[None, :], 1e-30)
    rmax = np.maximum(r, 1 / np.maximum(r, 1e-30))
    max_dev = float(np.max(np.abs(r[r >= 1] - 1)))
    pr = np.round(rmax)
    locks = int(np.sum((pr >= 2) & (np.abs(rmax - pr) < 1e-3)) // 2)
    verdict = ("仮説A支持（ユニゾン・ロックなし）" if max_dev < 1e-2 and locks == 0
               else ("仮説B支持（有理ロック出現）" if locks > 0 else "中間"))
    print(f"crossing={crossing} meta={meta} 比ずれ最大={max_dev:.3e} ロック={locks} → {verdict}")
    json.dump({"experiment": "paper8_matter_birth_N300_supplement_v1", "N": N, "M": int(M),
               "XMAX": XMAX, "subsample": 800, "crossing": crossing, "metastable_start": meta,
               "max_ratio_deviation": max_dev, "nontrivial_locks": locks, "verdict": verdict,
               "runtime_sec": time.time() - t0},
              open(HERE / "paper8_matter_birth_N300_result_v1.json", "w"),
              ensure_ascii=False, indent=2)

if __name__ == "__main__":
    main()
