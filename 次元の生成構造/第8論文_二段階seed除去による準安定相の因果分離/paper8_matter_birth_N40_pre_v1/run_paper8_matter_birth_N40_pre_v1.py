#!/usr/bin/env python3
"""第8論文v2予備実験 E-M3：分解能閾値仮説の判別（N=40）v1

対立仮説（測定前固定・判別実験）:
    仮説A（公理2説）: 円環性（周期的位相軸）がない限り、Nを上げても
        周波数は整数格子に乗らない。N=40でもユニゾン（比ずれ小）または
        非通約のまま、非自明整数比ロック L=0。
    仮説B（分解能閾値説・木原氏）: N=5(M=10)は位数124級の住所を載せる
        分解能が不足していただけ。N=40(M=780>124)では周波数が分化し、
        有理ロックが現れ始める。L>0 が持続的に出る。

観測量: E-M1と同一（辺位相の窓内傾き→非自明整数比ロック数L(t)、
        持続5窓、tol=1e-3、比≥2）＋ユニゾン診断（比ずれ最大）。

判定:
    ユニゾン（比ずれ<1e-2）かつ L=0 → 仮説A支持
    周波数分化（比ずれ>1e-1）かつ L>0 持続 → 仮説B支持
    分化するが整数比に乗らない（比ずれ大・L=0）→ 中間: 分化には分解能が
        効くが、整数格子化には円環性が必要（A/Bの合成）
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CODE = HERE.parent / "code" / "run_preliminary_seed_ablation_v1.py"
spec = importlib.util.spec_from_file_location("ablation_for_N40_v1", CODE)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)

N = 40
XMAX = 20000
WIN = 200
STRIDE = 100
RATIO_TOL = 1e-3
PERSIST = 5
FREQ_MIN = 1e-6


def main() -> None:
    t_start = time.time()
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = abl.build_init(N, initial_seed=False)
    M = sys_lr.m
    print(f"N={N}, M={M}, 構築 {time.time()-t_start:.1f}s")

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    phases = np.zeros((XMAX + 1, M), dtype=np.float32)
    fs = np.zeros(XMAX + 1)
    phases[0] = np.angle(Z); fs[0] = fval(Z)
    crossing = None
    for t in range(1, XMAX + 1):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z)
        fs[t] = fval(Z)
        if crossing is None and fs[t] > 0.05:
            crossing = t
        if t % 5000 == 0:
            print(f"  step {t} ({time.time()-t_start:.0f}s) f={fs[t]:.3e}")
    meta_start = crossing + abl.GUARD if crossing is not None else None
    print(f"crossing = {crossing}, 準安定開始 = {meta_start}")

    unwrapped = np.unwrap(phases.astype(np.float64), axis=0)
    centers, freqs = [], []
    x = np.arange(WIN)
    for s in range(0, XMAX - WIN, STRIDE):
        slope = np.polyfit(x, unwrapped[s:s + WIN], 1)[0]
        centers.append(s + WIN // 2)
        freqs.append(slope)
    centers = np.array(centers); freqs = np.abs(np.array(freqs))

    # ロック数（ベクトル化: 各窓で比行列→整数近接）
    L = np.zeros(len(centers), dtype=int)
    lock_sets = []
    for w in range(len(centers)):
        f = freqs[w]
        valid = f > FREQ_MIN
        fv = f[valid]
        if fv.size < 2:
            lock_sets.append(frozenset()); continue
        r = fv[:, None] / fv[None, :]
        hi = np.maximum(r, 1.0 / np.maximum(r, 1e-30))
        pr = np.round(hi)
        mask = (pr >= 2) & (np.abs(hi - pr) < RATIO_TOL)
        iu = np.triu_indices(fv.size, k=1)
        pairs = frozenset(
            (int(i), int(j), int(pr[i, j])) for i, j in zip(*iu) if mask[i, j])
        lock_sets.append(pairs)
    for w in range(len(lock_sets)):
        cnt = 0
        for key in lock_sets[w]:
            run = 1
            k = w - 1
            while k >= 0 and key in lock_sets[k]:
                run += 1; k -= 1
            k = w + 1
            while k < len(lock_sets) and key in lock_sets[k]:
                run += 1; k += 1
            if run >= PERSIST:
                cnt += 1
        L[w] = cnt

    meta_mask = centers > (meta_start or XMAX)
    tail = freqs[-40:]
    fbar = np.mean(tail, axis=0)
    fb = fbar[fbar > FREQ_MIN]
    rmat = fb[:, None] / np.maximum(fb[None, :], 1e-30)
    max_dev = float(np.max(np.abs(rmat[rmat >= 1] - 1))) if fb.size > 1 else 0.0
    Lmeta = int(np.max(L[meta_mask])) if np.any(meta_mask) else 0
    print(f"準安定末期: 周波数比ずれ最大 = {max_dev:.3e} / 準安定期の最大ロック数 L = {Lmeta}")

    if max_dev < 1e-2 and Lmeta == 0:
        verdict = "仮説A支持（ユニゾン継続・ロックなし——円環性の不在が原因、分解能でない）"
    elif Lmeta > 0:
        verdict = "仮説B支持（分解能増で有理ロック出現——物質誕生に分解能閾値）"
    else:
        verdict = "中間（周波数は分化するが整数比に乗らない——分化は分解能、整数格子化は円環性）"
    print(f"判定: {verdict}")

    payload = {
        "experiment": "paper8_matter_birth_N40_pre_v1",
        "N": N, "M": int(M), "XMAX": XMAX,
        "crossing": crossing, "metastable_start": meta_start,
        "metastable_max_ratio_deviation": max_dev,
        "metastable_max_locks": Lmeta,
        "verdict": verdict,
        "runtime_sec": time.time() - t_start,
    }
    (HERE / "paper8_matter_birth_N40_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({time.time()-t_start:.0f}s)")


if __name__ == "__main__":
    main()
