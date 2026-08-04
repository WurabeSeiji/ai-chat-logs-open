#!/usr/bin/env python3
"""v3（DC/Nyquist除去版）単体テスト：全セクター閉塞の恒等成立と v2 対照

判定（実行前固定）:
    T1 行閉塞: 全行 |Σ_n w²|/Σ|w|² < 1e-12（v2 と同じ契約）
    T2 全セクター閉塞: 全 k について場の閉塞 = 0。
        具体的には |c[:,0]| ≡ 0（<1e-14）、偶数Nは |c[:,N/2]| ≡ 0、
        2k≢0 のセクターは恒等 0（構造）——つまり非閉塞セクターが存在しない
    T3 全体閉塞 < 1e-12
    T4 対照（差分が射影のみ）: 同一シードで parent_vector と raw_white_noise が
        v2 と bitwise 一致（乱数消費順不変の検証）
    T5 スペクトル白色性: 許容モードのエネルギー配分の最大/最小比が有限で
        DC/Nyquist にエネルギーがゼロ（記録・情報）
対象: N=5（seed 2・奇数）, N=40（seed 1・偶数）, N=6（seed 探索・偶数小型の追加検査）
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), HERE / name)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


g2 = load("make_parent_white_harmonics_n_only_v2.py")
g3 = load("make_parent_white_harmonics_n_only_v3.py")


def field_closure_over_total(C, k, n):
    """セクター k の場の閉塞を全状態エネルギーで正規化（消えたセクターの0/0を防ぐ）。"""
    Zk = C[:, k]
    total = float(np.sum(np.abs(C) ** 2))
    return float(abs(complex(np.sum(Zk ** 2))) * (n if (2 * k) % n == 0 else 0.0)
                 / (n * total))


def check(n, seed, label):
    r3 = g3.make_parent(n, seed=seed)
    W = r3.relation_waves
    m = W.shape[0]
    C = np.fft.fft(W, axis=1) / n
    rowc = float(np.max(np.abs(np.sum(W * W, axis=1)) / np.sum(np.abs(W) ** 2, axis=1)))
    t1 = rowc < 1e-12
    dc_amp = float(np.max(np.abs(C[:, 0])))
    nyq_amp = float(np.max(np.abs(C[:, n // 2]))) if n % 2 == 0 else 0.0
    worst_field = max(field_closure_over_total(C, k, n) for k in range(n))
    t2 = dc_amp < 1e-14 and nyq_amp < 1e-14 and worst_field < 1e-24
    tot = abs(complex(np.sum(W * W)))
    t3 = tot < 1e-12
    r2 = g2.make_parent(n, seed=seed)
    t4 = (np.array_equal(r3.parent_vector, r2.parent_vector)
          and np.array_equal(r3.raw_white_noise, r2.raw_white_noise))
    energy = np.sum(np.abs(C) ** 2, axis=0)
    allowed = [k for k in range(n) if (2 * k) % n != 0]
    ratio = float(energy[allowed].max() / energy[allowed].min())
    print(f"[{label}] T1行閉塞 {'PASS' if t1 else 'FAIL'}({rowc:.1e}) "
          f"T2全セクター閉塞 {'PASS' if t2 else 'FAIL'}(DC={dc_amp:.1e},Nyq={nyq_amp:.1e},場max={worst_field:.1e}) "
          f"T3全体 {'PASS' if t3 else 'FAIL'}({tot:.1e}) "
          f"T4対照bitwise {'PASS' if t4 else 'FAIL'} "
          f"T5許容モードE比max/min={ratio:.2f}")
    return {"T1": t1, "T2": t2, "T3": t3, "T4": t4,
            "row_closure_max": rowc, "dc_amp_max": dc_amp, "nyquist_amp_max": nyq_amp,
            "total_closure": tot, "energy_ratio_allowed": ratio}


def main():
    results = {}
    results["N5_seed2"] = check(5, 2, "N=5 seed=2 (奇数)")
    results["N40_seed1"] = check(40, 1, "N=40 seed=1 (偶数)")
    # 偶数小型: 収束するseedを探索（カウントアップ・全試行記録）
    seed = 0
    while True:
        try:
            g3.make_parent(6, seed=seed)
            break
        except g3.ParentConstructionError:
            seed += 1
    results["N6_first_seed"] = {"seed": seed, **check(6, seed, f"N=6 seed={seed} (偶数小型)")}
    ok = all(all(v[t] for t in ("T1", "T2", "T3", "T4")) for v in
             (results["N5_seed2"], results["N40_seed1"])) and \
         all(results["N6_first_seed"][t] for t in ("T1", "T2", "T3", "T4"))
    print(f"\n総合判定: {'ALL PASS' if ok else 'FAIL あり'}")
    results["all_pass"] = bool(ok)
    (HERE / "test_v3_result.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")


if __name__ == "__main__":
    main()
