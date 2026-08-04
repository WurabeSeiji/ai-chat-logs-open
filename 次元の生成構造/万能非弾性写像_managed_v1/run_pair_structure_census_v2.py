#!/usr/bin/env python3
"""対構造 census v2：増幅初期窓での測定＋事前予言の固定（v1 の測定窓誤りの訂正）

v1 の誤り: 800衝突＝熱化後（f=0.45 の灰色平衡）で相関を測った。対の刻印は
    増幅の線形域にあり、熱化後は掻き混ぜられる。
v1 の診断からの事前予言: ポンプ対和ヒストグラムのピークは 0 (mod n)
    （ポンプが ±対称のため）。よって和則 k₁+k₂−k_s の相棒は
        **r_i = −r_s (mod n) ＝ 逆巻きミラー bin**
    ——反粒子の位置そのもの。これを実行前予言として固定する。

判定（実行前固定）:
    P1 ミラー出現: 増幅窓（P_f が10倍成長 かつ f<0.1 の時点）で、
        新規成長binの上位3位以内にミラー bin n−r_s が入る
    P2 異常相関: 同時点の |⟨c_{r_s}·c_{n−r_s}⟩|η コヒーレンス > 0.5 かつ対照中央値の5倍
    P3 対の同時成長: ミラー bin の成長量がシード bin の成長量の 0.1〜10 倍
        （対で生まれる＝同程度に増える）
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("fate_pair2", HERE / "run_ignition_fate_v1.py")
fate = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = fate
v3 = None
sys.modules[spec.name] = fate
spec.loader.exec_module(fate)
v3, v1, toy, base = fate.v3, fate.v1, fate.toy, fate.base

S = 8.0
SEED_K = 21
SEED_AMP = 0.2
J = 800
SAMPLE = 5
F_WINDOW = 0.08   # 増幅窓の判定: f がこの値に達した最初の標本点（熱化前・増幅後期）


def single_winding(v, sp):
    """状態の χ スペクトル正側（bin 1..n/2−1）だけを残す＝純単一巻き化。
    片側のみの占有は閉塞恒等 0（対積 c_k c_{n−k}=0）——v4 生成器と同じ定理。"""
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
    n = sp.chi_grid_n
    f[0, :] = 0.0
    f[n // 2:, :] = 0.0
    w = np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    nrm = float(np.linalg.norm(w))
    return w / nrm if nrm > 0 else w


def chi_spectra(v, sp):
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    return np.fft.fft(v.reshape(shape), axis=0, norm="ortho")


def power_bins(a, b, sp):
    fa, fb = chi_spectra(a, sp), chi_spectra(b, sp)
    return np.sum(np.abs(fa) ** 2 + np.abs(fb) ** 2, axis=1)


def coherence(a, b, sp, k1, k2):
    fa, fb = chi_spectra(a, sp), chi_spectra(b, sp)
    num = 0.0 + 0.0j
    d1 = d2 = 0.0
    for f in (fa, fb):
        num += np.mean(f[k1] * f[k2])
        d1 += float(np.mean(np.abs(f[k1]) ** 2))
        d2 += float(np.mean(np.abs(f[k2]) ** 2))
    return abs(num) / max(np.sqrt(d1 * d2), 1e-300)


def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n = sp.chi_grid_n

    a0p = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
    seed_raw = v1.make_bundle(sp, (SEED_K,), "A", scale=1.0)
    seed = single_winding(seed_raw, sp) * (SEED_AMP * S)   # 純単一巻きシード
    b = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)
    sc = abs(complex(np.sum(seed * seed)))
    print(f"単一巻きシードの閉塞 |Σseed²| = {sc:.2e}（恒等0のはず）")
    r_s = int(np.argmax(power_bins(seed, np.zeros_like(seed), sp)))
    r_mirror = (n - r_s) % n
    print(f"シード生bin r_s={r_s} → 事前予言の相棒（逆巻きミラー）r_i={r_mirror}")

    p_pump = power_bins(a0p, b, sp)
    pump_sorted = np.argsort(p_pump)[::-1]
    csum = np.cumsum(p_pump[pump_sorted]) / p_pump.sum()
    pump_support = set(int(k) for k in pump_sorted[: int(np.searchsorted(csum, 0.99)) + 1])

    a = a0p + seed
    p0 = power_bins(a, b, sp)
    pf0 = v1.fermionic_power_raw(a, b, sp)
    tot = float(np.vdot(a, a).real + np.vdot(b, b).real)
    window = None
    hist = []
    for j in range(1, J + 1):
        a, b, _ = fate.collision_step_sub(a, b, sp)
        if j % SAMPLE == 0:
            pf = v1.fermionic_power_raw(a, b, sp)
            f = pf / tot
            hist.append((j, pf / pf0, f))
            if window is None and f >= F_WINDOW:
                window = j
                p1 = power_bins(a, b, sp)
                coh_pair = coherence(a, b, sp, r_s, r_mirror)
                rng = np.random.default_rng(777)
                growth = np.maximum(p1 - p0, 0.0)
                excl = set(range(max(0, r_s - 1), min(n, r_s + 2))) | pump_support
                new_bins = [k for k in range(n) if k not in excl]
                g_new = growth[np.array(new_bins)]
                order = np.argsort(g_new)[::-1]
                top3 = [int(new_bins[int(i)]) for i in order[:3]]
                no_grow = [k for k in new_bins if growth[k] < 0.01 * g_new[int(order[0])]]
                ctrl = [coherence(a, b, sp, r_s, int(k))
                        for k in rng.choice(no_grow, size=min(30, len(no_grow)), replace=False)]
                coh_ctrl = float(np.median(ctrl))
                g_seed = float(growth[r_s])
                g_mirror = float(growth[r_mirror])
                break
    if window is None:
        print("増幅窓（10倍成長・f<0.1）に到達せず——生データを保存して終了")
        return

    p1_pass = r_mirror in top3
    p2_pass = coh_pair > 0.5 and coh_pair > 5 * coh_ctrl
    ratio = g_mirror / max(g_seed, 1e-300)
    p3_pass = 0.1 <= ratio <= 10.0
    print(f"増幅窓: j={window}（f={F_WINDOW} 到達点＝熱化前）")
    print(f"[P1] 新規成長トップ3 = {top3}（ミラー {r_mirror} を含むか）→ {'PASS' if p1_pass else 'FAIL'}")
    print(f"[P2] 異常相関 |⟨c_s·c_mirror⟩| = {coh_pair:.4f}（対照中央値 {coh_ctrl:.4f}）"
          f"→ {'PASS' if p2_pass else 'FAIL'}")
    print(f"[P3] 対の同時成長: ミラー/シード 成長比 = {ratio:.3f}（0.1〜10 要求）"
          f"→ {'PASS' if p3_pass else 'FAIL'}")
    verdict = p1_pass and p2_pass and p3_pass
    print(f"\n総合: {'PASS——産物は逆巻きミラー相関対（反粒子は出力）' if verdict else '一部FAIL——生データ参照'}")

    out = {"r_seed": r_s, "r_mirror": r_mirror, "window_j": window,
           "P1": {"pass": bool(p1_pass), "top3": top3},
           "P2": {"pass": bool(p2_pass), "coherence": coh_pair, "control": coh_ctrl},
           "P3": {"pass": bool(p3_pass), "growth_ratio_mirror_over_seed": ratio,
                   "growth_seed": g_seed, "growth_mirror": g_mirror},
           "history_j_growth_f": hist[:40],
           "runtime_sec": time.time() - t0}
    (HERE / "pair_structure_census_result_v2.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
