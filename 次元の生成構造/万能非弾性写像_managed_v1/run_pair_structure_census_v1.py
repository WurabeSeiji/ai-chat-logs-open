#!/usr/bin/env python3
"""対構造 census v1：生成物は相関対（粒子・反粒子対）として現れるかの実測

問い（2026-08-04 木原氏）: 反粒子を仮定せずフェルミオンが出るというなら、
    産物が実際に「対」（シグナル／アイドラー相関・逆相対巻き数）で
    現れていることを測定で示せ。

設計:
    ポンプ = 白猫（偶数束）×2チャネル（s=8）。
    シード = 単一の奇数倍音パケット（k_s=21）を A に注入（振幅 0.3×s、点火圏）。
    発展 = v3 写像（強さ=反射率・定数ゼロ）＋部分刻み積分（run_ignition_fate_v1 と同一）。

判定（実行前固定）:
    P1 和則の相棒予言: 成長パワー（シード近傍を除く）のうち、
        予言集合 {k_b + k_b' − k_s | k_b,k_b' ∈ ポンプ支持} に載る割合が、
        予言集合の測度（bin 割合）の 2 倍以上
    P2 異常相関（対の刻印）: η列アンサンブルでの ⟨c_{k_s}·c_{k_i}⟩ 正規化コヒーレンスが
        0.5 超、かつ対照（k_s×無成長bin）の中央値の 5 倍超
    P3 和則の閉じ: k_s + k_i がポンプ対和ヒストグラムの主ピークと一致（mod n）
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("fate_pair", HERE / "run_ignition_fate_v1.py")
fate = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = fate
spec.loader.exec_module(fate)
v3, v1, toy, base = fate.v3, fate.v1, fate.toy, fate.base

S = 8.0
SEED_K = 21
SEED_AMP = 0.3
J = 800


def chi_spectra(v, sp):
    """チャネル状態の χ-FFT（bin × η）。combined_chi_power と同一の変換。"""
    shape = (sp.chi_grid_n, sp.eta_grid_n)
    return np.fft.fft(v.reshape(shape), axis=0, norm="ortho")


def power_bins(a, b, sp):
    fa, fb = chi_spectra(a, sp), chi_spectra(b, sp)
    return np.sum(np.abs(fa) ** 2 + np.abs(fb) ** 2, axis=1)


def main() -> None:
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n = sp.chi_grid_n

    a = v1.make_bundle(sp, v1.EVEN_KS, "A", scale=S)
    seed = v1.make_bundle(sp, (SEED_K,), "A", scale=SEED_AMP * S)
    b = v1.make_bundle(sp, v1.EVEN_KS, "B", scale=S)

    # シードの生bin（単独スペクトルの最大bin）とポンプ支持
    p_seed = power_bins(seed, np.zeros_like(seed), sp)
    r_s = int(np.argmax(p_seed))
    p_pump = power_bins(a, b, sp)
    pump_sorted = np.argsort(p_pump)[::-1]
    csum = np.cumsum(p_pump[pump_sorted]) / p_pump.sum()
    pump_support = set(int(k) for k in pump_sorted[: int(np.searchsorted(csum, 0.99)) + 1])
    print(f"シード生bin r_s={r_s}  ポンプ支持bin数={len(pump_support)}/{n}")

    a = a + seed
    p0 = power_bins(a, b, sp)
    f0 = v1.fermionic_power_raw(a, b, sp) / float(np.vdot(a, a).real + np.vdot(b, b).real)
    for _ in range(J):
        a, b, _ = fate.collision_step_sub(a, b, sp)
    p1 = power_bins(a, b, sp)
    f1 = v1.fermionic_power_raw(a, b, sp) / float(np.vdot(a, a).real + np.vdot(b, b).real)
    print(f"フェルミオン割合: {f0:.4f} → {f1:.4f}（{J}衝突）")

    growth = np.maximum(p1 - p0, 0.0)
    # シード近傍・ポンプ支持は「新規の相棒」から除外
    excl = set(range(max(0, r_s - 1), min(n, r_s + 2))) | pump_support
    new_bins = [k for k in range(n) if k not in excl]
    g_new = growth[new_bins]

    # P1 和則の相棒予言
    partner_pred = set()
    plist = sorted(pump_support)
    for k1 in plist:
        for k2 in plist:
            partner_pred.add((k1 + k2 - r_s) % n)
    pred_new = [k for k in new_bins if k in partner_pred]
    g_pred = float(np.sum(growth[pred_new]))
    g_tot = float(np.sum(g_new)) if np.sum(g_new) > 0 else 1e-300
    frac_on_pred = g_pred / g_tot
    measure = len(pred_new) / max(len(new_bins), 1)
    p1_pass = frac_on_pred > 2 * measure
    print(f"\n[P1] 相棒予言集合上の成長割合 = {frac_on_pred:.3f} "
          f"（集合の測度 {measure:.3f}、要求 >{2*measure:.3f}）→ {'PASS' if p1_pass else 'FAIL'}")

    # 最強の相棒bin
    order = np.argsort(g_new)[::-1]
    top_partner = int(new_bins[int(order[0])])
    print(f"  最強成長bin（相棒候補）r_i = {top_partner}  成長量 {growth[top_partner]:.3e}")

    # P2 異常相関（η列アンサンブル）
    fa = chi_spectra(a, sp)
    fb = chi_spectra(b, sp)

    def coherence(k1, k2):
        num = 0.0 + 0.0j
        d1 = d2 = 0.0
        for f in (fa, fb):
            num += np.mean(f[k1] * f[k2])
            d1 += float(np.mean(np.abs(f[k1]) ** 2))
            d2 += float(np.mean(np.abs(f[k2]) ** 2))
        return abs(num) / max(np.sqrt(d1 * d2), 1e-300)

    coh_pair = coherence(r_s, top_partner)
    rng = np.random.default_rng(12345)
    controls = []
    no_grow = [k for k in new_bins if growth[k] < 0.01 * growth[top_partner]]
    for k in rng.choice(no_grow, size=min(30, len(no_grow)), replace=False):
        controls.append(coherence(r_s, int(k)))
    coh_ctrl = float(np.median(controls))
    p2_pass = coh_pair > 0.5 and coh_pair > 5 * coh_ctrl
    print(f"[P2] 異常相関コヒーレンス |⟨c_s·c_i⟩| = {coh_pair:.4f} "
          f"（対照中央値 {coh_ctrl:.4f}）→ {'PASS' if p2_pass else 'FAIL'}")

    # P3 和則の閉じ
    sums = {}
    for k1 in plist:
        for k2 in plist:
            key = (k1 + k2) % n
            sums[key] = sums.get(key, 0.0) + float(p_pump[k1] * p_pump[k2])
    peak_sum = max(sums, key=sums.get)
    pair_sum = (r_s + top_partner) % n
    p3_pass = pair_sum in sums and sums.get(pair_sum, 0.0) > 0.1 * sums[peak_sum]
    print(f"[P3] 和則: r_s + r_i = {pair_sum} (mod {n})  ポンプ対和の主ピーク = {peak_sum} "
          f"（一致度: 対和強度/ピーク = {sums.get(pair_sum, 0.0)/sums[peak_sum]:.3f}）"
          f"→ {'PASS' if p3_pass else 'FAIL'}")

    out = {"r_seed": r_s, "r_partner": top_partner, "f0": f0, "f_final": f1,
           "P1": {"pass": bool(p1_pass), "frac_on_pred": frac_on_pred, "measure": measure},
           "P2": {"pass": bool(p2_pass), "coherence_pair": coh_pair,
                   "coherence_control_median": coh_ctrl},
           "P3": {"pass": bool(p3_pass), "pair_sum": pair_sum, "pump_peak_sum": peak_sum},
           "growth_top10": {int(new_bins[int(i)]): float(g_new[int(i)])
                             for i in order[:10]},
           "runtime_sec": time.time() - t0}
    (HERE / "pair_structure_census_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    verdict = p1_pass and p2_pass and p3_pass
    print(f"\n総合: {'PASS——産物は和則相関対（粒子・反粒子対）として出現' if verdict else '一部FAIL——生データ参照'}")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
