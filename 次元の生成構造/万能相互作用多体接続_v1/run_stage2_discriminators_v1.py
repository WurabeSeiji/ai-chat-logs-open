#!/usr/bin/env python3
"""段階2 事前登録判定・残り3系統——D1 census移植／D2 E8特殊値／D3 二カイラリティ

前提: 頂点エンジン v1（run_stage2_vertex_engine_v1.py、修正済み検証 v2 ALL PASS）。
    ポンプ=スライス kp、種=スライス ks。頂点の和則（点ごと三次→k空間畳み込み）:
    パラメトリック項 w²w̄ は k_partner = 2kp − ks (mod Nreg) に相棒を生成。
    |w|²w 項は位相変調（内容の住所を変えない）。二次以降のカスケードで
    3+2-1=4 型の高次ビン、種²ポンプ̄ で DC (2·1-2=0) が到達可能。

D1 census移植（N=5, Nreg=5, kp=2, ks=1, δ=1e-2。判定は実行前固定）:
    D1-P1 和則排他性: T=30 で P_{k=3}（相棒）/P_{k=4}（カスケード）≥100
           かつ P_{k=3}/P_{k=0}（DC）≥100。
    D1-P2 対相関（位相ロック）: 生成された c₃ と予測積 c₂²c̄₁ の辺空間
           コヒーレンス ≥0.9（T=30）。二体census P2の多体版。
    D1-P3 同時性: P₃(0)=0（厳密）から立ち上がり、t∈[1,30] で単調増加
           （90%以上のstepで増加）。

D2 E8特殊値判別（レジスタ長の算術。判定は実行前固定）:
    グラフN=5固定・レジスタ長 Nreg ∈ 特殊値{120,124,128,137,144,240,248}
    と各±1近傍。kp=2, ks=1（差1→全ビン到達可能で全Nreg比較可能）。
    観測量: 点火係数 C=rate/f₀²、占有スペクトルエントロピーの半充填時間。
    判定: 各特殊値の C と半充填時間が±1近傍平均の [0.5,2] 倍以内なら null
    （線形相のnullに続きレジスタ力学でも盲目）、外れれば選択性の証拠。

D3 二カイラリティ検証（探索・方向仮説つき。pass/fail は課さない）:
    ポンプA=−枝 control親（点火する真空）／ポンプB=+枝一様安定解
    （make_parent(5,seed=2)、Kスペクトル{4,1,1,1,1}の暗い海）。
    同ノルム・同種 δ=1e-2、T=600。
    記録: 早期点火率の比（方向仮説: 読出しはパワー基底なので早期比~1）、
    相棒位相ロックの位相角（鏡映署名の探索）、長期 f_seed（方向仮説:
    −枝は自身のバーストに接続して成長持続、+枝は？——完全探索）。

使い方: python3 run_stage2_discriminators_v1.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

spec_e = importlib.util.spec_from_file_location(
    "s2v1d", HERE / "run_stage2_vertex_engine_v1.py")
s2 = importlib.util.module_from_spec(spec_e)
sys.modules[spec_e.name] = s2
spec_e.loader.exec_module(s2)
abl = s2.abl
gen3 = s2.gen3
VertexEngine = s2.VertexEngine


def slice_powers(eng):
    return np.sum(np.abs(eng.C) ** 2, axis=0)


def main() -> None:
    t0 = time.time()
    n, m = 5, 10
    out = {}
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r2 = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    vplus = r2.parent_vector / np.linalg.norm(r2.parent_vector)

    # ============ D1: census移植 ============
    print("=== D1 census移植（N=5, Nreg=5, kp=2, ks=1, δ=1e-2） ===")
    nreg = 5
    delta = 1e-2
    C0 = np.zeros((m, nreg), complex)
    C0[:, 2] = Z0c
    C0[:, 1] = delta * seed_state
    wps = {1: np.random.default_rng(92001).normal(size=m), 2: wp0.copy()}
    eng = VertexEngine(n, C0, wps, vertex_on=True)
    P_hist = [slice_powers(eng)]
    for t in range(300):
        eng.step()
        P_hist.append(slice_powers(eng))
    P_hist = np.array(P_hist)          # (T+1)×nreg
    P3_0 = P_hist[0, 3]
    T30 = 30
    p3, p4, p0 = P_hist[T30, 3], P_hist[T30, 4], P_hist[T30, 0]
    d1p1 = bool(p3 / max(p4, 1e-300) >= 100 and p3 / max(p0, 1e-300) >= 100)
    # 位相ロック: c3 vs c2²c̄1（T=30 の状態で）
    c1, c2, c3 = eng.C[:, 1], eng.C[:, 2], eng.C[:, 3]
    # T=30時点の状態を再取得（engはT=300まで走ったので再走行）
    eng2 = VertexEngine(n, C0, wps, vertex_on=True)
    for t in range(T30):
        eng2.step()
    c1, c2, c3 = eng2.C[:, 1], eng2.C[:, 2], eng2.C[:, 3]
    pred = c2 ** 2 * np.conj(c1)
    coh = abs(np.vdot(pred, c3)) / (np.linalg.norm(pred) * np.linalg.norm(c3))
    lock_phase = float(np.angle(np.vdot(pred, c3)))
    d1p2 = bool(coh >= 0.9)
    inc = np.diff(P_hist[1:T30 + 1, 3])
    mono_frac = float(np.mean(inc > 0))
    d1p3 = bool(P3_0 == 0.0 and mono_frac >= 0.9)
    out["D1"] = {"P3_T30": float(p3), "P4_T30": float(p4), "P0_T30": float(p0),
                  "ratio_P3_P4": float(p3 / max(p4, 1e-300)),
                  "ratio_P3_P0": float(p3 / max(p0, 1e-300)),
                  "lock_coherence": float(coh), "lock_phase": lock_phase,
                  "P3_initial": float(P3_0), "mono_frac": mono_frac,
                  "P1": d1p1, "P2": d1p2, "P3": d1p3}
    out["D1_pass"] = bool(d1p1 and d1p2 and d1p3)
    print(f"  P1 排他性: P3/P4={out['D1']['ratio_P3_P4']:.1e} P3/P0={out['D1']['ratio_P3_P0']:.1e} → {d1p1}")
    print(f"  P2 対相関: コヒーレンス={coh:.4f}（ロック位相={lock_phase:+.3f}rad） → {d1p2}")
    print(f"  P3 同時性: P3(0)={P3_0} 単調率={mono_frac:.2f} → {d1p3}")
    print(f"  D1 = {out['D1_pass']}")

    # ============ D2: E8特殊値判別 ============
    print("=== D2 E8特殊値判別（レジスタ長掃引・kp=2, ks=1, δ=1e-2） ===")
    specials = [120, 124, 128, 137, 144, 240, 248]
    all_nreg = sorted(set(sum([[s - 1, s, s + 1] for s in specials], [])))

    def run_reg(nreg_val, T=200):
        C0 = np.zeros((m, nreg_val), complex)
        C0[:, 2] = Z0c
        C0[:, 1] = delta * seed_state
        wps_l = {1: np.random.default_rng(92001).normal(size=m), 2: wp0.copy()}
        e = VertexEngine(n, C0, wps_l, vertex_on=True)
        f0 = e.diagnostics()["f_seed"]
        fs = []
        half_t = None
        Smax = np.log(nreg_val)
        for t in range(T):
            e.step()
            d = e.diagnostics()
            fs.append(d["f_seed"])
            if half_t is None:
                P = slice_powers(e)
                p = P / P.sum()
                ent = -np.sum(p[p > 0] * np.log(p[p > 0]))
                if ent >= 0.5 * Smax:
                    half_t = t + 1
        tt = np.arange(5, 60, dtype=float)
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, np.log(np.array(fs)[5:60]), rcond=None)
        return {"C_ign": float(coef[0] / f0 ** 2), "half_t": half_t}

    d2rows = {}
    for nreg_val in all_nreg:
        d2rows[nreg_val] = run_reg(nreg_val)
    verdicts = {}
    for s in specials:
        cs = d2rows[s]["C_ign"]
        cn = np.mean([d2rows[s - 1]["C_ign"], d2rows[s + 1]["C_ign"]])
        ratio_c = cs / cn
        hs = d2rows[s]["half_t"]
        hn = np.mean([x for x in (d2rows[s - 1]["half_t"], d2rows[s + 1]["half_t"])
                      if x is not None]) if any(
            x is not None for x in (d2rows[s - 1]["half_t"], d2rows[s + 1]["half_t"])) else None
        ratio_h = (hs / hn) if (hs is not None and hn) else None
        is_null = 0.5 <= ratio_c <= 2.0 and (ratio_h is None or 0.5 <= ratio_h <= 2.0)
        verdicts[s] = {"C_ratio": float(ratio_c),
                        "half_ratio": float(ratio_h) if ratio_h else None,
                        "null": bool(is_null)}
        print(f"  Nreg={s}: C比={ratio_c:.3f} 半充填比={ratio_h if ratio_h else '—'} "
              f"→ {'null' if is_null else '★選択性'}")
    out["D2"] = {"rows": {str(k): v_ for k, v_ in d2rows.items()},
                  "verdicts": {str(k): v_ for k, v_ in verdicts.items()}}
    out["D2_all_null"] = bool(all(v_["null"] for v_ in verdicts.values()))
    print(f"  D2 = {'全null（レジスタ力学も一次算術では盲目）' if out['D2_all_null'] else '選択性あり——要精査'}")

    # ============ D3: 二カイラリティ（探索） ============
    print("=== D3 二カイラリティ（−枝control vs +枝一様安定解、T=600） ===")
    nreg = 5
    d3 = {}
    for tag, pump in (("minus_branch", Z0c / np.linalg.norm(Z0c)),
                       ("plus_branch", vplus)):
        C0 = np.zeros((m, nreg), complex)
        C0[:, 2] = pump
        C0[:, 1] = delta * seed_state
        wps_l = {1: np.random.default_rng(92001).normal(size=m),
                 2: np.random.default_rng(96000).normal(size=m)}
        e = VertexEngine(n, C0, wps_l, vertex_on=True)
        fs = []
        for t in range(600):
            e.step()
            fs.append(e.diagnostics()["f_seed"])
        fs = np.array(fs)
        tt = np.arange(5, 40, dtype=float)
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, np.log(fs[5:40]), rcond=None)
        # 相棒位相ロック（T=30 再走行）
        e2 = VertexEngine(n, C0, wps_l, vertex_on=True)
        for t in range(30):
            e2.step()
        pred = e2.C[:, 2] ** 2 * np.conj(e2.C[:, 1])
        c3 = e2.C[:, 3]
        coh3 = abs(np.vdot(pred, c3)) / max(np.linalg.norm(pred) * np.linalg.norm(c3), 1e-300)
        ph3 = float(np.angle(np.vdot(pred, c3)))
        d3[tag] = {"early_rate": float(coef[0]), "f_seed_final": float(fs[-1]),
                    "f_seed_T300": float(fs[299]), "lock_coherence": float(coh3),
                    "lock_phase": ph3}
        print(f"  {tag}: 早期rate={coef[0]:.3e} f_seed(600)={fs[-1]:.3e} "
              f"ロック(コヒーレンス={coh3:.3f}, 位相={ph3:+.3f})")
    ratio_rate = d3["plus_branch"]["early_rate"] / d3["minus_branch"]["early_rate"]
    ratio_final = d3["plus_branch"]["f_seed_final"] / d3["minus_branch"]["f_seed_final"]
    dphase = d3["plus_branch"]["lock_phase"] - d3["minus_branch"]["lock_phase"]
    out["D3"] = {**d3, "early_rate_ratio_plus_minus": float(ratio_rate),
                  "final_ratio_plus_minus": float(ratio_final),
                  "lock_phase_diff": float(dphase)}
    print(f"  早期rate比(+/−)={ratio_rate:.3f}  最終f_seed比(+/−)={ratio_final:.3e}  "
          f"ロック位相差={dphase:+.3f}rad")

    out["runtime_sec"] = time.time() - t0
    (HERE / "stage2_discriminators_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
