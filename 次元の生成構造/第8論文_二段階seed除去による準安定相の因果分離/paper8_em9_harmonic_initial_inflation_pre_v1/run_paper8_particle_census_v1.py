#!/usr/bin/env python3
"""粒子的グループ census v1：N=5 / N=40 の海の粒子構成の分類台帳

分類体系（2026-08-04 整理・木原氏承認）:
    第0層 塊の同定: 振幅偏差重み付き巻き数スペクトル A_k = |Σw e^{ikθ}|/Σw の
        ピーク（>0.5）→ 波長 1/k のグループ。成分帰属は偏差閾値
        |dev| > max(1e-3, 0.1·max|dev|)。等振幅状態は塊なし＝一様海。
    第1層 統計性: 行動的半周期テスト（E-M11 倍音時計 micro=段n、72 tick 進めた
        重なりの実部符号。偶=ボゾン型/奇=フェルミオン型）＋巻き数偶奇。
    第2層 倍音構成と位相: (段 n, 巻き数 k)、フレーム方位 arg(R_k)/k、
        同一状態内の塊間相対角（ゲージ共有なので物理量）。
    第3層 束縛・電荷: 塊の周波数と海（残余）の周波数の生の比を記録
        （一定にならなくてもそのまま記録——木原氏指示）。整数比ロックは
        E-M4 基準で判定。
    第4層 質量的指標: 塊の振幅割れ RMS(dev)（辞書未確定と明記）。
    第5層 周回符号（毛の類似量・相対のみ）: s = Σ_S dev·sin(2k(θ−φ))。
        同一状態内の相対符号のみ物理。同 k・逆符号の塊＝対候補。
        絶対ラベル（粒子/反粒子）は規約であり付与しない。

対象: N=5（対照＋倍音海8段）、N=40（対照＋4段）。緩和 6000 tick 後に判定。
規約: 種ラベル IF なし・read-only import・SHA-256 記録・空欄も記録。
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
MPH = REPO / "次元の生成構造" / "make_parent_harmonic_unit_v1" / "make_parent_harmonic_v1.py"

spec = importlib.util.spec_from_file_location("abl_cen", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("mph_cen", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

T_RELAX = 6000
WIN = 500
T_HALF = 72
K_MAX = 12
A_PEAK = 0.5
TOL_LOCK = 1e-3
FREQ_MIN = 1e-8

CFG = {5: {"H": 8, "seed": 40260801},
       40: {"H": 4, "seed": 40260802}}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def census_state(n, level, v0, wp, label):
    """1状態の census。level=倍音段（対照は1）。"""
    micro = level
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(v0))
    m = sys_lr.m
    Z = v0.copy()
    for _ in range(T_RELAX):
        for _ in range(micro):
            Z, wp = abl.evolve(sys_lr, Z, wp)
    Zref = Z.copy()
    phases = np.zeros((WIN, m))
    for t in range(WIN):
        for _ in range(micro):
            Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z)
    # 半周期状態（緩和点から T_HALF tick）
    Zh = Zref.copy(); wp2 = wp.copy()
    for _ in range(T_HALF):
        for _ in range(micro):
            Zh, wp2 = abl.evolve(sys_lr, Zh, wp2)

    amp = np.abs(Zref)
    th = np.angle(Zref)
    med = float(np.median(amp))
    dev = amp - med
    w = np.abs(dev)
    freqs = np.abs(np.polyfit(np.arange(WIN), np.unwrap(phases, axis=0), 1)[0])

    state = {"label": label, "level": level, "M": m,
             "amp_range": [float(amp.min()), float(amp.max())],
             "amp_split_rms_total": float(np.std(amp)),
             "omega_mean": float(np.mean(freqs))}

    # 全体の行動的パリティ
    ov_all = complex(np.vdot(Zref, Zh))
    state["parity_behavioral"] = "odd" if ov_all.real < 0 else "even"
    state["parity_overlap_re"] = float(ov_all.real)

    # 塊の同定
    clumps = []
    uniform = w.sum() < 1e-8 * max(1.0, amp.sum())
    k_hi = max(1, min(K_MAX, m // 2))          # エイリアス排除: k ≤ M/2
    if not uniform:
        Ak = {k: abs(np.sum(w * np.exp(1j * k * th))) / w.sum()
              for k in range(1, k_hi + 1)}
        peaks = [k for k, a in Ak.items() if a > A_PEAK]
        state["winding_spectrum"] = {str(k): float(a) for k, a in Ak.items()}
        base_sel = w > max(1e-3, 0.1 * w.max())
        for k in peaks:
            zc = np.sum(w[base_sel] * np.exp(1j * k * th[base_sel]))
            c = float(np.angle(zc))            # k·θ の残差中心
            phi = c / k
            resid = np.abs(np.angle(np.exp(1j * (k * th - c))))
            member_mask = base_sel & (resid < 0.8)
            S = np.nonzero(member_mask)[0]
            if S.size < 2:
                continue
            rest = ~member_mask
            energy = float(np.sum(amp[S] ** 2))
            ovS = complex(np.sum(np.conj(Zref[S]) * Zh[S]))
            fS = freqs[S]
            fsea = freqs[rest] if rest.any() else np.array([])
            if fsea.size and np.mean(fsea) > FREQ_MIN and np.mean(fS) > FREQ_MIN:
                ratio = float(np.mean(fS) / np.mean(fsea))
                rr = max(ratio, 1 / ratio)
                locked = bool(round(rr) >= 2 and abs(rr - round(rr)) < TOL_LOCK)
            else:
                ratio = float("nan")           # 残余なし（塊=全体）: 比は定義不能と記録
                locked = False
            s_circ = float(np.sum(dev[S] * np.sin(2 * k * (th[S] - phi))))
            clumps.append({
                "k": k, "wavelength": f"1/{k}", "n_members": int(S.size),
                "energy": energy, "alignment": float(Ak[k]),
                "parity_winding": "odd" if k % 2 else "even",
                "parity_behavioral": "odd" if ovS.real < 0 else "even",
                "frame_deg": float(np.degrees(phi)),
                "freq_ratio_vs_rest_raw": ratio,
                "freq_spread_rel": float(np.std(fS) / max(np.mean(fS), 1e-30)),
                "nontrivial_lock": locked,
                "mass_proxy_amp_split_rms": float(np.std(dev[S])),
                "circulation_s": s_circ,
                "circulation_sign": int(np.sign(s_circ)) if s_circ != 0 else 0,
            })
    state["uniform_sea"] = bool(uniform)
    state["clumps"] = clumps

    # 対候補（同 k・逆周回符号）
    pairs = []
    for i in range(len(clumps)):
        for j in range(i + 1, len(clumps)):
            if (clumps[i]["k"] == clumps[j]["k"]
                    and clumps[i]["circulation_sign"] * clumps[j]["circulation_sign"] < 0):
                pairs.append((i, j))
    state["antiparticle_pair_candidates"] = pairs
    return state


def print_state(st):
    kind = "一様海（塊なし）" if st["uniform_sea"] else f"塊 {len(st['clumps'])} 群"
    print(f"  [{st['label']}] 段n={st['level']} 全体パリティ={st['parity_behavioral']} "
          f"(Re⟨⟩={st['parity_overlap_re']:+.4f}) 振幅割れRMS={st['amp_split_rms_total']:.2e} → {kind}")
    for c in st["clumps"]:
        print(f"    k={c['k']}（波長{c['wavelength']}）: 成分{c['n_members']} E={c['energy']:.3f} "
              f"整列度={c['alignment']:.3f} 巻きパリティ={c['parity_winding']}/行動={c['parity_behavioral']} "
              f"フレーム={c['frame_deg']:+.1f}° 周波数比(生)={c['freq_ratio_vs_rest_raw']:.6f} "
              f"ロック={c['nontrivial_lock']} 質量指標={c['mass_proxy_amp_split_rms']:.2e} "
              f"周回符号={c['circulation_sign']:+d}(s={c['circulation_s']:+.2e})")
    if st["antiparticle_pair_candidates"]:
        print(f"    対候補: {st['antiparticle_pair_candidates']}")


def main() -> None:
    t0 = time.time()
    print("粒子的グループ census v1 実行")
    print(f"  import: ABL {sha256(ABL)[:16]}…  MPH {sha256(MPH)[:16]}…")
    results = {"imports": {"abl": sha256(ABL), "mph": sha256(MPH),
                            "engine": mph.ENGINE_SHA256},
               "params": {"T_RELAX": T_RELAX, "WIN": WIN, "T_HALF": T_HALF,
                           "K_MAX": K_MAX, "A_PEAK": A_PEAK, "TOL_LOCK": TOL_LOCK}}

    for n, cfg in CFG.items():
        print(f"\n===== N={n}（M={n*(n-1)//2}） =====")
        states = []
        _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(n, False)
        st = census_state(n, 1, Z0, wp0.copy(), "control")
        st["family"] = "control"; print_state(st); states.append(st)
        Zh, info = mph.make_parent_harmonic(n, cfg["H"], cfg["seed"],
                                             iters=2000, restarts=10, tol=1e-12)
        for h in range(1, cfg["H"] + 1):
            lv = info["levels"][h - 1]
            fam = "N-1" if abs(lv["sigma1"] - (n - 1)) < 1e-9 else "broken"
            v0 = Zh[:, h - 1] * np.sqrt(cfg["H"])
            wp = np.random.default_rng(90000 + (h - 1)).normal(size=len(v0))
            st = census_state(n, h, v0, wp, f"harmonic n={h}")
            st["family"] = fam; st["sigma1"] = lv["sigma1"]
            print_state(st); states.append(st)
        results[f"N{n}"] = states

        nb = sum(1 for s in states if s["parity_behavioral"] == "even")
        nf = sum(1 for s in states if s["parity_behavioral"] == "odd")
        nc = sum(len(s["clumps"]) for s in states)
        npair = sum(len(s["antiparticle_pair_candidates"]) for s in states)
        nlock = sum(1 for s in states for c in s["clumps"] if c["nontrivial_lock"])
        print(f"  -- N={n} 集計: ボゾン型状態={nb} フェルミオン型状態={nf} "
              f"塊総数={nc} 対候補={npair} 非自明ロック={nlock}")

    results["runtime_sec"] = time.time() - t0
    (HERE / "paper8_particle_census_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
