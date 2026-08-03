#!/usr/bin/env python3
"""第8論文v2予備実験 E-M1：物質はどこで生まれるか——閉鎖完走のタイミング測定 v1

背景:
    第9論文本体は「粒子=閉鎖周期を完走した稀な事象」とし、完走には安定な
    レジスタ（凍結した周波数構造）が要ると論じた。ならば物質の誕生は、
    インフレーション的拡大の最中ではなく、拡大停止（準安定化）の瞬間で
    なければならない。本実験は第7/8論文の無seed系（条件A・N=5）で、
    この予言を実測する。

操作的定義（測定前固定）:
    - 各関係辺 e の瞬間周波数 ω_e(t) = 位相の窓内傾き（窓W=200, 刻み50）
    - 光の海（自明）: 全辺が単一周波数のユニゾン（比1）
    - 物質候補（非自明ロック）: 辺対 (e,f) が |ω̄_e/ω̄_f| ≈ 整数 p≥2 に
      tol=1e-3 でロックし、連続5窓以上持続するもの

予言（測定前固定）:
    P1（海の斉唱）: 潜伏期（f<0.05）は全辺ほぼ単一周波数（分散小）で、
        非自明ロック数 L(t) ≈ 0
    P2（拡大中は生まれない）: crossing〜準安定開始（crossing+3000）の
        急拡大期は周波数が滑り、持続的な非自明ロックは形成されない
    P3（物質誕生=停止）: 準安定開始後に周波数クラスタが凍結し、
        非自明ロック L(t) が立ち上がり以後持続する（化石化）
    反証条件: 急拡大中に持続的非自明ロックが出る、または準安定後も
        L=0 のままなら、P2/P3 は反証として記録する。

実測結果（v1・記録）:
    P1成立・P2成立・【P3反証】。診断: 準安定期の系は全辺・全方向射影が
    単一周波数 ω≈0.04363 のユニゾン（周波数比のずれ最大 5e-4、
    非自明整数比なし）。N=5 の無seed関係系では、停止後も系は
    「三方向＋単一周波数の光の海」のままであり、物質（調和族＝
    整数比周波数ラダー）は生まれない。
    帰結: 拡大停止は物質誕生の必要条件（P1/P2: 拡大中は生まれない）
    だが十分条件ではない。物質にはレジスタ（周波数の整数ラダー＝
    スケール層）の分化が別途必要であり、これは幾何辞書のスケール層
    （距離・質量・固有時の三辞書）と同一の未解決問題である。
    N=5 は共形層（方向）だけを生み、スケール層を生まなかった——
    層分離の実測的裏付け。次段: N=40/300 での周波数分化の有無。

規約: 第7/8論文コードは read-only import（無変更）。条件Aは乱数を
      seed 生成に消費しない（第8論文と同一の初期化）。
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CODE = HERE.parent / "code" / "run_preliminary_seed_ablation_v1.py"
spec = importlib.util.spec_from_file_location("ablation_for_matter_birth_v1", CODE)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)

N = 5
XMAX = 20000            # crossing+GUARD を十分に超える長さ（N=5 の既知時間尺度）
WIN = 200
STRIDE = 50
RATIO_TOL = 1e-3
PERSIST = 5
FREQ_MIN = 1e-6


def main() -> None:
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = abl.build_init(N, initial_seed=False)
    M = sys_lr.m

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    phases = np.zeros((XMAX + 1, M))
    fs = np.zeros(XMAX + 1)
    phases[0] = np.angle(Z); fs[0] = fval(Z)
    crossing = None
    for t in range(1, XMAX + 1):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z)
        fs[t] = fval(Z)
        if crossing is None and fs[t] > 0.05:
            crossing = t
    meta_start = crossing + abl.GUARD if crossing is not None else None
    print(f"crossing(f>0.05) = {crossing}, 準安定開始 = {meta_start}")

    unwrapped = np.unwrap(phases, axis=0)
    centers, freqs = [], []
    for s in range(0, XMAX - WIN, STRIDE):
        seg = unwrapped[s:s + WIN]
        x = np.arange(WIN)
        slope = np.polyfit(x, seg, 1)[0]        # 辺ごとの位相傾き
        centers.append(s + WIN // 2)
        freqs.append(slope)
    centers = np.array(centers); freqs = np.abs(np.array(freqs))

    # 非自明整数比ロックの数え上げ（持続性つき）
    lock_now = []
    for w in range(len(centers)):
        f = freqs[w]
        locked = set()
        for i in range(M):
            for j in range(i + 1, M):
                hi, lo = max(f[i], f[j]), min(f[i], f[j])
                if lo < FREQ_MIN:
                    continue
                r = hi / lo
                pr = round(r)
                if pr >= 2 and abs(r - pr) < RATIO_TOL:
                    locked.add((i, j, pr))
        lock_now.append(locked)
    L = np.zeros(len(centers), dtype=int)
    for w in range(len(centers)):
        cnt = 0
        for key in lock_now[w]:
            run = 1
            k = w - 1
            while k >= 0 and key in lock_now[k]:
                run += 1; k -= 1
            k = w + 1
            while k < len(centers) and key in lock_now[k]:
                run += 1; k += 1
            if run >= PERSIST:
                cnt += 1
        L[w] = cnt

    disp = np.std(freqs, axis=1) / np.maximum(np.mean(freqs, axis=1), 1e-30)

    def phase_of(step):
        if crossing is None or step < crossing:
            return "latent"
        if step < meta_start:
            return "expansion"
        return "metastable"

    seg_L = {"latent": [], "expansion": [], "metastable": []}
    for w, c in enumerate(centers):
        seg_L[phase_of(c)].append(int(L[w]))
    stats = {k: {"max": (max(v) if v else 0), "mean": (float(np.mean(v)) if v else 0.0),
                 "windows": len(v)} for k, v in seg_L.items()}
    for k, s in stats.items():
        print(f"{k:10s}: 窓数={s['windows']:4d} 非自明ロック 最大={s['max']} 平均={s['mean']:.2f}")

    p1 = stats["latent"]["max"] == 0
    p2 = stats["expansion"]["max"] == 0
    meta_L = seg_L["metastable"]
    p3 = len(meta_L) > 0 and max(meta_L) >= 1 and (np.mean(meta_L[len(meta_L)//2:]) >= 1)

    # 診断: 準安定末期の周波数構造（ユニゾン判定）
    tail = freqs[-60:]
    fbar = np.mean(tail, axis=0)
    rmat = fbar[:, None] / np.maximum(fbar[None, :], 1e-30)
    max_ratio_dev = float(np.max(np.abs(rmat[rmat >= 1] - 1)))
    print(f"診断: 準安定末期の辺周波数 平均={np.mean(fbar):.6f} "
          f"比のずれ最大={max_ratio_dev:.2e} → {'ユニゾン（単一周波数の海）' if max_ratio_dev < 1e-2 else '分化あり'}")
    print(f"\nP1 海の斉唱（潜伏期 L=0）: {'PASS' if p1 else 'FAIL'}")
    print(f"P2 拡大中は生まれない（急拡大期 L=0）: {'PASS' if p2 else 'FAIL'}")
    print(f"P3 物質誕生=停止（準安定期に L≥1 が立ち持続）: {'PASS' if p3 else 'FAIL'}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True, constrained_layout=True)
    axes[0].plot(fs, lw=0.8); axes[0].set_yscale("log")
    axes[0].set_ylabel("f outside parent"); axes[0].set_title(f"N={N} seedless (condition A)")
    for ax in axes:
        if crossing: ax.axvline(crossing, color="tab:orange", ls=":", label="crossing")
        if meta_start: ax.axvline(meta_start, color="tab:red", ls=":", label="metastable start")
    for i in range(M):
        axes[1].plot(centers, freqs[:, i], lw=0.7)
    axes[1].set_ylabel("edge |frequency|")
    axes[2].plot(centers, L, lw=1.2, color="tab:green")
    axes[2].set_ylabel("nontrivial integer locks L(t)"); axes[2].set_xlabel("step")
    axes[0].legend(loc="lower right")
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"paper8_matter_birth_v1.{ext}", dpi=160)
    plt.close(fig)

    np.savetxt(HERE / "paper8_matter_birth_L_v1.csv",
               np.column_stack([centers, L, disp]), delimiter=",",
               header="window_center,nontrivial_locks,freq_dispersion", comments="")
    payload = {
        "experiment": "paper8_matter_birth_timing_pre_v1",
        "engine": "第7/8論文エンジン read-only import（条件A・無seed・N=5）",
        "params": {"N": N, "XMAX": XMAX, "WIN": WIN, "STRIDE": STRIDE,
                   "RATIO_TOL": RATIO_TOL, "PERSIST": PERSIST},
        "crossing": crossing, "metastable_start": meta_start,
        "phase_stats": stats,
        "P1_latent_silence": bool(p1), "P2_no_birth_during_expansion": bool(p2),
        "P3_birth_at_arrest": bool(p3),
        "diagnosis": {"metastable_mean_freq": float(np.mean(fbar)),
                       "max_ratio_deviation_from_unison": max_ratio_dev},
        "conclusion": (
            f"P1={'成立' if p1 else '反証'}, P2={'成立' if p2 else '反証'}, "
            f"P3={'成立' if p3 else '反証'}。P3反証の診断: 準安定期は単一周波数"
            f"ω≈{np.mean(fbar):.5f} のユニゾン（比ずれ最大 {max_ratio_dev:.1e}）。"
            "拡大停止は物質誕生の必要条件だが十分条件ではない——N=5では停止後も"
            "三方向＋光の海のままで、物質（整数比周波数ラダー＝レジスタ）は"
            "生まれない。スケール層の分化が別途必要（共形層/スケール層分離の実測）。"
            "次段: N=40/300での周波数分化の有無"),
    }
    (HERE / "paper8_matter_birth_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print("\nsaved: paper8_matter_birth_result_v1.json")


if __name__ == "__main__":
    main()
