#!/usr/bin/env python3
"""E-M9r 補正 v2：閉塞帳簿の訂正（係数閉塞→場の閉塞）

v1 の誤り（2026-08-04 木原氏指摘で判明・訂正対象は帳簿のみ）:
    v1 はセクター場 w_k[m,n] = c[m,k]·e^{2πikn/N} を係数ベクトル C[:,k]∈C^M に
    潰し、係数の二次形式 |Σ_m c²| を「初期閉塞」として報告した。これは誤った
    二次形式である。物理的な閉塞は (m,n) 場の全成分二乗和
        Σ_{m,n} w_k² = (Σ_m c[m,k]²)·N·δ_{2k≡0 mod N}
    であり、2k≢0 (mod N) のセクターは恒等的に零閉塞する。

軌道は訂正不要:
    セクター場の任意の標本スライスは同一の C^M 射線（n 依存は大域位相=ゲージ、
    K は位相差のみに依存）なので、v1 の注入軌道はセクター力学として正しい。
    さらに Cayley は Σ_m Z² を厳密保存するため、場の閉塞
    (Σ_m Z(t)²)·N·δ_{2k≡0} は発展を通じて不変——2k≢0 セクターは全時刻で厳密閉塞。

本スクリプトの処理:
    (1) 各セクターの正しい場の閉塞を計算し、v1 の力学結果と結合した訂正表を出力
    (2) 保存検証: 代表セクターを短時間発展させ、Σ_m Z² の保存（場の閉塞の不変性）を実測
    (3) 自己対セクター（k=0、偶数Nの k=N/2）は生成器欠陥 (b) の対象として区分表示

使い方: python3 run_paper8_em9r_field_closure_fix_v2.py <N>
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
EM9S = (PAPER8 / "paper8_em9_harmonic_initial_inflation_pre_v1"
        / "run_paper8_em9_N40_N300_supplement_v1.py")
CENSUS = REPO / "次元の生成構造" / "standalone_parent_census_v1"
PARENT_DIR = {5: CENSUS / "parent_white_harmonics_N5_v2",
              40: CENSUS / "parent_white_harmonics_N40_v2"}

spec = importlib.util.spec_from_file_location("em9s_fx", EM9S)
em9s = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = em9s
spec.loader.exec_module(em9s)
abl = em9s.abl

VERIFY_STEPS = 500


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    t0 = time.time()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    pdir = PARENT_DIR[n]
    waves = np.load(pdir / "relation_waves.npy")
    m, nn = waves.shape
    v1 = json.loads((HERE / f"paper8_em9r_result_N{n:05d}_v1.json").read_text(encoding="utf-8"))
    C = np.fft.fft(waves, axis=1) / n
    print(f"E-M9r 閉塞帳簿訂正 v2  N={n}（M={m}）")
    print(f"  正本 sha256: {sha256(pdir / 'relation_waves.npy')[:16]}…  v1結果を帳簿訂正")

    rows = []
    for k in range(n):
        Zk = C[:, k]
        nk2 = float(np.sum(np.abs(Zk) ** 2))
        coef2 = complex(np.sum(Zk ** 2))
        selfpair = (2 * k) % n == 0
        nsum = n if selfpair else 0.0
        field_closure = abs(coef2) * nsum          # 場の閉塞（絶対値）
        field_closure_rel = (abs(coef2) / nk2) if selfpair else 0.0
        d = v1["components"][f"k{k}"]
        rows.append({
            "k": k, "self_paired": selfpair,
            "field_closure_rel": field_closure_rel,
            "coef_closure_rel_v1_ARTIFACT": d["initial_closure_rel"],
            "energy_share": d["energy_share"],
            "crossing": d["crossing"],
            "rank_Q_metastable_mode": d["rank_Q_metastable_mode"],
            "f_late_mean": d["f_late_mean"],
            "axiom_status": ("非閉塞（生成器欠陥(b)対象）" if selfpair
                              else "厳密零閉塞（正当な波）"),
        })

    # (2) 保存検証: 2k≢0 の代表セクターと自己対セクターを短時間発展
    verify = {}
    for k in ([1, 0] if n % 2 else [1, 0, n // 2]):
        Zk = C[:, k] / np.linalg.norm(C[:, k])
        sys_lr = abl.LowRankSystem(n)
        sys_lr.set_theta(np.angle(Zk))
        wp = np.random.default_rng(93000 + k).normal(size=m)
        Z = Zk.copy()
        c0 = complex(Z @ Z)
        for _ in range(VERIFY_STEPS):
            Z, wp = abl.evolve(sys_lr, Z, wp)
        c1 = complex(Z @ Z)
        drift = abs(c1 - c0)
        selfpair = (2 * k) % n == 0
        f0 = abs(c0) * (n if selfpair else 0.0)
        f1 = abs(c1) * (n if selfpair else 0.0)
        verify[f"k{k}"] = {"coef_bilinear_drift": drift,
                            "field_closure_t0": f0, "field_closure_t500": f1}
        print(f"  保存検証 k={k}: |Σ Z²| ドリフト={drift:.2e} "
              f"場の閉塞 t=0→{VERIFY_STEPS}: {f0:.3e} → {f1:.3e}")

    # 訂正表
    print(f"\n  k | 自己対 | 場の閉塞(相対) | v1誤帳簿(係数閉塞) | E比 | crossing | rank_Q | 公理判定")
    show = rows if n <= 6 else rows[:4] + [r for r in rows if r["self_paired"]]
    for r in show:
        print(f"  {r['k']:>2} | {'是' if r['self_paired'] else '否'} | "
              f"{r['field_closure_rel']:.3e} | {r['coef_closure_rel_v1_ARTIFACT']:.3f} | "
              f"{r['energy_share']:.4f} | {r['crossing']} | {r['rank_Q_metastable_mode']} | "
              f"{r['axiom_status']}")

    legit = [r for r in rows if not r["self_paired"]]
    defect = [r for r in rows if r["self_paired"]]
    print(f"\n==== 訂正後まとめ N={n} ====")
    print(f"正当な閉塞波（2k≢0）: {len(legit)} 本——全て厳密零閉塞・"
          f"crossing {min(r['crossing'] for r in legit)}〜{max(r['crossing'] for r in legit)}・"
          f"rank_Q=4 が {sum(1 for r in legit if r['rank_Q_metastable_mode'] == 4)} 本")
    dk = [r["k"] for r in defect]
    dv = [f"{r['field_closure_rel']:.3e}" for r in defect]
    print(f"自己対セクター（非閉塞・生成器欠陥(b)）: {len(defect)} 本 (k={dk}, 場の閉塞相対 {dv})")

    out = {"N": n, "correction": "係数閉塞（誤）→場の閉塞（正）。軌道はセクター力学として有効",
           "rows": rows, "conservation_verify": verify,
           "v1_artifact_note": ("v1のinitial_closure_relは係数二次形式であり物理閉塞ではない。"
                                 "2k≢0セクターは場として恒等零閉塞・発展で厳密保存"),
           "runtime_sec": time.time() - t0}
    (HERE / f"paper8_em9r_fix_N{n:05d}_v2.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
