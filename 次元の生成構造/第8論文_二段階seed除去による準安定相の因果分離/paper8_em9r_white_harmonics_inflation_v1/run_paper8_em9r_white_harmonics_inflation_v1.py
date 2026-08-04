#!/usr/bin/env python3
"""E-M9r：新 make_parent（白色・零閉塞・倍音事後読出し版）によるインフレーション実験のやり直し v1

背景（2026-08-04 木原氏指示）:
    倍音対応 make_parent（H段・事前配置）は廃止され、Codex セッションで
    make_parent_white_harmonics_n_only_v2 に置き換えられた:
      - 理論入力は N のみ（H 廃止）。M=N(N-1)/2 本の関係波が各 N 点の内部波形を持つ
      - 各関係波は等ノルム直交二成分 (q1+iq2)/√2 により恒等的に零閉塞
      - 倍音は事前配置せず、N点DFTで事後読出し
      - シードは整数を順に試し最初に収束したものを採用（全試行監査）
    本実験は E-M9（倍音海注入インフレーション再現）を新生成器でやり直す。

初期値（正本を read-only 読込・再生成しない）:
    standalone_parent_census_v1/parent_white_harmonics_N5_v2（seed 2）
    standalone_parent_census_v1/parent_white_harmonics_N40_v2（seed 1）
    ——安定波分類実験（52c83f13）と同一の状態。

注入系列（各 N）:
    (1) parent_vector v ∈ C^M（古典的親＝振幅位相階層）
    (2) 倍音成分 k=0..N−1: relation_waves の標本軸 DFT
        C[:,k] = (1/N)Σ_n relation_waves[:,n]·e^{−2πikn/N} ∈ C^M を正規化して注入。
        成分は個別には零閉塞しない——初期 |Z·Z| を観測量として記録
        （Cayley発展は Z·Z を厳密保存するため、この値は不変量として運ばれる）。
    時計は全系列 micro=1（標本軸の倍音番号 k は読出しラベル。
    段別時計 n·ω₀ の力学実装は E-M11 で確立済みだが、本実験では導入しない）。

測定: E-M9 と同一（abl 無改変 read-only、crossing=f>0.05・準安定 rank_Q 最頻値、
    XMAX=12000、run_injected は E-M9 補足スクリプトから import）。
    P0: 対照（既存 build_init）の crossing 再現（N=5:1166 / N=40:2011）。

使い方: python3 run_paper8_em9r_white_harmonics_inflation_v1.py <N>
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
GENERATOR = CENSUS / "make_parent_white_harmonics_n_only_v2.py"
EXPECTED_CONTROL_CROSSING = {5: 1166, 40: 2011}

spec = importlib.util.spec_from_file_location("em9s_r", EM9S)
em9s = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = em9s
spec.loader.exec_module(em9s)
abl = em9s.abl


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    t0 = time.time()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    pdir = PARENT_DIR[n]
    manifest = json.loads((pdir / "manifest.json").read_text(encoding="utf-8"))
    waves = np.load(pdir / "relation_waves.npy")
    parent = np.load(pdir / "parent_vector.npy")
    m, nn = waves.shape
    assert nn == n and m == n * (n - 1) // 2

    print(f"E-M9r 白色倍音 make_parent 版インフレーションやり直し N={n}")
    print(f"  正本: {pdir.name} seed={manifest['accepted_seed']} shape={waves.shape}")
    print(f"  sha256: relation_waves={sha256(pdir / 'relation_waves.npy')[:16]}… "
          f"generator={sha256(GENERATOR)[:16]}…")

    results = {"N": n, "source": pdir.name,
               "accepted_seed": manifest["accepted_seed"],
               "imports": {"relation_waves": sha256(pdir / "relation_waves.npy"),
                            "parent_vector": sha256(pdir / "parent_vector.npy"),
                            "generator": sha256(GENERATOR),
                            "em9_supplement": sha256(EM9S)},
               "params": {"XMAX": em9s.XMAX, "clock": "micro=1（kは読出しラベル）",
                           "dft": "C[:,k]=(1/N)Σ_n w[:,n]e^{-2πikn/N}"}}

    # ---- P0 対照 ----
    print(f"\n[対照] abl.build_init({n}, initial_seed=False)")
    _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(n, False)
    ctrl = em9s.run_injected(n, Z0, wp0.copy(), "control")
    p0 = ctrl["crossing"] == EXPECTED_CONTROL_CROSSING[n]
    print(f"  P0 駆動検証（crossing={EXPECTED_CONTROL_CROSSING[n]} 再現）: "
          f"{'PASS' if p0 else 'FAIL'}")
    results["control"] = {k: ctrl[k] for k in
                           ("crossing", "metastable_start", "rank_Q_metastable_mode",
                            "f_late_mean")}
    results["P0"] = bool(p0)

    # ---- (1) 古典的親 ----
    print(f"\n[白色版 親ベクトル v]（‖v‖={np.linalg.norm(parent):.6f} "
          f"|v·v|={abs(complex(parent @ parent)):.1e}）")
    wp = np.random.default_rng(91000).normal(size=m)
    r = em9s.run_injected(n, parent / np.linalg.norm(parent), wp, "parent_vector")
    results["parent_vector"] = {**{k: r[k] for k in
                                    ("crossing", "metastable_start",
                                     "rank_Q_metastable_mode", "f_late_mean")},
                                 "initial_closure_abs": abs(complex(parent @ parent))}

    # ---- (2) DFT 倍音成分 k=0..N−1 ----
    C = np.fft.fft(waves, axis=1) / n            # C[:,k] ∈ C^M
    energy = np.sum(np.abs(C) ** 2, axis=0)
    energy_share = energy / energy.sum()
    print(f"\n[倍音成分] エネルギー配分: "
          f"{np.array2string(energy_share, precision=4, suppress_small=True)}")
    comps = {}
    for k in range(n):
        Zk = C[:, k]
        nk = float(np.linalg.norm(Zk))
        closure0 = abs(complex(Zk @ Zk)) / max(nk ** 2, 1e-300)
        wp = np.random.default_rng(92000 + k).normal(size=m)
        r = em9s.run_injected(n, Zk / nk, wp, f"k={k} (E比={energy_share[k]:.4f})")
        comps[f"k{k}"] = {**{kk: r[kk] for kk in
                              ("crossing", "metastable_start",
                               "rank_Q_metastable_mode", "f_late_mean")},
                           "energy_share": float(energy_share[k]),
                           "initial_closure_rel": float(closure0)}
    results["components"] = comps

    # ---- まとめ ----
    ncross = sum(1 for c in comps.values() if c["crossing"] is not None)
    r4 = sum(1 for c in comps.values() if c["rank_Q_metastable_mode"] == 4)
    pv = results["parent_vector"]
    print(f"\n==== まとめ N={n} ====")
    print(f"親ベクトル: crossing={pv['crossing']} rank_Q={pv['rank_Q_metastable_mode']}")
    print(f"倍音成分 {n} 本中: crossing 発生={ncross} 本, 準安定 rank_Q=4 が {r4} 本")
    results["summary"] = {"components_total": n, "components_crossed": ncross,
                           "components_rank4": r4}
    results["runtime_sec"] = time.time() - t0
    (HERE / f"paper8_em9r_result_N{n:05d}_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float),
        encoding="utf-8")
    print(f"saved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
