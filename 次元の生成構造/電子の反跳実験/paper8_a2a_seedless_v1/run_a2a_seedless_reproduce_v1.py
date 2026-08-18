#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
シードなし系列（第8論文 Stage A2a）の再現プログラム。

  python3 run_a2a_seedless_reproduce_v1.py

【この系列の位置づけ】
シードあり系列（第7論文 δ=1e-15、`paper7_f_projection_v1/`）とは**別系列**である。
親状態 v は両者で同一（PRNG seed 40265722、iters=1200、tol=1e-12）だが、
初期状態が違う。

    第7論文系列 : Z0 = (v + 1e-15 * g) / |v + 1e-15 * g|    → f(0) = 1.066e-30
    本系列 A2a  : Z0 = v.copy()（明示シードなし）           → f(0) = 3.275e-33

混同を避けるため、データ・図・プログラムをフォルダごと分離する。

【原本】
  次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/
    paper8_stage_A2a_seedless_N5/run_seedless.py
  凍結条件は同フォルダ config_locked.json、原本SHAは expected_hashes.json。
  `delta` / `zero_closure_kernel_seed` の使用は config の forbidden に明記されている。

【方針】
  力学（LowRankSystem, make_parent, sigma_max_power, cayley_step）は原本を
  import してそのまま使う。転写するのは run_seedless.py の走行ループのうち
  f 系列を作る部分だけで、変更点は出力先を本フォルダに変えることのみ。
  原本 run_seedless.py は raw/ が空でないと実行を拒否する設計なので直接は呼べない。

【判定】
  公開済み
    paper8_stage_A2a_seedless_N5/raw/A2a_N5_seedless_f64_e1/f_timeseries.csv
  の f 列（5001 step）と、本環境の再計算値が .17e 文字列で全行一致すること。
  許容誤差は置かない。一致しなければ図を作ってはならない。
"""

import csv
import hashlib
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
PKG = os.path.join(REPO, "次元の生成構造",
                   "第8論文_二段階seed除去による準安定相の因果分離",
                   "paper8_stage_A2a_seedless_N5")
ENGINE = os.path.join(REPO, "時間軸Q軸とフェルミオンの生成構造", "検証_対照実験",
                      "第5論文原本_自発的分裂予備実験_v1",
                      "run_n_scaling_lowrank_v1.py")
FLOAT_FMT = ".17e"
RUN_ID = "A2a_N5_seedless_f64_e1"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    cfg = json.load(open(os.path.join(PKG, "config_locked.json"), encoding="utf-8"))
    exp = json.load(open(os.path.join(PKG, "expected_hashes.json"), encoding="utf-8"))

    print("=== 原本の照合 ===")
    want = exp["sources"]["run_n_scaling_lowrank_v1.py"]["sha256"]
    got = sha256(ENGINE)
    print(f"  {'OK ' if got == want else '★NG'} run_n_scaling_lowrank_v1.py  {got[:16]}…")
    if got != want:
        raise SystemExit("SOURCE_MISMATCH")

    print("\n=== 凍結条件（config_locked.json） ===")
    for k in ("n", "dtype", "parent_prng_seed", "parent_iters", "parent_tolerance",
              "initial_state", "max_step"):
        print(f"  {k:<18} {cfg[k]}")
    print(f"  forbidden          delta={cfg['forbidden']['delta']}, "
          f"zero_closure_kernel_seed={cfg['forbidden']['zero_closure_kernel_seed']}, "
          f"high_precision={cfg['forbidden']['high_precision']}")

    spec = importlib.util.spec_from_file_location("lowrank_a2a", ENGINE)
    lowrank = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lowrank)

    n = cfg["n"]
    t0 = time.time()
    rng = np.random.default_rng(cfg["parent_prng_seed"])
    sys_lr = lowrank.LowRankSystem(n)
    v, parent_residual, parent_sigma = lowrank.make_parent(
        sys_lr, rng, iters=cfg["parent_iters"], tol=cfg["parent_tolerance"])
    print(f"\n親残差 = {parent_residual:.6e}   親σスペクトル本数 = {len(parent_sigma)}")

    Z = v.copy()                       # 凍結条件: Z0 = v.copy()（シードなし）
    if not np.array_equal(Z, v):
        raise SystemExit("Z0 is not bitwise equal to v")
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=sys_lr.m)     # 状態には加えない（原本の規約）

    def f_value(state):
        tr = state - p * (p @ state) - q * (q @ state)
        return float(np.real(np.vdot(tr, tr))) / float(np.real(np.vdot(state, state)))

    initial_norm = float(np.real(np.vdot(Z, Z)))
    rows = []
    for step in range(cfg["max_step"] + 1):
        norm_sq = float(np.real(np.vdot(Z, Z)))
        f = f_value(Z)
        rows.append([str(step),
                     format(f, FLOAT_FMT),
                     format(math.log10(f) if f > 0.0 else float("nan"), FLOAT_FMT),
                     format(abs(norm_sq - 1.0), FLOAT_FMT),
                     format(abs(complex(Z @ Z)), FLOAT_FMT),
                     format(abs(norm_sq - initial_norm), FLOAT_FMT)])
        if step == cfg["max_step"]:
            break
        sys_lr.set_theta(np.angle(Z))
        sigma_estimate, wp = sys_lr.sigma_max_power(wp)
        if not np.isfinite(sigma_estimate) or sigma_estimate <= 0.0:
            raise SystemExit(f"invalid sigma at step {step}")
        Z = sys_lr.cayley_step(Z, sigma_estimate)
    elapsed = time.time() - t0

    # --- 判定: 公開CSVの f 列と全行一致 ---
    ref_path = os.path.join(PKG, "raw", RUN_ID, "f_timeseries.csv")
    ref = list(csv.DictReader(open(ref_path, encoding="utf-8")))
    mism = [(i, ref[i]["f"], rows[i][1]) for i in range(min(len(ref), len(rows)))
            if ref[i]["f"] != rows[i][1]]
    exact = (len(ref) == len(rows)) and not mism
    print(f"\n公開CSV照合: {len(ref)} 行 vs {len(rows)} 行、不一致 {len(mism)} → "
          f"{'一致' if exact else '★不一致'}")
    for i, a, b in mism[:3]:
        print(f"  step {i}: 公開 {a}  本環境 {b}")
    if not exact:
        raise SystemExit("REPRODUCTION_FAILED")

    out = os.path.join(HERE, f"a2a_seedless_f_timeseries_N00005.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "f", "log10_f", "norm_error", "closure_error", "conservation_error"])
        w.writerows(rows)

    fv = np.array([float(r[1]) for r in rows])
    meta = {"series": "seedless (第8論文 Stage A2a)", "N": n,
            "initial_state": cfg["initial_state"],
            "parent_prng_seed": cfg["parent_prng_seed"],
            "parent_iters": cfg["parent_iters"], "parent_tolerance": cfg["parent_tolerance"],
            "max_step": cfg["max_step"], "rows": len(rows),
            "parent_residual": parent_residual,
            "reference_csv": os.path.relpath(ref_path, REPO),
            "reproduction_exact": bool(exact),
            "f_first": float(fv[0]), "f_min": float(fv.min()), "f_max": float(fv.max()),
            "elapsed_sec": elapsed}
    json.dump(meta, open(os.path.join(HERE, "a2a_seedless_meta_N00005.json"), "w",
                         encoding="utf-8"), indent=2, ensure_ascii=False)
    print(f"\nf: 初期 {fv[0]:.17e}  最小 {fv.min():.6e}  最大 {fv.max():.6e}  ({elapsed:.1f}s)")
    print(f"→ {os.path.basename(out)}")
    print("\n=== PASS。公開CSVと全 5001 step が一致。図化に進んでよい。 ===")


if __name__ == "__main__":
    main()
