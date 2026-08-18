#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
対照実験 CTRL-1: 第8論文が実際に使ったプログラムで固定点を取る。

  python3 run_control_paper8_a2a_v1.py

【CTRL-0 の誤り】
第8論文の λ=0.04937 を、第7論文の `onset_probe`（delta=1e-6、f=1-h1/htot）が
出したものと推定したが、**誤り**。論文が使ったのは

  次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/
      paper8_stage_A2a_seedless_N5/run_seedless.py

であり、条件は config_locked.json に凍結されている。

  親    : make_parent(sys, default_rng(40265722), iters=1200, tol=1e-12)
  初期  : Z0 = v.copy()（明示 seed なし・再正規化なし。delta 使用は forbidden）
  step  : 0..5000
  f     : transverse = Z - p(p·Z) - q(q·Z);  f = |transverse|² / |Z|²

onset_probe の `f = 1 - h1/htot` とは**別の式**である。1 から引く形は桁落ちで
1e-12 付近が床になるが、射影差分の形は桁落ちしないため 3.27e-33 まで下りる。
論文が数値床 1e-30 台から 18桁を測れたのはこの違いによる。

【本対照の判定】
論文の実出力 raw/A2a_N5_seedless_f64_e1/f_timeseries.csv（5001 step）の
f 列を、本環境で再計算した f と **.17e 文字列で全行一致**するか比較する。
許容誤差は置かない。全行一致して初めて σ₁ 計装の資格ができる。

論文成果物には一切書き込まない（読むだけ）。
"""

import csv
import hashlib
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
PKG = os.path.join(REPO, "次元の生成構造",
                   "第8論文_二段階seed除去による準安定相の因果分離",
                   "paper8_stage_A2a_seedless_N5")
# 論文の expected_hashes.json が指定する原本の所在（第7論文の作業コピーではない）
ENGINE = os.path.join(REPO, "時間軸Q軸とフェルミオンの生成構造", "検証_対照実験",
                      "第5論文原本_自発的分裂予備実験_v1",
                      "run_n_scaling_lowrank_v1.py")
FLOAT_FMT = ".17e"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    cfg = json.load(open(os.path.join(PKG, "config_locked.json"), encoding="utf-8"))
    exp = json.load(open(os.path.join(PKG, "expected_hashes.json"), encoding="utf-8"))

    want = exp["sources"]["run_n_scaling_lowrank_v1.py"]["sha256"]
    got = sha256(ENGINE)
    print(f"原本(論文指定パス): {ENGINE}")
    print(f"SHA-256 {'一致' if got == want else '★不一致'}  {got[:16]}…")
    if got != want:
        raise SystemExit("SOURCE_MISMATCH")

    spec = importlib.util.spec_from_file_location("lowrank_paper8", ENGINE)
    lowrank = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lowrank)

    n = cfg["n"]
    print(f"凍結条件: N={n} seed={cfg['parent_prng_seed']} "
          f"iters={cfg['parent_iters']} tol={cfg['parent_tolerance']} "
          f"max_step={cfg['max_step']} Z0={cfg['initial_state']}")

    rng = np.random.default_rng(cfg["parent_prng_seed"])
    sys_lr = lowrank.LowRankSystem(n)
    v, parent_residual, parent_sigma = lowrank.make_parent(
        sys_lr, rng, iters=cfg["parent_iters"], tol=cfg["parent_tolerance"])
    print(f"親残差 = {parent_residual:.6e}   親σスペクトル本数 = {len(parent_sigma)}")

    Z = v.copy()
    assert np.array_equal(Z, v)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=sys_lr.m)     # 状態には加えない

    def f_value(state):
        tr = state - p * (p @ state) - q * (q @ state)
        return float(np.real(np.vdot(tr, tr))) / float(np.real(np.vdot(state, state)))

    # 計装はここだけ: σ₁ を毎 step 保存する（力学には一切戻さない）
    f_mine, sigma_hist = [], []
    for step in range(cfg["max_step"] + 1):
        f_mine.append(f_value(Z))
        if step == cfg["max_step"]:
            break
        sys_lr.set_theta(np.angle(Z))
        sigma_estimate, wp = sys_lr.sigma_max_power(wp)
        sigma_hist.append(float(sigma_estimate))
        Z = sys_lr.cayley_step(Z, sigma_estimate)

    ref_path = os.path.join(PKG, "raw", "A2a_N5_seedless_f64_e1", "f_timeseries.csv")
    ref = [row["f"] for row in csv.DictReader(open(ref_path, encoding="utf-8"))]
    print(f"\n論文出力: {os.path.relpath(ref_path, REPO)}  ({len(ref)} 行)")

    mism = [(i, ref[i], format(f_mine[i], FLOAT_FMT))
            for i in range(min(len(ref), len(f_mine)))
            if ref[i] != format(f_mine[i], FLOAT_FMT)]
    exact = (len(ref) == len(f_mine)) and not mism

    print(f"step 0   論文 {ref[0]}   本環境 {format(f_mine[0], FLOAT_FMT)}")
    print(f"step 100 論文 {ref[100]}   本環境 {format(f_mine[100], FLOAT_FMT)}")
    print(f"step 5000論文 {ref[5000]}   本環境 {format(f_mine[5000], FLOAT_FMT)}")
    print(f"\n不一致 step 数: {len(mism)} / {len(ref)}")
    for i, a, b in mism[:5]:
        print(f"  step {i}: 論文 {a}  本環境 {b}")

    s = np.array(sigma_hist)
    print(f"\n[計装] σ₁: 初期 {s[0]:.10e}  最終 {s[-1]:.10e}  "
          f"最小 {s.min():.6e}  最大 {s.max():.6e}  最大/最小 {s.max()/s.min():.6f}")

    out = {
        "control": "CTRL-1 第8論文 Stage A2a 再現",
        "engine": ENGINE, "engine_sha256": got,
        "config": {k: cfg[k] for k in
                   ("n", "dtype", "parent_prng_seed", "parent_iters",
                    "parent_tolerance", "initial_state", "max_step")},
        "parent_residual": parent_residual,
        "reference_csv": os.path.relpath(ref_path, REPO),
        "rows_reference": len(ref), "rows_reproduced": len(f_mine),
        "mismatch_count": len(mism),
        "first_mismatches": [{"step": i, "paper": a, "here": b} for i, a, b in mism[:20]],
        "f_timeseries_exact": bool(exact),
        "sigma1": {"n": len(sigma_hist), "first": s[0], "last": s[-1],
                   "min": float(s.min()), "max": float(s.max()),
                   "max_over_min": float(s.max() / s.min())},
    }
    op = os.path.join(HERE, "result_control_paper8_a2a_v1.json")
    json.dump(out, open(op, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    np.save(os.path.join(HERE, "sigma1_timeseries_A2a_N5_v1.npy"), s)
    print(f"\n出力: {op}")
    print("\n" + ("=== PASS。論文出力と全 5001 step が .17e で厳密一致。固定点を確保した。 ==="
                  if exact else
                  "=== FAIL。論文出力と一致しない。環境差を先に潰すこと。 ==="))
    sys.exit(0 if exact else 1)


if __name__ == "__main__":
    main()
