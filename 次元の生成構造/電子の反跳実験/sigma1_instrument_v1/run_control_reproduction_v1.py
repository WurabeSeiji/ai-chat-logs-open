#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
対照実験 CTRL-0: 計装前に固定点を確保する。

  python3 run_control_reproduction_v1.py

原本エンジン `次元の生成構造/自発的分裂予備実験_v1/run_n_scaling_lowrank_v1.py`
（SHA-256 ba0fc19b…、第7論文リリースノートで「不変更」と宣言済み）を
read-only import し、公開済み出力 summary_N00005_seed0.json（2026-07-21 21:45）
を厳密再現できるかを確認する。

二つの経路を同時に検定する。

  経路1（import 経路の検定）: 原本の onset_probe / relax_probe を**そのまま**呼ぶ。
      公開値と全数値フィールドが厳密一致しなければ、import 経路自体が信用できない。

  経路2（転写の検定）: onset_probe を本ファイル内に一行ずつ転写した
      `onset_probe_transcribed` を呼ぶ。σ₁ 計装はまだ入れない。
      経路1 と厳密一致しなければ、転写が間違っている。

この二つが通って初めて、σ₁ 計装を足す資格ができる。
判定は「厳密一致（== による float 比較）」。許容誤差は置かない。

【CTRL-0 実行で判明したこと（2026-08-18）】
公開 JSON の onset セクションは **stale** であり、固定点として使えない。
現行原本（SHA ba0fc19b）は onset_rate=0.049378 を返すが、JSON は 0.052973。
親残差が 1.38e-5 → 1.82e-9（4桁改善）、parent_rank_planes が 4 → 5 に変わっている。
JSON（2026-07-21 21:45）は make_parent 改良より前の遺物である。
一方 relax セクションは 13/13 フィールドが厳密一致する（relax は make_parent を使わない）。

したがって固定点は **論文値**に取る。
  第8論文 (21614402):        λ = 0.04937,  R² = 0.9999993
  開始様式判別 (21798854):   0.0494/step,  R² = 0.9999998
現行原本 N=5 は 0.049378 / R²=0.9999995 を返し、かつ **seed 非依存**
（seed=0: 0.049378、seed=2: 0.049380）。rate は親平衡の不安定固有値であり
seed の性質ではないので、これが正しい振る舞いである。

N=40 は使わない。現行 onset_probe は seed 依存（seed=0: 0.03612、seed=1: 0.02902）で
論文値 0.0350 と合わない。第8論文の N=40 は白色親（make_parent_white_harmonics）で
構成されており、本 onset_probe の親構成とは別物である。
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.normpath(os.path.join(
    HERE, "..", "..", "自発的分裂予備実験_v1", "run_n_scaling_lowrank_v1.py"))
REF = os.path.join(HERE, "reference_summary_N00005_seed0.json")

# 固定点（論文値）。CTRL-0 の合否はこちらで判定する。
PAPER_RATE = 0.04937          # 第8論文 λ、開始様式判別 0.0494/step
PAPER_RATE_TOL = 5e-5         # 論文の有効桁（4桁）
PAPER_R2_MIN = 0.99999        # 論文 R²=0.9999993 / 0.9999998

# 原本の SHA-256（第7論文リリースノート記載値の完全形）
ORIG_SHA256 = "ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d"

N, SEED = 5, 0
DELTA, CAP = 1e-6, 4000

# 時刻依存のため比較から除くフィールド
TIME_FIELDS = {"t_parent_sec", "t_run_sec", "sec_per_step"}


def load_original():
    import hashlib
    raw = open(ORIG, "rb").read()
    got = hashlib.sha256(raw).hexdigest()
    if got != ORIG_SHA256:
        raise SystemExit(f"原本の SHA-256 不一致\n  期待 {ORIG_SHA256}\n  実際 {got}")
    spec = importlib.util.spec_from_file_location("orig_engine", ORIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # main() は __main__ ガードで走らない
    return mod, got


def onset_probe_transcribed(mod, n, seed, delta=DELTA, cap=CAP):
    """原本 onset_probe（278–335行）の逐語転写。計装は入れない。"""
    sys_lr = mod.LowRankSystem(n)
    rng = np.random.default_rng(40260721 + 1000 * n + seed)
    t0 = time.time()
    v, residual, sig = mod.make_parent(sys_lr, rng)
    t_parent = time.time() - t0

    g = mod.zero_closure_kernel_seed(sys_lr, rng)
    Z = v + delta * g
    Z = Z / np.linalg.norm(Z)
    ztz0 = complex(Z @ Z)

    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    f0 = None
    f_hist = []
    max_ztz = 0.0
    wp = rng.normal(size=sys_lr.m)
    t0 = time.time()
    for t in range(cap + 1):
        h1 = abs(p @ Z) ** 2 + abs(q @ Z) ** 2
        htot = float(np.real(np.conj(Z) @ Z))
        f = 1.0 - h1 / htot
        f_hist.append(f)
        if f0 is None:
            f0 = max(f, 1e-300)
        if f > 0.05:
            break
        if t % 200 == 0 and t > 0:
            mod.progress(f"開始率走行 τ={t} f={f:.3e}")
        max_ztz = max(max_ztz, abs(complex(Z @ Z)))
        sys_lr.set_theta(np.angle(Z))
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
    t_run = time.time() - t0

    f_arr = np.array(f_hist)
    lo, hi = max(10.0 * f0, 1e-13), 1e-3
    mask = (f_arr > lo) & (f_arr < hi)
    idx = np.where(mask)[0]
    rate = None
    if len(idx) >= 5:
        rate = float(np.polyfit(idx, np.log(f_arr[idx]), 1)[0])
    return {
        "parent_residual": residual,
        "parent_sigma2_over_sigma1": float(sig[1] / sig[0]) if len(sig) > 1 else 0.0,
        "parent_rank_planes": int(len(sig)),
        "abs_ztz_initial": abs(ztz0),
        "abs_ztz_max": max_ztz,
        "f_initial": float(f0),
        "steps_run": len(f_hist) - 1,
        "onset_rate_per_step": rate,
        "t_parent_sec": t_parent,
        "t_run_sec": t_run,
    }


def compare(label, got, ref):
    """厳密一致の判定。時刻フィールドは除く。"""
    rows, ok = [], True
    for k in ref:
        if k in TIME_FIELDS:
            continue
        a, b = got.get(k, "<欠落>"), ref[k]
        same = (a == b)
        ok &= same
        rows.append((k, b, a, same))
    print(f"\n--- {label} ---")
    print(f"{'フィールド':<32}{'公開値':<26}{'再現値':<26}判定")
    for k, b, a, same in rows:
        print(f"{k:<32}{str(b):<26}{str(a):<26}{'一致' if same else '★不一致'}")
    return ok


def main():
    mod, sha = load_original()
    print(f"原本: {ORIG}")
    print(f"SHA-256 照合: OK ({sha[:12]}…)")
    print(f"GAMMA = {mod.GAMMA!r}   (= tan(pi/144), 1step の回転角 = pi/72)")
    ref = json.load(open(REF, encoding="utf-8"))
    assert (ref["n"], ref["seed"]) == (N, SEED)
    if ref["gamma"] != mod.GAMMA:
        raise SystemExit("GAMMA が公開値と違う")

    print("\n[経路1] 原本の関数をそのまま呼ぶ")
    a_onset = mod.onset_probe(N, SEED, delta=DELTA, cap=CAP)
    a_relax = mod.relax_probe(N, SEED, cap=3000)

    print("\n[経路2] 転写した onset_probe を呼ぶ（計装なし）")
    b_onset = onset_probe_transcribed(mod, N, SEED)

    ok1 = compare("経路1 onset vs 公開JSON（stale。参考表示のみ）", a_onset, ref["onset"])
    ok2 = compare("経路1 relax vs 公開JSON（本判定に含む）", a_relax, ref["relax"])
    ok3 = compare("経路2 onset vs 経路1 onset（転写の一致。本判定に含む）",
                  b_onset, {k: v for k, v in a_onset.items()})

    # 固定点判定: 論文値との一致
    r1 = a_onset["onset_rate_per_step"]
    r2_ = b_onset["onset_rate_per_step"]
    ok4 = abs(r1 - PAPER_RATE) <= PAPER_RATE_TOL and abs(r2_ - PAPER_RATE) <= PAPER_RATE_TOL
    print(f"\n--- 固定点 vs 論文値 ---")
    print(f"論文 λ = {PAPER_RATE} (±{PAPER_RATE_TOL})")
    print(f"経路1 rate = {r1!r}  差 {abs(r1-PAPER_RATE):.2e}  {'一致' if abs(r1-PAPER_RATE)<=PAPER_RATE_TOL else '★不一致'}")
    print(f"経路2 rate = {r2_!r}  差 {abs(r2_-PAPER_RATE):.2e}  {'一致' if abs(r2_-PAPER_RATE)<=PAPER_RATE_TOL else '★不一致'}")

    out = {
        "control": "CTRL-0 固定点確保",
        "original_path": ORIG,
        "original_sha256": sha,
        "gamma": mod.GAMMA,
        "n": N, "seed": SEED, "delta": DELTA, "cap": CAP,
        "path1_onset_vs_stale_json_exact": bool(ok1),
        "stale_json_note": "公開JSONのonsetはmake_parent改良前の遺物。固定点に使わない",
        "path1_relax_vs_json_exact": bool(ok2),
        "path2_transcription_matches_path1": bool(ok3),
        "rate_matches_paper_lambda": bool(ok4),
        "paper_rate": PAPER_RATE,
        "all_pass": bool(ok2 and ok3 and ok4),
        "path1_onset": a_onset, "path1_relax": a_relax,
        "path2_onset": b_onset,
    }
    p = os.path.join(HERE, "result_control_reproduction_v1.json")
    json.dump(out, open(p, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    print(f"\n出力: {p}")
    print("\n" + ("=== 全判定 PASS（relax厳密一致・転写一致・論文λ一致）。固定点を確保した。 ==="
                  if out["all_pass"] else
                  "=== 判定 FAIL。固定点が取れていない。計装に進んではならない。 ==="))
    sys.exit(0 if out["all_pass"] else 1)


if __name__ == "__main__":
    main()
