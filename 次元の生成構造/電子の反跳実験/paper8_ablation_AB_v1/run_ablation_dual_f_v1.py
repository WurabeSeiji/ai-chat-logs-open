#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第8論文 予備実験（条件A=無シード / 条件B=δ=1e-15）の「両形式出力」版。

  python3 run_ablation_dual_f_v1.py 40 A B
  python3 run_ablation_dual_f_v1.py 5 A B          # N=5 でも同手順

【なぜ1本で両系列を扱えるか】
原本 `run_preliminary_seed_ablation_v1.py` の `build_init(n, initial_seed)` が
シードあり／なしを同一関数で切り替える。docstring より:

    条件A: initial seed OFF（Z0 = v。kernel seed 生成で乱数を消費しない）
    条件B: initial seed ON （Z0 = (v+δg)/‖·‖, δ=1e-15）

親状態 v は両条件で同一（rng = default_rng(40260722 + 1000*n)、iters=1200、tol=1e-12）。
条件B は kernel seed 生成で乱数を消費するため、以後の warm-start `wp` が条件A と異なる。
これはシードを置くことの不可避な帰結であって、実装差ではない。

【原本】
  次元の生成構造/第8論文_二段階seed除去による準安定相の因果分離/code/
    run_preliminary_seed_ablation_v1.py
  sha256 75a10a5b951302bef2ba77cf363066fe8cbdce9a7f1a0af70ccd2df4e520b1d8
  （開始様式判別論文 リリースノート記載の「力学コード」SHA と一致。実行時に照合する）

【変更点はこれだけ】
  1. 出力先を本フォルダへ（原本は PAPER8/raw/ を上書きするため直接は呼べない）
  2. **末尾に1列追加**: f_projection = |Z-p(p·Z)-q(q·Z)|²/|Z|²
     この式は原本 run() 内に `fval` として既に存在し、crossing 判定に使われている。
     CSV に記録されていなかっただけで、新しい量ではない。
     既存27列の `f_outside_parent` は `1 - E_P1/totZ`（桁落ちする引き算形）のまま。

  条件D（準安定シード注入）は本版では扱わない（A/B のみ）。

【判定】
  公開済み
    第8論文…/raw/N{n:05d}/condition_{A,B}_*.csv
  の**既存27列**と全行が文字列で一致すること。許容誤差は置かない。
"""

import csv
import hashlib
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
PAPER8 = os.path.join(REPO, "次元の生成構造", "第8論文_二段階seed除去による準安定相の因果分離")
ORIG = os.path.join(PAPER8, "code", "run_preliminary_seed_ablation_v1.py")
ORIG_SHA = "75a10a5b951302bef2ba77cf363066fe8cbdce9a7f1a0af70ccd2df4e520b1d8"

BASE_COLUMNS = ["step", "time", "N", "condition", "initial_seed_enabled",
                "metastable_seed_enabled", "initial_seed_amplitude",
                "metastable_seed_amplitude", "parent_plane_occupation",
                "f_outside_parent", "q1", "q2", "q3", "q4", "rank_Q",
                "dominant_plane_occupation", "non_dominant_occupation",
                "kernel_occupation", "residual_occupation", "norm_Z",
                "dagger_norm_error", "zero_square_real", "zero_square_imag",
                "zero_square_abs", "projection_closure_error",
                "crossing_detected", "metastable_start_detected"]
NEW_COLUMN = "f_projection"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_original():
    got = sha256(ORIG)
    print(f"原本 SHA-256 {'一致' if got == ORIG_SHA else '★不一致'}  {got[:16]}…")
    if got != ORIG_SHA:
        raise SystemExit("SOURCE_MISMATCH")
    spec = importlib.util.spec_from_file_location("abl_orig", ORIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # run() は __main__ ガードで走らない
    return mod


def run(mod, n, cond, se_ev=None, xmax=None, tag=""):
    """原本 run(n, cond) の逐語転写。出力先と、末尾1列の追加のみ。
    se_ev / xmax を指定すると記録間隔・終了stepだけを変える（毎ステップ版で使う）。"""
    assert cond in ("A", "B"), "本版は条件A/Bのみ"
    c = mod.CONDITIONS[cond]
    t0 = time.time()
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = mod.build_init(n, c["initial_seed"])
    M = sys_lr.m
    occ, qsv4 = mod.occ, mod.qsv4
    gram_reduce, dominant_plane, evolve = mod.gram_reduce, mod.dominant_plane, mod.evolve
    XMAX = mod.XMAX if xmax is None else xmax
    GUARD, Q_REL_TAU, D_EPS, DELTA = mod.GUARD, mod.Q_REL_TAU, mod.D_EPS, mod.DELTA

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    rows = []
    fmt = "%.10e"
    if se_ev is None:
        se_ev = mod.SAMPLE[n]
    init_amp = DELTA if c["initial_seed"] else 0.0
    meta_amp = D_EPS if c["metastable_seed"] else 0.0
    dg = {"crossing_step": None, "metastable_start_step": None}

    crossing = None
    t = 0
    while True:
        f = fval(Z)
        if crossing is None and f > 0.05:
            crossing = t; dg["crossing_step"] = t
        t1 = (crossing + GUARD) if crossing is not None else None
        if t1 is not None and dg["metastable_start_step"] is None:
            dg["metastable_start_step"] = t1

        if t % se_ev == 0 or t == XMAX:
            totZ = float(np.real(np.conj(Z) @ Z))
            E_P1 = occ(B_p1, Z); E_other = occ(B_rot, Z); E_ker = totZ - E_P1 - E_other
            gr = gram_reduce(sys_lr, Z); _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
            E_dom = occ(Bdom, Z)
            qs = qsv4(B0, Bdom)
            rankQ = int(np.sum(qs > Q_REL_TAU * qs[0]))
            ztz = complex(Z @ Z)
            closure = abs(totZ - E_P1 - E_other - E_ker)
            residual = closure
            met_start = int(dg["metastable_start_step"] is not None
                            and t >= dg["metastable_start_step"])
            rows.append([str(t), str(t), str(n), cond,
                         str(int(c["initial_seed"])), str(int(c["metastable_seed"])),
                         fmt % init_amp, fmt % meta_amp, fmt % (E_P1 / totZ),
                         fmt % (1 - E_P1 / totZ),
                         fmt % qs[0], fmt % qs[1], fmt % qs[2], fmt % qs[3], str(rankQ),
                         fmt % (E_dom / totZ), fmt % (E_other / totZ), fmt % (E_ker / totZ),
                         fmt % (residual / totZ), fmt % np.sqrt(totZ),
                         fmt % abs(totZ - 1.0), fmt % ztz.real, fmt % ztz.imag,
                         fmt % abs(ztz), fmt % (residual / totZ),
                         str(int(crossing is not None)), str(met_start),
                         "%.17e" % fval(Z)])          # ★ 追加列
        if t >= XMAX:
            break
        Z, wp = evolve(sys_lr, Z, wp); t += 1
    elapsed = time.time() - t0

    # --- 判定: 既存27列が公開CSVと一致するか（記録間隔が公開と同じときのみ全行、
    #            変えたときは公開の step にある行だけを部分照合） ---
    ref_path = os.path.join(PAPER8, "raw", f"N{n:05d}", f"{c['file']}.csv")
    verified, checked, bad = None, 0, []
    if os.path.isfile(ref_path):
        ref = {r[0]: r for r in list(csv.reader(open(ref_path, encoding="utf-8")))[1:]}
        for r in rows:
            if r[0] in ref:
                checked += 1
                if ref[r[0]] != r[:27]:
                    bad.append(r[0])
        verified = (checked > 0) and not bad
        print(f"  [N={n} 条件{cond}{tag}] 公開CSV照合 {checked} 行 × 27列: "
              f"不一致 {len(bad)} → {'一致' if verified else '★不一致'}")
        if not verified:
            raise SystemExit(f"REPRODUCTION_FAILED: N={n} 条件{cond}")
    else:
        print(f"  [N={n} 条件{cond}{tag}] 公開CSV なし: {ref_path}")

    stem = f"ablation_N{n:05d}_cond{cond}{tag}"
    out = os.path.join(HERE, f"{stem}.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(BASE_COLUMNS + [NEW_COLUMN]); w.writerows(rows)

    fp = np.array([float(r[27]) for r in rows])
    fs = np.array([float(r[9]) for r in rows])
    meta = {"N": n, "condition": cond, "M": M, "initial_seed": c["initial_seed"],
            "initial_seed_amplitude": init_amp, "sample_every": se_ev, "xmax": XMAX,
            "rows": len(rows), "crossing": dg["crossing_step"],
            "metastable_start": dg["metastable_start_step"],
            "reference_csv": os.path.relpath(ref_path, REPO),
            "verified_rows": checked, "verified_27cols": verified,
            "f_projection_first": float(fp[0]), "f_projection_min": float(fp.min()),
            "f_subtraction_first": float(fs[0]), "f_subtraction_min": float(fs.min()),
            "elapsed_sec": elapsed}
    json.dump(meta, open(os.path.join(HERE, f"{stem}_meta.json"), "w", encoding="utf-8"),
              indent=2, ensure_ascii=False)
    print(f"        crossing={dg['crossing_step']} 標本={len(rows)} {elapsed:.1f}s")
    print(f"        引き算形 初期={fs[0]:.6e} 最小={fs.min():.6e}")
    print(f"        射影形   初期={fp[0]:.6e} 最小={fp.min():.6e}")
    print(f"        → {os.path.basename(out)}")
    return meta


def main():
    args = sys.argv[1:]
    n = int(args[0]) if args else 40
    conds = [a.upper() for a in args[1:]] or ["A", "B"]
    print("=== 原本の照合 ===")
    mod = load_original()
    print(f"DELTA={mod.DELTA} XMAX={mod.XMAX} SAMPLE={mod.SAMPLE} GUARD={mod.GUARD}\n")
    for cond in conds:
        run(mod, n, cond)


if __name__ == "__main__":
    main()
