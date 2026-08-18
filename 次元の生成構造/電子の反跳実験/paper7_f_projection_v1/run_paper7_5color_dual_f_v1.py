#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第7論文 5色占有時系列の「両形式出力」版。

  python3 run_paper7_5color_dual_f_v1.py 5 [40] [300]

【原本】
  時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/
    exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_5color_timeseries.py
  sha256 fe5c7cbc33437890a5f50944cbbae1594e5f647739d4955402b578f515658503
  （第8論文 Stage A2a expected_hashes.json 登録値。実行時に照合する）

【変更点はこれだけ】
  1. 出力先を本フォルダへ（原本の raw/ を上書きしないため）
  2. **末尾に1列追加**: splitting_fraction_projection
         = |Z - p(p·Z) - q(q·Z)|² / |Z|²
     この式は原本 run() 内に `fval` として既に存在し、crossing 判定に使われている。
     記録されていなかっただけで、新しい量ではない。

  既存16列は式も書式も原本のまま（fmt="%.10e"）。列の順序も変えない。
  したがって元の図はこの CSV からそのまま描ける（列名参照のため末尾追加は無害）。
  追加列だけ %.17e にしてある。10桁では 1e-30 の情報が落ちるため。

【自己検証】
  走行後、既存16列を公開CSV
    paper7_longtime/raw/N{n:05d}/paper7_long_timeseries.csv
  と全行照合する。1行でも違えば異常終了し、追加列も破棄する。
  力学・親構成・サンプリングは一切変更していないので、一致しなければ環境側の問題である。
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
P7 = os.path.join(REPO, "時間軸Q軸とフェルミオンの生成構造", "検証_対照実験",
                  "第5論文原本_自発的分裂予備実験_v1", "exact_lowN_eigenspectrum_v2",
                  "paper7_longtime")
ORIG = os.path.join(P7, "code", "run_paper7_5color_timeseries.py")
HASHES = os.path.join(REPO, "次元の生成構造", "第8論文_二段階seed除去による準安定相の因果分離",
                      "paper8_stage_A2a_seedless_N5", "expected_hashes.json")

BASE_COLUMNS = ["step", "time", "crossing_flag", "splitting_fraction",
                "direction_1_occupation", "direction_2_occupation",
                "direction_3_occupation", "direction_4_occupation",
                "other_rotating_occupation", "kernel_occupation", "occupation_sum",
                "plane_1_occupation", "plane_2_occupation",
                "norm_error", "conservation_error", "projection_closure_error"]
NEW_COLUMN = "splitting_fraction_projection"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def verify_sources():
    exp = json.load(open(HASHES, encoding="utf-8"))
    ok = True
    for group in ("sources", "dependencies"):
        for name, item in exp.get(group, {}).items():
            p = os.path.join(REPO, item["path"])
            got = sha256(p) if os.path.isfile(p) else "<欠落>"
            hit = (got == item["sha256"])
            ok &= hit
            print(f"  {'OK ' if hit else '★NG'} {name}  {got[:16]}…")
    return ok


def load_original():
    spec = importlib.util.spec_from_file_location("p7_5color_orig", ORIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(mod, n):
    """原本 run(n) の逐語転写＋末尾1列。出力先のみ変更。"""
    t0 = time.time()
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = mod.build(n)
    M = sys_lr.m
    occ, s4_new_dirs, align_2d = mod.occ, mod.s4_new_dirs, mod.align_2d
    gram_reduce, dominant_plane = mod.gram_reduce, mod.dominant_plane
    XMAX, SAMPLE = mod.XMAX, mod.SAMPLE

    def fval(Z):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z) @ Z))

    Zc = Z.copy(); wpc = wp.copy(); crossing = None; t = 0
    while True:
        if fval(Zc) > 0.05:
            crossing = t; break
        sys_lr.set_theta(np.angle(Zc)); se, wpc = sys_lr.sigma_max_power(wpc)
        Zc = sys_lr.cayley_step(Zc, se); t += 1

    rows = []
    fmt = "%.10e"
    se_ev = SAMPLE[n]
    f_prev = None
    max_close = 0.0
    Zr = Z.copy(); wpr = wp.copy(); t = 0
    while True:
        if t % se_ev == 0 or t == XMAX:
            totZ = float(np.real(np.conj(Zr) @ Zr))
            E_P1 = occ(B_p1, Zr)
            E_other = occ(B_rot, Zr)
            E_ker = totZ - E_P1 - E_other
            f = 1.0 - E_P1 / totZ
            gr = gram_reduce(sys_lr, Zr)
            _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
            e34 = s4_new_dirs(B0, Bdom)
            proj = B_rot @ (B_rot.T @ e34)
            fq, _ = np.linalg.qr(proj)
            f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
            f34 = align_2d(f_prev, f34); f_prev = f34
            E_d3 = occ(f34[:, [0]], Zr)
            E_d4 = occ(f34[:, [1]], Zr) if f34.shape[1] > 1 else 0.0
            E_rem_other = max(0.0, E_other - E_d3 - E_d4)
            E_a1 = occ(B_p1[:, [0]], Zr); E_a2 = occ(B_p1[:, [1]], Zr)
            osum = (E_P1 + E_d3 + E_d4 + E_rem_other + E_ker) / totZ
            close = abs(osum - 1.0)
            max_close = max(max_close, close)
            plane1 = E_P1 / totZ
            plane2 = (E_d3 + E_d4) / totZ
            rows.append([str(t), str(t), str(int(t >= crossing)), fmt % f,
                         fmt % (E_a1 / totZ), fmt % (E_a2 / totZ),
                         fmt % (E_d3 / totZ), fmt % (E_d4 / totZ),
                         fmt % (E_rem_other / totZ), fmt % (E_ker / totZ), fmt % osum,
                         fmt % plane1, fmt % plane2,
                         fmt % abs(totZ - 1.0), fmt % abs(totZ - 1.0), fmt % close,
                         "%.17e" % fval(Zr)])          # ← 追加列
        if t >= XMAX:
            break
        sys_lr.set_theta(np.angle(Zr)); se, wpr = sys_lr.sigma_max_power(wpr)
        Zr = sys_lr.cayley_step(Zr, se); t += 1

    # --- 自己検証: 既存16列が公開CSVと全行一致するか ---
    ref_path = os.path.join(P7, "raw", f"N{n:05d}", "paper7_long_timeseries.csv")
    verified = None
    if os.path.isfile(ref_path):
        ref = list(csv.reader(open(ref_path, encoding="utf-8")))[1:]
        bad = [i for i, (a, b) in enumerate(zip(ref, rows)) if a != b[:16]]
        verified = (len(ref) == len(rows)) and not bad
        print(f"[N={n}] 公開CSV照合: {len(ref)}行 vs {len(rows)}行, 不一致 {len(bad)} → "
              f"{'一致' if verified else '★不一致'}")
        if not verified:
            raise SystemExit(f"REPRODUCTION_FAILED: N={n} の既存16列が公開CSVと一致しない")
    else:
        print(f"[N={n}] 公開CSV が無いため照合スキップ: {ref_path}")

    out = os.path.join(HERE, f"dual_f_timeseries_N{n:05d}.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(BASE_COLUMNS + [NEW_COLUMN]); w.writerows(rows)

    meta = {"N": n, "M": M, "crossing": crossing, "xmax": XMAX, "sample_every": se_ev,
            "rows": len(rows), "max_projection_closure_error": max_close,
            "reproduction_verified_16cols": verified,
            "f_subtraction_first": float(rows[0][3]), "f_projection_first": float(rows[0][16]),
            "f_subtraction_min": min(float(r[3]) for r in rows),
            "f_projection_min": min(float(r[16]) for r in rows),
            "elapsed_sec": time.time() - t0}
    json.dump(meta, open(os.path.join(HERE, f"dual_f_meta_N{n:05d}.json"), "w",
                         encoding="utf-8"), indent=2, ensure_ascii=False)
    print(f"[N={n}] crossing={crossing} 標本={len(rows)} {meta['elapsed_sec']:.1f}s")
    print(f"        引き算形 初期={meta['f_subtraction_first']:.6e} 最小={meta['f_subtraction_min']:.6e}")
    print(f"        射影形   初期={meta['f_projection_first']:.6e} 最小={meta['f_projection_min']:.6e}")
    print(f"        → {os.path.basename(out)}")
    return meta


def main():
    ns = [int(a) for a in sys.argv[1:]] or [5]
    print("=== 原本の照合 ===")
    if not verify_sources():
        raise SystemExit("SOURCE_MISMATCH")
    mod = load_original()
    print(f"\nDELTA={mod.DELTA} XMAX={mod.XMAX} SAMPLE={mod.SAMPLE}\n")
    for n in ns:
        run(mod, n)


if __name__ == "__main__":
    main()
