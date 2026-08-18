#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
対照実験 CTRL-2: 第7論文 5色占有時系列（N=5）の厳密再現。

  python3 run_control_paper7_5color_v1.py

【原本の特定】
第7論文リリースノートが公開プログラムとして挙げる `run_paper7_5color_timeseries.py`。
所在と SHA-256 は第8論文 Stage A2a の expected_hashes.json に登録済み。

  path   : 時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1/
           exact_lowN_eigenspectrum_v2/paper7_longtime/code/run_paper7_5color_timeseries.py
  sha256 : fe5c7cbc33437890a5f50944cbbae1594e5f647739d4955402b578f515658503
  role   : occ、s4_new_dirs、align_2d

【方針】
`occ` / `build` / `s4_new_dirs` / `align_2d` および力学（LowRankSystem, make_parent,
cayley_step 等）は**原本を import してそのまま使う**。転写するのは `run(n)` 本体だけで、
変更点は出力先を本フォルダに変えることのみ（原本は raw/ を上書きするため呼べない）。

【判定】
公開済み `paper7_longtime/raw/N00005/paper7_long_timeseries.csv`（2202行・16列）と、
**全行全列が文字列として厳密一致**すること。許容誤差は置かない。

【計装（判定には含めない）】
原本 `run()` 内に既に定義されている射影形

    fval(Z) = |Z - p(p·Z) - q(q·Z)|^2 / |Z|^2

は crossing の判定にしか使われず、CSV に記録されていない。CSV の列
`splitting_fraction` は `f = 1.0 - E_P1/totZ`（固定親基底からの引き算＝桁落ちする形）
である。本対照では同じ走行から fval も同時に記録し、別ファイルへ出す。
再現判定は 16列の一致のみで行い、追加記録は判定に影響しない。
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
ORIG = os.path.join(
    REPO, "時間軸Q軸とフェルミオンの生成構造", "検証_対照実験",
    "第5論文原本_自発的分裂予備実験_v1", "exact_lowN_eigenspectrum_v2",
    "paper7_longtime", "code", "run_paper7_5color_timeseries.py")
HASHES = os.path.join(
    REPO, "次元の生成構造", "第8論文_二段階seed除去による準安定相の因果分離",
    "paper8_stage_A2a_seedless_N5", "expected_hashes.json")
REF = os.path.join(
    REPO, "時間軸Q軸とフェルミオンの生成構造", "検証_対照実験",
    "第5論文原本_自発的分裂予備実験_v1", "exact_lowN_eigenspectrum_v2",
    "paper7_longtime", "raw", "N00005", "paper7_long_timeseries.csv")

N = 5
COLUMNS = ["step", "time", "crossing_flag", "splitting_fraction",
           "direction_1_occupation", "direction_2_occupation",
           "direction_3_occupation", "direction_4_occupation",
           "other_rotating_occupation", "kernel_occupation", "occupation_sum",
           "plane_1_occupation", "plane_2_occupation",
           "norm_error", "conservation_error", "projection_closure_error"]


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


def main():
    print("=== 原本の照合（第8論文 expected_hashes.json 登録値） ===")
    if not verify_sources():
        raise SystemExit("SOURCE_MISMATCH")

    spec = importlib.util.spec_from_file_location("p7_5color", ORIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # run() は __main__ ガードで走らない
    print(f"\n原本 import 完了: DELTA={mod.DELTA}  XMAX={mod.XMAX}  SAMPLE={mod.SAMPLE}")

    # ---- 以下、原本 run(n) の逐語転写。出力先のみ本フォルダへ ----
    t0 = time.time()
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = mod.build(N)
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
    print(f"crossing = {crossing}")

    rows, extra = [], []
    fmt = "%.10e"
    se_ev = SAMPLE[N]
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
                         fmt % abs(totZ - 1.0), fmt % abs(totZ - 1.0), fmt % close])
            # 計装（判定外）: 射影形。原本 run() 内に既にある式をそのまま記録するだけ
            extra.append([str(t), "%.17e" % fval(Zr), "%.17e" % f])
        if t >= XMAX:
            break
        sys_lr.set_theta(np.angle(Zr)); se, wpr = sys_lr.sigma_max_power(wpr)
        Zr = sys_lr.cayley_step(Zr, se); t += 1
    print(f"走行 {time.time()-t0:.1f}s、{len(rows)} 標本、閉鎖誤差max={max_close:.1e}")

    out = os.path.join(HERE, "reproduced_paper7_long_timeseries_N00005.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(COLUMNS); w.writerows(rows)

    # ---- 判定: 公開CSVと全行全列の文字列一致 ----
    ref = list(csv.reader(open(REF, encoding="utf-8")))
    ref_head, ref_rows = ref[0], ref[1:]
    mism = []
    if ref_head != COLUMNS:
        mism.append(("header", ref_head, COLUMNS))
    for i, (a, b) in enumerate(zip(ref_rows, rows)):
        if a != b:
            mism.append((i, a, b))
    same_len = (len(ref_rows) == len(rows))
    exact = same_len and not mism

    print(f"\n公開CSV {len(ref_rows)} 行 / 再現 {len(rows)} 行")
    print(f"不一致行数: {len(mism)}")
    for i, a, b in mism[:3]:
        print(f"  行{i}\n    公開 {a}\n    再現 {b}")

    ext = os.path.join(HERE, "instrumented_f_projection_vs_subtraction_N00005.csv")
    with open(ext, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["step", "f_projection", "f_subtraction"])
        w.writerows(extra)

    fp = np.array([float(r[1]) for r in extra])
    fs = np.array([float(r[2]) for r in extra])
    print(f"\n[計装] 射影形 f: 最小 {fp.min():.6e}  初期 {fp[0]:.6e}")
    print(f"[計装] 引き算形 f: 最小 {fs.min():.6e}  初期 {fs[0]:.6e}")

    json.dump({
        "control": "CTRL-2 第7論文 5色占有時系列 N=5 再現",
        "original": ORIG, "original_sha256": sha256(ORIG),
        "rows_reference": len(ref_rows), "rows_reproduced": len(rows),
        "mismatch_count": len(mism), "exact": bool(exact),
        "crossing": crossing, "max_projection_closure_error": max_close,
        "f_projection_min": float(fp.min()), "f_projection_first": float(fp[0]),
        "f_subtraction_min": float(fs.min()), "f_subtraction_first": float(fs[0]),
    }, open(os.path.join(HERE, "result_control_paper7_5color_v1.json"), "w",
            encoding="utf-8"), indent=2, ensure_ascii=False)

    print("\n" + ("=== PASS。公開CSVと全行全列が厳密一致。固定点を確保した。 ==="
                  if exact else
                  "=== FAIL。再現していない。計装結果は使ってはならない。 ==="))
    sys.exit(0 if exact else 1)


if __name__ == "__main__":
    main()
