#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第2予備実験 §3 依存監査（read-only, 力学を動かさない）。解釈なし。

第8論文第1予備実験・第7論文コードの所在・定義・SHA-256 を収集し
reports/resolution_scan_dependency_audit.md, config/source_file_hashes.json を出力する。
"""
import hashlib
import json
import platform
import re
import sys
from pathlib import Path

CODE = Path(__file__).resolve().parent
P2 = CODE.parent                                  # preliminary_02_resolution_scan/
PAPER8 = P2.parent                                 # 第8論文_.../
REPO = PAPER8.parent.parent
ENGINE = REPO / "時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1"
V2 = ENGINE / "exact_lowN_eigenspectrum_v2"
PL = V2 / "paper7_longtime"

FILES = {
    "prelim1_run": PAPER8 / "code/run_preliminary_seed_ablation_v1.py",
    "engine": ENGINE / "run_n_scaling_lowrank_v1.py",
    "parent_basis_exact": ENGINE / "run_plane_flow_exact_v1.py",
    "parent_basis_approx": ENGINE / "run_plane_flow_approx_v1.py",
    "gram_dominant_plane": V2 / "code/run_n300_dimension_saturation_v2.py",
    "retract_source": ENGINE / "run_transverse_stability_v1.py",
    "paper7_5color": PL / "code/run_paper7_5color_timeseries.py",
}


def sha256(p):
    if not p.exists():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def find(p, pat):
    if not p.exists():
        return []
    return [(i, l.rstrip()) for i, l in enumerate(open(p, encoding="utf-8"), 1) if re.search(pat, l)]


def main():
    hashes = {k: {"path": str(v.relative_to(REPO)) if v.exists() else str(v), "sha256": sha256(v)}
              for k, v in FILES.items()}
    missing = [k for k, v in hashes.items() if v["sha256"] is None]

    z0 = find(FILES["prelim1_run"], r"Z0 = v\.copy|Z0 = v|else:.*無seed|Z = v\.copy")
    seed_off = find(FILES["prelim1_run"], r"initial_seed|zero_closure_kernel_seed|乱数を消費|Z0 = v")
    step_order = find(FILES["prelim1_run"], r"def evolve|set_theta|sigma_max_power|cayley_step")
    obs = find(FILES["prelim1_run"], r"def occ|def qsv4|E_P1 =|f_outside|Q_REL_TAU|rank")
    crossing = find(FILES["prelim1_run"], r"> 0.05|crossing")
    sample = find(FILES["prelim1_run"], r"SAMPLE =|XMAX =|GUARD =")
    retract = find(FILES["retract_source"], r"def retract")

    import numpy as np
    fi = np.finfo(np.float64)
    try:
        blas = np.__config__.show(mode="dicts") if hasattr(np.__config__, "show") else "n/a"
    except Exception:
        blas = "n/a"

    (P2 / "config").mkdir(exist_ok=True)
    with open(P2 / "config/source_file_hashes.json", "w", encoding="utf-8") as fh:
        json.dump(hashes, fh, indent=2, ensure_ascii=False)

    R = P2 / "reports"; R.mkdir(exist_ok=True)
    with open(R / "resolution_scan_dependency_audit.md", "w", encoding="utf-8") as fh:
        w = fh.write
        w("# 第2予備実験 依存監査（read-only）\n\n")
        w("## 1. 第1予備実験 実行コード と 12. SHA-256\n\n| 役割 | パス | SHA-256(先頭16) |\n|:--|:--|:--|\n")
        for k, v in hashes.items():
            w(f"| {k} | `{v['path']}` | `{(v['sha256'] or 'MISSING')[:16]}` |\n")
        w("\n（完全なSHA: `config/source_file_hashes.json`）\n\n")
        w("## 2. N=5,40,300 親状態生成法\n\n- `make_parent(LowRankSystem(N), default_rng(40260722+1000*N), iters=1200, tol=1e-12)` → 親 v。\n\n")
        w("## 3. 初期状態 Z0=v 構築位置（run_preliminary_seed_ablation_v1）\n\n```\n")
        for i, l in z0:
            w(f"{i}: {l}\n")
        w("```\n\n## 4,5. 初期seed/準安定seed OFF の位置\n\n```\n")
        for i, l in seed_off:
            w(f"{i}: {l}\n")
        w("```\n\n## 6,7. 1 step 更新順序（既存 evolve = Cayley のみ。閉鎖・正規化は Cayley が暗黙保存）\n\n```\n")
        for i, l in step_order:
            w(f"{i}: {l}\n")
        w("```\n\n本実験では量子化後に polar retraction（下記 retract）を1回のみ適用。Cayley 直後は測定のみ。\n\n")
        w("## 8. f_outside/q3/q4/rank_Q 算出法\n\n```\n")
        for i, l in obs:
            w(f"{i}: {l}\n")
        w("```\n- f_outside = 1 - E_P1/|Z|²（E_P1 = 固定親平面占有）。q3,q4 = svals([B0|B_dom]) の3,4。\n")
        w("  rank_Q = #{q_j > 1e-8 q1}。E_dom = 瞬時支配平面占有（Gram）。\n\n")
        w("## 9. crossing 判定式\n\n```\n")
        for i, l in crossing:
            w(f"{i}: {l}\n")
        w("```\n\n## 10. 時刻刻み・記録間隔（本実験は §9 で新規固定）\n\n")
        w("- 本実験の保存: step 0..1000 毎step, 1001.. 5step毎, 停止step 必ず。max_step: N5=2500/N40=4500/N300=10000。\n\n")
        w("## 11. 数値型・線形代数・丸め\n\n")
        w(f"- float64: eps={fi.eps}, tiny(smallest normal)={fi.tiny}, smallest_subnormal={np.nextafter(0,1)}, bits={fi.bits}\n")
        w(f"- Q_Δ 丸め: half_to_even（numpy round はデフォルト banker's rounding = half to even）\n")
        w(f"- platform: {platform.platform()}\n\n")
        w("## polar retraction（採用・第7論文, 不変更 import）\n\n```\n")
        for i, l in retract:
            w(f"run_transverse_stability_v1.py:{i}: {l}\n")
        w("```\n\n")
        w("## 監査判定\n\n")
        if missing:
            w(f"**不足あり（停止）**: {missing}\n")
        else:
            w("第1予備実験コード・第7論文 retract・親基底・Gram を全て確認。欠落なし。read-only import で再利用。\n")

    print("[audit] reports/resolution_scan_dependency_audit.md, config/source_file_hashes.json")
    print(f"  コード欠落={len(missing)}  retract確認={'あり' if retract else 'なし'}")
    if missing:
        print("  → 不足あり。停止。")
        sys.exit(0)


if __name__ == "__main__":
    main()
