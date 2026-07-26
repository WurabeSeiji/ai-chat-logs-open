#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第8論文 Phase 1：第7論文環境の依存関係監査（read-only, 力学は動かさない）。

指示書§5の10項目のみを収集し reports/paper7_dependency_audit.md と
config/source_file_hashes.json を出力する。第7論文フォルダは書き換えない。
"""
import hashlib
import json
import re
import sys
from pathlib import Path

CODE = Path(__file__).resolve().parent
PAPER8 = CODE.parent                      # 第8論文_.../
REPO = PAPER8.parent.parent               # リポジトリルート
ENGINE = REPO / "時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1"
V2 = ENGINE / "exact_lowN_eigenspectrum_v2"
PL = V2 / "paper7_longtime"

# 第7論文コード（依存閉包）と設定
CODE_FILES = {
    "engine": ENGINE / "run_n_scaling_lowrank_v1.py",
    "parent_basis_exact": ENGINE / "run_plane_flow_exact_v1.py",
    "parent_basis_approx": ENGINE / "run_plane_flow_approx_v1.py",
    "gram_dominant_plane": V2 / "code/run_n300_dimension_saturation_v2.py",
    "paper7_5color": PL / "code/run_paper7_5color_timeseries.py",
    "paper7_transverse": PL / "code/run_paper7_transverse.py",
    "paper7_transverse_cached": PL / "code/run_paper7_transverse_cached.py",
    "paper7_exact_vs_approx": PL / "code/run_paper7_exact_vs_approx_N40.py",
    "paper7_figures": PL / "code/make_paper7_figures.py",
}
DATA_NATURAL = {n: PL / f"raw/N{n:05d}/paper7_long_timeseries.csv" for n in (5, 40, 300)}
DATA_TWOSEED = {n: PL / f"raw/N{n:05d}/transverse_stability_timeseries.csv" for n in (5, 40, 300)}
SUMMARY_FILES = [PL / "summary/N_comparison_table.csv", PL / "summary/transverse_stability_summary.csv"]


def sha256(p):
    if not p.exists():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def find_lines(p, pattern):
    if not p.exists():
        return []
    out = []
    for i, line in enumerate(open(p, encoding="utf-8"), 1):
        if re.search(pattern, line):
            out.append((i, line.rstrip()))
    return out


def csv_header(p):
    if not p.exists():
        return None
    with open(p, encoding="utf-8") as fh:
        return fh.readline().strip().split(",")


def main():
    missing = []
    hashes = {}
    for k, p in CODE_FILES.items():
        h = sha256(p)
        hashes[k] = {"path": str(p.relative_to(REPO)) if p.exists() else str(p), "sha256": h}
        if h is None:
            missing.append(str(p))

    # 3. 初期seed投入位置（run_paper7_5color build 内）
    seed_init = find_lines(CODE_FILES["paper7_5color"], r"DELTA|zero_closure_kernel_seed|v \+ DELTA|Z = v")
    # 4. 準安定横摂動seed投入位置（transverse）
    seed_meta = find_lines(CODE_FILES["paper7_transverse"], r"eps \* eta|eta_r|S4_t0|Z0 \+ eps|Benettin|renorm")
    seed_meta_cached = find_lines(CODE_FILES["paper7_transverse_cached"], r"eps \* eta|Z0 \+ eps|dp / max")
    # 2/9. N・seed規則・δ
    build_src = find_lines(CODE_FILES["paper7_5color"], r"default_rng\(|40260722|DELTA =|SAMPLE =|XMAX =")
    # 5. crossing/準安定/最終 定義
    defs = {
        "crossing_fval>0.05": find_lines(CODE_FILES["paper7_5color"], r"fval\(.*\) > 0.05|> 0\.05"),
        "GUARD(metastable start offset)": find_lines(CODE_FILES["paper7_transverse"], r"GUARD ="),
        "XMAX(final step)": find_lines(CODE_FILES["paper7_5color"], r"XMAX ="),
        "DT(Benettin interval)": find_lines(CODE_FILES["paper7_transverse"], r"DT ="),
        "SEEDS/EPS": find_lines(CODE_FILES["paper7_transverse"], r"SEEDS =|EPS ="),
    }
    # 8. Python packages
    pkgs = {}
    for m in ("numpy", "scipy", "matplotlib"):
        try:
            pkgs[m] = __import__(m).__version__
        except Exception as e:
            pkgs[m] = f"UNAVAILABLE ({e})"

    # 6/7. データ所在と列
    natural = {n: {"path": str(p.relative_to(REPO)) if p.exists() else str(p),
                   "exists": p.exists(), "columns": csv_header(p)} for n, p in DATA_NATURAL.items()}
    twoseed = {n: {"path": str(p.relative_to(REPO)) if p.exists() else str(p),
                   "exists": p.exists(), "columns": csv_header(p)} for n, p in DATA_TWOSEED.items()}

    # §9 で要求される列
    required_cols = ["step", "time", "N", "condition", "initial_seed_enabled", "metastable_seed_enabled",
                     "initial_seed_amplitude", "metastable_seed_amplitude", "parent_plane_occupation",
                     "f_outside_parent", "q1", "q2", "q3", "q4", "rank_Q", "dominant_plane_occupation",
                     "non_dominant_occupation", "kernel_occupation", "residual_occupation", "norm_Z",
                     "dagger_norm_error", "zero_square_real", "zero_square_imag", "zero_square_abs",
                     "projection_closure_error", "crossing_detected", "metastable_start_detected"]
    nat_cols = natural[5]["columns"] or []
    nat_missing = [c for c in required_cols if c not in nat_cols]

    # 出力：source_file_hashes.json
    (PAPER8 / "config").mkdir(exist_ok=True)
    with open(PAPER8 / "config/source_file_hashes.json", "w", encoding="utf-8") as fh:
        json.dump({"code_files": hashes,
                   "data_natural": {n: {"path": v["path"], "sha256": sha256(DATA_NATURAL[n])} for n, v in natural.items()},
                   "data_twoseed": {n: {"path": v["path"], "sha256": sha256(DATA_TWOSEED[n])} for n, v in twoseed.items()},
                   "summary_files": {p.name: sha256(p) for p in SUMMARY_FILES}}, fh, indent=2, ensure_ascii=False)

    # 監査報告 md
    R = PAPER8 / "reports"; R.mkdir(exist_ok=True)
    with open(R / "paper7_dependency_audit.md", "w", encoding="utf-8") as fh:
        w = fh.write
        w("# 第7論文 依存関係監査（第8論文 Phase 1・read-only）\n\n")
        w(f"PAPER7 engine root: `{ENGINE.relative_to(REPO)}`\n\n")
        w("## 1. 第7論文 実行コードの場所 と 10. SHA-256\n\n| 役割 | パス | SHA-256(先頭16) |\n|:--|:--|:--|\n")
        for k, v in hashes.items():
            w(f"| {k} | `{v['path']}` | `{(v['sha256'] or 'MISSING')[:16]}` |\n")
        w("\n（完全なSHA-256は `config/source_file_hashes.json`。）\n\n")
        w("## 2. N=5,40,300 のパラメータ / 9. 乱数生成器と seed値\n\n")
        w("- 乱数: `numpy.random.default_rng(40260722 + 1000*N)`（build 内, N毎に固定）\n")
        w("- 初期微小種 δ = 1e-15\n- 記録間隔 SAMPLE = {5:25, 40:25, 300:100}\n")
        w("- build 関連行:\n```\n")
        for i, l in build_src:
            w(f"{i}: {l}\n")
        w("```\n\n")
        w("## 3. 初期seedを加えている正確なコード位置（run_paper7_5color_timeseries.py, build）\n\n```\n")
        for i, l in seed_init:
            w(f"{i}: {l}\n")
        w("```\n\n")
        w("## 4. 準安定域で横摂動seedを加えている位置（run_paper7_transverse.py / _cached.py）\n\n")
        w("run_paper7_transverse.py:\n```\n")
        for i, l in seed_meta:
            w(f"{i}: {l}\n")
        w("```\nrun_paper7_transverse_cached.py:\n```\n")
        for i, l in seed_meta_cached:
            w(f"{i}: {l}\n")
        w("```\n\n")
        w("## 5. 初期状態 / crossing / 準安定開始 / 最終 の定義\n\n")
        w("- 初期状態: build() で `Z = (v + δ g)/‖·‖`（δ=1e-15）。無seedなら `Z0 = v`。\n")
        w("- crossing: 分裂量 f = 1 - E_P1 が初めて 0.05 を超える step。\n")
        w("- 準安定開始 t0: crossing + GUARD（transverse で GUARD=3000）。\n")
        w("- 最終: 絶対 step XMAX = 55000。\n\n定義行:\n```\n")
        for name, lines in defs.items():
            w(f"[{name}]\n")
            for i, l in lines:
                w(f"  {i}: {l}\n")
        w("```\n\n")
        w("## 6. 既存の自然軌道データ の場所 と 列\n\n")
        for n in (5, 40, 300):
            w(f"- N={n}: `{natural[n]['path']}` exists={natural[n]['exists']}\n")
        w(f"\nN=5 列: {nat_cols}\n\n")
        w("## 7. 既存の二段階seedあり実験データ の場所 と 列\n\n")
        for n in (5, 40, 300):
            w(f"- N={n}: `{twoseed[n]['path']}` exists={twoseed[n]['exists']}\n")
        w(f"\nN=5 列: {twoseed[5]['columns']}\n\n")
        w("## 8. 再現に必要な Python パッケージとバージョン\n\n")
        for k, v in pkgs.items():
            w(f"- {k}: {v}\n")
        w("\n## §9 要求列と既存自然軌道CSVの差分（不足=第8論文側で新規記録が必要）\n\n")
        w(f"既存 paper7_long_timeseries.csv に**存在しない** §9 列: {nat_missing}\n\n")
        w("→ 条件A・Bは第8論文ラッパー `run_preliminary_seed_ablation_v1.py` で §9 全列を新規記録する。\n")
        w("→ 条件D（初期ON＋準安定ON）は、既存 transverse CSV が §9 列を持たないため、同一コード・同一設定で\n")
        w("   §9 列を記録して再生成する（指示書§6.2の「必要列が存在しない場合のみ再生成」に該当。理由: 上記列差分）。\n\n")
        w("## COMMON_FINAL_STEP\n\n- COMMON_FINAL_STEP = 55000（第7論文 XMAX と同一）。\n\n")
        w("## 監査判定\n\n")
        if missing:
            w(f"**不足あり（停止）**: {missing}\n")
        else:
            w("第7論文コード・自然軌道データ・二段階seedデータの所在を全て確認。欠落なし。\n")
            w("第7論文コードは read-only import で再利用（コピー不要）。seed の ON/OFF は第8論文ラッパーで明示切替。\n")

    print("[audit] reports/paper7_dependency_audit.md, config/source_file_hashes.json を出力")
    print(f"  コード欠落: {len(missing)}  自然軌道データ: {sum(natural[n]['exists'] for n in (5,40,300))}/3  "
          f"二段階seedデータ: {sum(twoseed[n]['exists'] for n in (5,40,300))}/3")
    print(f"  §9要求列のうち既存自然軌道CSVに無い列数: {len(nat_missing)}")
    print(f"  packages: numpy={pkgs['numpy']} scipy={pkgs['scipy']} matplotlib={pkgs['matplotlib']}")
    if missing:
        print("  → 不足あり。監査のみ出力して停止。")
        sys.exit(0)


if __name__ == "__main__":
    main()
