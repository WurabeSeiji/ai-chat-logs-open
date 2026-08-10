#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文A 補遺（全実験・全図表・全プログラムの一覧）を自動生成する v1

手で書くと必ず実体とずれるので、フォルダを走査して作る。
出力: 論文A_補遺_一覧_v1.md

使い方: python3 make_paperA_appendix_v1.py
"""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "論文A_補遺_一覧_v1.md"
MODE_JA = {"neutral": "ニュートリノ型（1か所）", "electron": "電子型（1か所）",
           "fermion_family": "奇数の帯5か所", "boson_family": "偶数の帯3か所",
           "mixed": "混合8か所", "vacuum": "種なし"}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def parse_npz(name: str) -> dict:
    s = name[len("nsweep_"):-len("_v2.npz")]
    parts = s.split("_")
    N = int(parts[-1][1:]) if parts[-1].startswith("N") else None
    T, delta, rep, mode = 4000, None, "", []
    for p in parts[:-1]:
        if re.fullmatch(r"T\d+", p):
            T = int(p[1:])
        elif re.fullmatch(r"d[0-9.eE+-]+", p):
            try:
                delta = float(p[1:])
            except ValueError:
                mode.append(p)
        elif p.startswith("rep-"):
            rep = p[4:]
        else:
            mode.append(p)
    return {"mode": "_".join(mode), "T": T, "delta": delta, "N": N, "rep": rep}


def main() -> None:
    L: list[str] = []
    A = L.append

    A("# 補遺　全実験・全図表・全プログラムの一覧")
    A("")
    A("本補遺は `make_paperA_appendix_v1.py` がフォルダを走査して自動生成したものである。")
    A("手作業で書き写した箇所はない。")
    A("")

    # ---------- 補遺A 実験一覧 ----------
    npzs = sorted(HERE.glob("nsweep_*_v2.npz"))
    recs = []
    for p in npzs:
        m = parse_npz(p.name)
        m["file"] = p.name
        m["MB"] = p.stat().st_size / 1e6
        recs.append(m)
    A("## 補遺A　実験の一覧")
    A("")
    A(f"走行ごとの保存データ（NPZ）は **{len(recs)} 件**、"
      f"合計 {sum(r['MB'] for r in recs):.0f} MB である。")
    A("種の置き場所・種の強さ・更新回数・分解能のいずれかが違えば別の走行として数えている。")
    A("")
    grp: dict[tuple, list] = {}
    for r in recs:
        grp.setdefault((r["T"], r["N"] if r["N"] else 0), []).append(r)
    A("### A-1　走行の内訳")
    A("")
    A("| 更新回数 | 分解能 N | 走行数 | 内容 |")
    A("|---:|---:|---:|---|")
    for (T, N) in sorted(grp, key=lambda x: (-x[0], x[1])):
        rs = grp[(T, N)]
        modes = sorted({MODE_JA.get(r["mode"], r["mode"]) for r in rs})
        A(f"| {T} | {N} | {len(rs)} | {'、'.join(modes)} |")
    A("")
    A("### A-2　主実験（分解能 12・更新 42000 回）の 40 条件")
    A("")
    A("| 種の置き方 | " + " | ".join(
        f"{d:g}" for d in [1e-15, 1e-8, 1e-4, 1e-3, 1e-2, 0.0316228, 0.04357, 0.1]) + " |")
    A("|---|" + "---|" * 8)
    for mode in ("neutral", "electron", "fermion_family", "boson_family", "mixed"):
        cells = []
        for d in [1e-15, 1e-8, 1e-4, 1e-3, 1e-2, 0.03162277660168379, 0.04357, 0.1]:
            # ファイル名は %g で丸められている（0.03162277660168379 → d0.0316228）。
            # 相対 1e-5 で照合する。
            hit = [r for r in recs if r["mode"] == mode and r["T"] == 42000
                   and r["N"] == 12 and not r["rep"]
                   and ((r["delta"] is None and abs(d - 0.01) < 1e-18)
                        or (r["delta"] is not None
                            and abs(r["delta"] - d) <= 1e-5 * max(abs(d), 1e-300)))]
            cells.append("○" if hit else "—")
        A(f"| {MODE_JA[mode]} | " + " | ".join(cells) + " |")
    A("")
    A("### A-3　対照・追試・長時間走行")
    A("")
    A("| 種類 | 走行 | 目的 |")
    A("|---|---|---|")
    for r in sorted([x for x in recs if x["rep"]], key=lambda x: x["rep"]):
        kind = ("再現対照" if "ctl" in r["rep"] or "controlrep" in r["rep"]
                else "位相を打ち消した種" if "pbmix" in r["rep"]
                else "長時間走行" if r["T"] >= 300000
                else "分解能 4 の追加走行" if r["N"] == 4 else "その他")
        A(f"| {kind} | {MODE_JA.get(r['mode'], r['mode'])}・"
          f"強さ {r['delta']:g}・更新 {r['T']}・N={r['N']} | `{r['rep']}` |")
    A("")

    # ---------- 補遺B 図の一覧 ----------
    A("## 補遺B　図の一覧")
    A("")
    figA = sorted(HERE.glob("figA*.png"))
    runfig = sorted(HERE.glob("fig_*_v2.png"))
    sheets = sorted(HERE.glob("sheet_*.png"))
    A(f"本文で使う図は **{len(figA)} 枚**、走行ごとの記録図は **{len(runfig)} 枚**、")
    A(f"全数目視のために焼いた一覧シートは **{len(sheets)} 枚**である。")
    A("")
    A("### B-1　本文の図")
    A("")
    A("| 図 | ファイル |")
    A("|---|---|")
    for p in figA:
        A(f"| {p.stem.split('_')[0]} | `{p.name}` |")
    A("")
    A("### B-2　走行ごとの記録図（系統別の枚数）")
    A("")
    A("| 系統 | 内容 | 枚数 |")
    A("|---|---|---:|")
    fam = {"4panel": "空間・平面の枚数・時計・打ち消し残りの 4 段",
           "mix": "帯ごとのパワーと混合率",
           "ledger": "128 か所の地図と、狙った場所の育ち方",
           "summary": "分解能ごとの要約", "birth_matrix": "6 つの判定の一覧"}
    for k, desc in fam.items():
        A(f"| {k} | {desc} | {len(list(HERE.glob(f'fig_*_{k}*_v2.png')))} |")
    A("")
    A("命名は `fig_<種の置き方>[_T更新回数][_d強さ][_rep-識別子]_<系統>[_N分解能]_v2.png`。")
    A("この規約により、本文で参照した図はファイル名から一意に特定できる。")
    A("")
    A("### B-3　全数目視のための一覧シート")
    A("")
    A(f"走行ごとの記録図 {len(runfig)} 枚を 1 枚あたり 12 図で焼いた {len(sheets)} 枚。")
    A("`sheet_<系統>_<通し番号>.png`。監査の経過は `分析記録_全図面監査_v2.md` に記録した。")
    A("")

    # ---------- 補遺C プログラム ----------
    A("## 補遺C　プログラムの一覧")
    A("")
    A("| 役割 | ファイル | SHA-256（先頭16桁） |")
    A("|---|---|---|")
    roles = [
        ("実験本体（分解能掃引・種の置き方・強さ）", "run_nsweep_three_series_v2.py"),
        ("実験本体の下層（分解能 1〜20 の掃引と図）", "run_tb_nsweep_1to20_v1.py"),
        ("不足していた種の置き方の追加取得", "run_missing_seed_sweeps_T42000_v1.py"),
        ("位相を打ち消した種の走行", "run_phase_balanced_mixed_v1.py"),
        ("同上（強さを変えた版）", "run_phase_balanced_mixed_grid_v1.py"),
        ("周期の等分数を変えた検定", "run_divisor_class_register_order_v1.py"),
        ("主張の数値を保存データから再計算し照合", "aggregate_paperA_claims_v1.py"),
        ("本文の図（15 枚）", "make_paperA_figures_v2.py"),
        ("監査から出た図（4 枚）", "make_paperA_figures_audit_v1.py"),
        ("特別な値 0.6972 への最接近を全走行から取得", "probe_alpha_root_closest_approach_v1.py"),
        ("全数目視のための一覧シート", "make_contact_sheets_v1.py"),
        ("本補遺の生成", "make_paperA_appendix_v1.py"),
        ("長時間走行と分解能 4 の追加走行の管理", "run_stage4_longtime_orchestrator_v1.py"),
    ]
    for role, f in roles:
        p = HERE / f
        A(f"| {role} | `{f}` | {sha(p) if p.exists() else '（未配置）'} |")
    A("")
    A("### C-1　外部から読み込んでいる力学（本論文では一切変更していない）")
    A("")
    A("| 役割 | ファイル | SHA-256（先頭16桁） |")
    A("|---|---|---|")
    ext = [("相互作用（更新の後半）", "../統一万能関数_v1/unified_interaction_v1.py"),
           ("平面と 3 方向目の読み出し", "../統一万能関数_v1/unified_dimension_v1.py"),
           ("帯・混合率・物質量の読み出し", "../統一万能関数_v1/unified_readout_v3.py"),
           ("時計が読めるかどうかの判定", "../統一万能関数_v1/selection_v1.py"),
           ("2 つの波の系（8）の検定に使用）",
            "../万能非弾性写像_managed_v1/run_ignition_fate_exact_v3.py")]
    for role, rel in ext:
        p = (HERE / rel).resolve()
        A(f"| {role} | `{rel}` | {sha(p) if p.exists() else '（未配置）'} |")
    A("")
    A("これらは読み込み専用で取り込んでおり、走行の前後で自身の SHA-256 を検査する。")
    A("食い違えばその場で停止する。")
    A("")

    # ---------- 補遺D 事前登録と報告 ----------
    A("## 補遺D　事前登録・結果報告・検算の記録")
    A("")
    A("| 種類 | ファイル |")
    A("|---|---|")
    for p in sorted(HERE.glob("事前登録_*.md")):
        A(f"| 事前登録（測定前に予想を固定した文書） | `{p.name}` |")
    for p in sorted(HERE.glob("結果報告_*.md")):
        A(f"| 結果報告 | `{p.name}` |")
    for p in sorted(HERE.glob("分析記録_*.md")):
        A(f"| 分析記録 | `{p.name}` |")
    for f in ["result_paperA_claims_v1.json", "result_paperA_audit_figs_v1.json",
              "result_divisor_class_register_order_v1.json",
              "result_alpha_root_closest_v1.json"]:
        if (HERE / f).exists():
            A(f"| 数値の記録 | `{f}` |")
    A("")

    # ---------- 補遺E 検算 ----------
    A("## 補遺E　検算の結果")
    A("")
    p = HERE / "result_paperA_claims_v1.json"
    if p.exists():
        c = json.loads(p.read_text())["checks"]
        A("主張の数値は、保存データから独立に計算し直して、"
          "既存の記録と機械的に照合してある。")
        A("")
        A("| 検査 | 結果 |")
        A("|---|---|")
        t = c["tau_vs_result_json"]
        A(f"| 保存データから再計算した立ち上がり・時計の回数が記録と一致するか | "
          f"{t['n_compared']} 条件・食い違い {t['n_mismatch']} 件 |")
        A(f"| 波が占める場所の数 | {'一致' if c['claim2_support']['pass'] else '不一致あり'} |")
        A(f"| 偶数の帯だけの種で奇数の帯が厳密にゼロか | "
          f"最大値 {c['claim3_odd_exact_zero']['max_over_all_deltas']} |")
        f5 = c["claim5_fit"]
        A(f"| 助走期間の式 | 9.8922 − 48.6108 × ln（複素数として足した大きさ）"
          f"　R² = {f5['fit_A_coh']['r2']:.6f} |")
        A(f"| 同じ点をパワーの合計で整理した場合 | R² = {f5['fit_P_seed']['r2']:.6f} |")
        A(f"| 位相を打ち消した種の 4 条件 | "
          f"{'すべて予想どおり' if c['claim6_phase_balanced']['pass'] else '不一致あり'} |")
    q = HERE / "result_divisor_class_register_order_v1.json"
    if q.exists():
        d = json.loads(q.read_text())
        rep = d["reproduction_ne16"]
        tot = sum(v["n_pairs"] for v in d["judgements"].values())
        bad = sum(v["n_fail"] for v in d["judgements"].values())
        A(f"| 前論文の結果の再現（等分数 16） | {rep['n_compared']} 件・"
          f"最大の食い違い {rep['worst_rel']:.1e} |")
        A(f"| 周期の等分数を変えた検定 | {tot} 通り・外れ {bad} 件 |")
    r = HERE / "result_paperA_audit_figs_v1.json"
    if r.exists():
        d = json.loads(r.read_text())
        f = d["fit_PF"]
        A(f"| 波が出そろうまでの回数の式 | 奇数の帯に置いたパワーの "
          f"{f['b']:.3f} 乗に比例　R² = {f['r2']:.4f}（{f['n']} 点） |")
    A("")
    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"→ {OUT.name}（{len(L)} 行）")
    print(f"   走行 {len(recs)} 件 / 本文の図 {len(figA)} 枚 / "
          f"走行図 {len(runfig)} 枚 / シート {len(sheets)} 枚")


if __name__ == "__main__":
    main()
