#!/usr/bin/env python3
"""v4（純進行波）単体テスト＋カタログ閉塞監査＋修正分類器での再分析 v1

Part 1  v4 単体テスト（実行前固定）:
    T1 行閉塞 < 1e-12（全行）
    T2 DC・Nyquist 消滅（|c|<1e-14）・全セクター場の閉塞 0
    T3 全体閉塞 < 1e-12
    T4 親ベクトルが v2/v3 と bitwise 一致（乱数消費順の保存）
    T5 対積負債 全 |k| で 0（純進行波の要件——v4 の存在理由）
    T6 方向均衡（記録）: 海全体の +k/−k 占有比

Part 2  カタログ閉塞監査（閉塞禁止の一般則）:
    基底 qλ₀ の bin は k=N/q。自己対 ⟺ 2k≡0 (mod N) ⟺ q∈{1,2}（q=2は偶数Nのみ実現）
    倍音次数 j（波長 q/j λ₀、bin jN/q）の自己対 ⟺ q | 2j
    → 各族について: 基底が禁止 or 選択倍音数が「非禁止の利用可能次数」の数を超える
      場合、その族は単独閉塞不能（存在資格なし）と判定
    → 監査済みカタログ（catalog_closure_audited_N{5,40}.csv）を管理下に出力
      （原本カタログは無改変）

Part 3  修正分類器 × 監査済みカタログで再分析:
    白色v4 / 単一波版（v·e^{2πin/N}/√N）を分類（禁止族は照合対象から除外）
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CENSUS = HERE.parent / "standalone_parent_census_v1"


def load(name, alias):
    spec = importlib.util.spec_from_file_location(alias, HERE / name)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


orig = load("classifier_original_v2.py", "clf_o4")
orig.HERE = CENSUS
g2 = load("make_parent_white_harmonics_n_only_v2.py", "g2_4")
g4 = load("make_parent_white_harmonics_n_only_v4.py", "g4_4")
fixed = load("stable_wave_classification_fixed_v1.py", "fx_4")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def parse_orders(text: str) -> list[int]:
    return [int(x) for x in re.findall(r"\d+", text)] if text and text != "なし" else []


def audit_catalog(n: int):
    """カタログの各族に閉塞監査列を追加した管理下コピーを作る。"""
    path = orig.find_catalog(n)
    with path.open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    audited = []
    for row in rows:
        q = orig.parse_base_q(row["基底波長"])
        base_forbidden = q in (1, 2) and (q == 1 or n % 2 == 0)
        odd_avail = parse_orders(row["利用可能な奇数倍音次数"])
        even_avail = parse_orders(row["利用可能な偶数倍音次数"])
        odd_ok = [j for j in odd_avail if (2 * j) % q != 0]
        even_ok = [j for j in even_avail if (2 * j) % q != 0]
        odd_forbidden = [j for j in odd_avail if (2 * j) % q == 0]
        even_forbidden = [j for j in even_avail if (2 * j) % q == 0]
        need_odd = int(row["選択する奇数倍音数"])
        need_even = int(row["選択する偶数倍音数"])
        harmonics_infeasible = need_odd > len(odd_ok) or need_even > len(even_ok)
        forbidden = base_forbidden or harmonics_infeasible
        reason = []
        if base_forbidden:
            reason.append(f"基底{q}λ₀が自己対(2k≡0)")
        if harmonics_infeasible:
            reason.append(
                f"選択倍音が禁止次数を要求(奇:禁止{odd_forbidden}/可{odd_ok}, "
                f"偶:禁止{even_forbidden}/可{even_ok})")
        row2 = dict(row)
        row2["閉塞監査"] = "禁止" if forbidden else "可"
        row2["閉塞監査理由"] = "; ".join(reason) if reason else "全構成が単独零閉塞可能"
        audited.append(row2)
    out = HERE / f"catalog_closure_audited_N{n}.csv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audited[0].keys()))
        w.writeheader(); w.writerows(audited)
    ok_keys = {(orig.parse_base_q(r["基底波長"]), int(r["選択する奇数倍音数"]),
                int(r["選択する偶数倍音数"])): r
               for r in audited if r["閉塞監査"] == "可"}
    n_forbidden = sum(1 for r in audited if r["閉塞監査"] == "禁止")
    print(f"  [カタログ監査 N={n}] 全{len(audited)}族 → 可 {len(audited)-n_forbidden} / 禁止 {n_forbidden}")
    for r in audited:
        if r["閉塞監査"] == "禁止":
            print(f"    禁止: {r['族ID']}（{r['基底波長']}）: {r['閉塞監査理由']}")
    return ok_keys, audited


def main() -> None:
    results = {"generator_v4_sha256": sha256(HERE / "make_parent_white_harmonics_n_only_v4.py")}
    for n, seed in ((5, 2), (40, 1)):
        print(f"\n===== N={n} seed={seed} =====")
        r4 = g4.make_parent(n, seed=seed)
        W = r4.relation_waves
        m = W.shape[0]
        C = np.fft.fft(W, axis=1) / n

        rowc = float(np.max(np.abs(np.sum(W * W, axis=1)) / np.sum(np.abs(W) ** 2, axis=1)))
        t1 = rowc < 1e-12
        dc = float(np.max(np.abs(C[:, 0])))
        nyq = float(np.max(np.abs(C[:, n // 2]))) if n % 2 == 0 else 0.0
        t2 = dc < 1e-14 and nyq < 1e-14
        tot = abs(complex(np.sum(W * W)))
        t3 = tot < 1e-12
        r2 = g2.make_parent(n, seed=seed)
        t4 = bool(np.array_equal(r4.parent_vector, r2.parent_vector))
        debts = np.array([[abs(2.0 * C[i, k] * C[i, n - k]) for k in range(1, (n - 1) // 2 + 1)]
                          for i in range(m)])
        t5 = float(debts.max()) < 1e-14   # FFT丸め水準（未占有binの残渣~1e-17との積）
        n_half = (n - 1) // 2
        plus_frac = float(np.mean([np.sum(np.abs(C[i, 1:n_half + 1]) ** 2)
                                    / np.sum(np.abs(C[i]) ** 2) for i in range(m)]))
        print(f"  T1行閉塞 {'PASS' if t1 else 'FAIL'}({rowc:.1e}) "
              f"T2 DC/Nyq {'PASS' if t2 else 'FAIL'}({dc:.1e}/{nyq:.1e}) "
              f"T3全体 {'PASS' if t3 else 'FAIL'}({tot:.1e}) "
              f"T4親bitwise {'PASS' if t4 else 'FAIL'} "
              f"T5対積負債 {'PASS' if t5 else 'FAIL'}(max={debts.max():.1e}) "
              f"T6 +k側パワー比={plus_frac:.3f}")
        results[f"N{n}_tests"] = {"T1": t1, "T2": t2, "T3": t3, "T4": t4,
                                    "T5": t5, "row_closure": rowc,
                                    "pair_debt_max": float(debts.max()),
                                    "plus_fraction": plus_frac}

        ok_keys, audited = audit_catalog(n)
        results[f"N{n}_catalog"] = {"total": len(audited),
                                      "allowed": len(ok_keys),
                                      "forbidden": [r["族ID"] for r in audited
                                                     if r["閉塞監査"] == "禁止"]}

        # 修正分類器 ×監査済みカタログ
        catalog = {k: {"族ID": v["族ID"]} for k, v in ok_keys.items()}
        v = r4.parent_vector
        single = np.outer(v, np.exp(2j * np.pi * np.arange(n) / n) / math.sqrt(n))
        for label, waves in (("白色v4", W), ("単一波版", single)):
            agg, _ = fixed.classify_state(waves, n, catalog, label)
            results[f"N{n}_{label}"] = agg

    (HERE / "test_audit_v4_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print("\nsaved: test_audit_v4_result_v1.json / catalog_closure_audited_N5.csv / N40.csv")


if __name__ == "__main__":
    main()
