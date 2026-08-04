#!/usr/bin/env python3
"""安定波分類器 修正版 v1：監査4欠陥（A〜D）の修正と初期値分析

原本: classifier_original_v2.py（= standalone_parent_census_v1/
    particle_table_white_harmonics_n_only_v2.py、SHA 15ebb7e8… 一致確認済み・read-only）

監査で確定した欠陥と修正:
    A: bin0 を「1λ₀ の N点節エイリアス」として許容波長 q=1 に算入していた。
       → 修正: どちらの解釈（真のDC / k=N の λ₀ 波）でも 2k≡0 (mod N) の自己対で
         単独零閉塞が不可能（存在資格なし）。曖昧さを裁く必要はなく一律
         「closure_forbidden」に分類。
    B: 偶数 N の k=N/2（ナイキスト）を通常の許容波長 2λ₀ としていた。
       → 同じ定理により closure_forbidden。
    C: 定常波表との照合なのに ±k の関係を無検査（位相は記録のみ）。
       → 発見: 同一関係波内の ±k 均衡（定常対）は閉塞と両立しない
         （Σ(c₊e^{+}+c₋e^{−})² = 2N c₊c₋ ≠ 0）。よって「関係波内の定常波」は
         公理的に禁止され、定常性は関係横断でしか実現できない。
         修正: 各 |k| の対積負債 debt_k = |2 c₊ c₋|/全パワー を測定・記録
         （閉塞を他の |k| との相殺で借金している度合い）。
         分類の硬い基準にはせず記録（規則化は木原氏の裁定事項）。
    D: 族キーが（基底波長, 奇数倍音数, 偶数倍音数）の個数のみで、
       異なる倍音組が衝突しうる。→ 実際の倍音組を出力し、衝突可能性を警告。

分析対象:
    (0) 対照再現: 原本分類器で白色v2状態を再分類し、コミット済み census と一致確認
    (1) 白色v2（parent_white_harmonics_N5_v2 / N40_v2 正本）
    (2) 白色v3（管理下 DC/Nyquist 除去生成器・同一シード）
    (3) 単一波版: w[m,n] = v[m]·e^{2πin/N}/√N（リング基本波 k=+1・進行波。
        閉塞定理により ±均衡定常波と k=N（λ₀）は単独存在不可のため、
        公理適合の「単一振動数・単一波長」はこの形が唯一）
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CENSUS = HERE.parent / "standalone_parent_census_v1"
COMMITTED = CENSUS / "stable_wave_classification_after_parent_N5_N40_N300_v1"

spec = importlib.util.spec_from_file_location("clf_orig", HERE / "classifier_original_v2.py")
orig = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = orig
spec.loader.exec_module(orig)
orig.HERE = CENSUS   # 分類表カタログの探索パスを原本フォルダへ（ファイルは無改変）
spec2 = importlib.util.spec_from_file_location("gen3_c", HERE / "make_parent_white_harmonics_n_only_v3.py")
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def fixed_classify(spectrum: np.ndarray, n: int, catalog) -> dict:
    """修正版分類（A/B/C'/D 適用）。spectrum = fft(samples)/sqrt(n)。"""
    power = np.abs(spectrum) ** 2
    total = float(power.sum())
    eps = float(np.finfo(spectrum.real.dtype).eps)
    floor = n * eps * max(float(power.max()), np.finfo(float).tiny)
    present = power > floor

    allowed_wl: set[int] = set()
    outside_orders: list[int] = []
    forbidden_power = 0.0
    outside_power = 0.0
    allowed_power = 0.0
    for b in range(n):
        if not present[b]:
            continue
        order = orig.signed_order(b, n)
        if (2 * b) % n == 0:                       # A/B: 自己対（bin0, N/2）
            forbidden_power += float(power[b])
            continue
        if n % abs(order) == 0:
            allowed_wl.add(n // abs(order))
            allowed_power += float(power[b])
        else:
            outside_orders.append(order)
            outside_power += float(power[b])

    # C': 対積負債（|k| ごとの ±k 閉塞借金）
    debts = {}
    for k in range(1, (n - 1) // 2 + 1):
        cp, cm = spectrum[k], spectrum[n - k]
        debts[k] = float(abs(2.0 * cp * cm) / max(total, 1e-300))
    max_debt = max(debts.values()) if debts else 0.0

    status = "分類外"
    family = None
    base_q = None
    harmonic_set: list[int] = []
    if forbidden_power > 0:
        reason = "存在資格なし成分（自己対 bin0/N2）を含む"
    elif outside_orders:
        reason = "安定波分類表外の波長を含む"
    elif not allowed_wl:
        reason = "数値誤差床より上の波長成分がない"
    else:
        base_q = max(allowed_wl)
        if any(base_q % q != 0 for q in allowed_wl):
            reason = "単一基底波の倍音束ではない"
        else:
            harmonic_set = sorted(base_q // q for q in allowed_wl if q != base_q)
            odd_c = sum(o % 2 == 1 for o in harmonic_set)
            even_c = len(harmonic_set) - odd_c
            family = catalog.get((base_q, odd_c, even_c))
            if family is None:
                reason = "波長構成に対応する安定族が分類表にない"
            else:
                status = "安定波分類表に一致"
                reason = "基底波長と倍音波長の組が一致（注意: 族キーは倍音個数のみで縮退可能性あり）"

    return {"status": status, "reason": reason,
            "base_q": base_q, "harmonic_set": harmonic_set,
            "family_id": family["族ID"] if family else None,
            "allowed_power_frac": allowed_power / total,
            "outside_power_frac": outside_power / total,
            "forbidden_power_frac": forbidden_power / total,
            "max_pair_debt": max_debt, "pair_debts": debts}


def classify_state(waves: np.ndarray, n: int, catalog, label: str):
    m = waves.shape[0]
    rows = [fixed_classify(np.fft.fft(waves[i]) / math.sqrt(n), n, catalog)
            for i in range(m)]
    counts = defaultdict(int)
    for r in rows:
        counts[r["status"]] += 1
    reasons = defaultdict(int)
    for r in rows:
        if r["status"] != "安定波分類表に一致":
            reasons[r["reason"]] += 1
    fam = defaultdict(int)
    for r in rows:
        if r["family_id"]:
            fam[r["family_id"]] += 1
    agg = {
        "label": label, "M": m,
        "counts": dict(counts), "reasons": dict(reasons),
        "families": dict(fam),
        "mean_allowed_power": float(np.mean([r["allowed_power_frac"] for r in rows])),
        "mean_outside_power": float(np.mean([r["outside_power_frac"] for r in rows])),
        "mean_forbidden_power": float(np.mean([r["forbidden_power_frac"] for r in rows])),
        "max_pair_debt": float(max(r["max_pair_debt"] for r in rows)),
        "mean_max_pair_debt": float(np.mean([r["max_pair_debt"] for r in rows])),
    }
    print(f"  [{label}] M={m} 判定={dict(counts)}")
    print(f"    理由={dict(reasons)}")
    if fam:
        print(f"    一致族={dict(fam)}")
    print(f"    パワー配分: 許容={agg['mean_allowed_power']:.4f} 表外={agg['mean_outside_power']:.4f} "
          f"存在資格なし={agg['mean_forbidden_power']:.4f}  対積負債max={agg['max_pair_debt']:.4f}"
          f"（平均max={agg['mean_max_pair_debt']:.4f}）")
    return agg, rows


def control_reproduction(n: int, waves: np.ndarray, catalog) -> bool:
    """原本分類器で committed census の summary を再現。"""
    committed = json.loads((COMMITTED / f"N{n}" / "census.json").read_text(encoding="utf-8"))
    s = committed["summary"]
    stat = defaultdict(int)
    reasons = defaultdict(int)
    a_pow = o_pow = t_pow = 0.0
    for i in range(waves.shape[0]):
        spectrum = np.fft.fft(waves[i]) / math.sqrt(n)
        res, comps = orig.classify_spectrum(spectrum, n, waves.shape[0], catalog)
        stat[res["classification_status"]] += 1
        if res["classification_status"] != "安定波分類表に一致":
            reasons[res["classification_reason"]] += 1
        for c in comps:
            t_pow += c["power"]
            if c["present_above_numerical_floor"] and c["stationary_wavelength_allowed"]:
                a_pow += c["power"]
            elif c["present_above_numerical_floor"]:
                o_pow += c["power"]
    ok = (stat.get("安定波分類表に一致", 0) == s["stable_classified_wave_count"]
          and stat.get("分類外", 0) == s["unclassified_wave_count"]
          and dict(reasons) == s["classification_reason_counts"]
          and abs(a_pow / t_pow - s["allowed_stationary_wavelength_power"]) < 1e-9)
    print(f"  [対照 N={n}] 一致={stat.get('安定波分類表に一致',0)} 分類外={stat.get('分類外',0)} "
          f"許容パワー比={a_pow/t_pow:.12f}（committed {s['allowed_stationary_wavelength_power']:.12f}）"
          f" → {'PASS' if ok else 'FAIL'}")
    return ok


def main() -> None:
    out = {"classifier_original_sha256": sha256(HERE / "classifier_original_v2.py"),
           "generator_v3_sha256": sha256(HERE / "make_parent_white_harmonics_n_only_v3.py")}
    all_rows = {}
    for n, seed in ((5, 2), (40, 1)):
        print(f"\n===== N={n} =====")
        catalog_path, catalog = orig.load_catalog(n)
        print(f"  分類表: {catalog_path.name}（族 {len(catalog)}）")
        w2 = np.load(CENSUS / f"parent_white_harmonics_N{n}_v2" / "relation_waves.npy")
        ctrl_ok = control_reproduction(n, w2, catalog)
        r3 = gen3.make_parent(n, seed=seed)
        w3 = r3.relation_waves
        v = r3.parent_vector
        single = np.outer(v, np.exp(2j * np.pi * np.arange(n) / n) / math.sqrt(n))
        rowc = float(np.max(np.abs(np.sum(single * single, axis=1))
                            / np.sum(np.abs(single) ** 2, axis=1)))
        print(f"  単一波版の行閉塞検査 max={rowc:.2e}（進行波 k=+1: 恒等0のはず）")
        res = {"control_reproduction": bool(ctrl_ok),
               "single_wave_row_closure_max": rowc}
        for label, w in (("白色v2", w2), ("白色v3", w3), ("単一波版", single)):
            agg, rows = classify_state(w, n, catalog, label)
            res[label] = agg
            all_rows[f"N{n}_{label}"] = [
                {k: r[k] for k in ("status", "reason", "base_q", "harmonic_set",
                                    "family_id", "allowed_power_frac",
                                    "outside_power_frac", "forbidden_power_frac",
                                    "max_pair_debt")} for r in rows]
        out[f"N{n}"] = res
    out["per_wave_rows"] = all_rows
    (HERE / "stable_wave_classification_fixed_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print("\nsaved: stable_wave_classification_fixed_result_v1.json")


if __name__ == "__main__":
    main()
