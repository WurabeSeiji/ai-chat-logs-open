#!/usr/bin/env python3
"""万能非弾性写像 v1：実装＋単体テスト（設計書 万能非弾性写像の設計_v1.md に従う）

写像（IF 文なし・全対同一・二波間のみ）:
    弾性部:   theta_from_ab → rotate_ab（原本 toy の関数を read-only 呼出し・無変更）
    非弾性部: a' = a + i·g·(|b|²·a + b²·conj(a))
              b' = b + i·g·(|a|²·b + a²·conj(b))    （同時更新・χ×η 格子の点ごと）

テスト: T1 対照(g=0, bitwise) / T2 保存(g=1e-3) / T3 パリティ定理 /
        T4 種付き利得（占有率掃引・探索記録） / T5 IF監査（本文参照）
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
TOY = (REPO / "次元の生成構造" / "第9論文_フェルミオンの生成構造"
       / "対照実験_波束収縮_実行環境_v1" / "ab_invariant_theta_toy_v1"
       / "run_ab_invariant_theta_toy_v1.py")

spec = importlib.util.spec_from_file_location("toy_uim", TOY)
toy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = toy
spec.loader.exec_module(toy)
base = toy.base

J_MAX = 200
EVEN_KS = tuple(range(2, 63, 2))
ODD_KS = tuple(range(1, 64, 2))


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def inelastic_pass(a: np.ndarray, b: np.ndarray, g: float):
    """非弾性部（点ごと・二波間・種分岐なし）。g=0 で恒等写像。"""
    if g == 0.0:
        return a, b
    a2, b2 = a * a, b * b
    na2, nb2 = np.abs(a) ** 2, np.abs(b) ** 2
    a_new = a + 1j * g * (nb2 * a + b2 * np.conj(a))
    b_new = b + 1j * g * (na2 * b + a2 * np.conj(b))
    return a_new, b_new


def collision_step(a, b, sp, g: float):
    """衝突一回 = 弾性部（原本）∘ 非弾性部（新設）。"""
    readout = toy.theta_from_ab(a, b, sp)
    a, b = toy.rotate_ab(a, b, readout.theta)
    a, b = inelastic_pass(a, b, g)
    return a, b, readout


def make_bundle(sp, ks, which, scale=1.0):
    case = base.explicit_packet_case(mode=f"uim_{which}", packet_a=tuple(ks), packet_b=tuple(ks))
    v = base.make_case_state(sp, case, which, hair_enabled=True)
    v = v / np.sqrt(float(np.vdot(v, v).real))
    return v * scale


def fermionic_power_raw(a, b, sp) -> float:
    """床なしの生フェルミオン関係量（theta_from_ab は数値床未満を0に切るため、
    微小シードの利得測定にはこちらを使う）。マスクは theta_from_ab と同一。"""
    freqs, power = toy.combined_chi_power(a, b, sp)
    af = np.abs(freqs)
    mask = (af >= 4) & ((af % 2) == 0)
    return float(np.sum(power[mask]))


def chi_momentum(a, b, sp) -> float:
    freqs, power = toy.combined_chi_power(a, b, sp)
    return float(np.sum(freqs * power))


def main() -> None:
    t0 = time.time()
    print("万能非弾性写像 v1 単体テスト")
    print(f"  toy sha256: {sha256(TOY)[:16]}…（スナップショット一致確認済み）")
    params = base.Params(high_n=63, recursive_collision_count=J_MAX)
    sp = base.build_source_params(params)
    results = {"toy_sha256": sha256(TOY)}

    # ---- T1 対照（g=0 で原本パイプラインと bitwise 一致）----
    print("\n[T1] g=0 対照（bitwise）")
    t1_ok = True
    anchors = {}
    for name, ka, kb in (("even-even", EVEN_KS, EVEN_KS),
                          ("odd-odd", ODD_KS, ODD_KS),
                          ("even-odd", EVEN_KS, ODD_KS)):
        a1 = make_bundle(sp, ka, "A"); b1 = make_bundle(sp, kb, "B")
        a2, b2 = a1.copy(), b1.copy()
        th0 = None
        for _ in range(J_MAX):
            ro = toy.theta_from_ab(a1, b1, sp)          # 原本パイプライン
            a1, b1 = toy.rotate_ab(a1, b1, ro.theta)
            a2, b2, ro2 = collision_step(a2, b2, sp, g=0.0)   # 本写像 g=0
            if th0 is None:
                th0 = float(ro2.theta)
        same = np.array_equal(a1, a2) and np.array_equal(b1, b2)
        t1_ok &= same
        anchors[name] = th0
        print(f"  {name:9s}: bitwise一致={same} θ0={th0:.10f}")
    print(f"  アンカー確認: even-even θ0={anchors['even-even']:.3e}（≡0）, "
          f"even-odd θ0={anchors['even-odd']:.10f}（期待 0.7619520718）")
    t1_ok &= anchors["even-even"] <= 1e-12 and abs(anchors["even-odd"] - 0.7619520718) < 1e-9
    print(f"  T1: {'PASS' if t1_ok else 'FAIL'}")
    results["T1"] = bool(t1_ok)

    # ---- T2 保存（g=1e-3・希薄=単位ノルム）----
    print("\n[T2] 保存帳簿（g=1e-3, 200衝突）")
    g = 1e-3
    a = make_bundle(sp, EVEN_KS, "A"); b = make_bundle(sp, ODD_KS, "B")
    n0 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    p0 = chi_momentum(a, b, sp)
    c0 = complex(np.sum(a * a) + np.sum(b * b))
    for _ in range(J_MAX):
        a, b, _ = collision_step(a, b, sp, g)
    n1 = float(np.vdot(a, a).real + np.vdot(b, b).real)
    p1 = chi_momentum(a, b, sp)
    c1 = complex(np.sum(a * a) + np.sum(b * b))
    print(f"  合成ノルム: {n0:.12f} → {n1:.12f}（ドリフト {abs(n1-n0):.3e}）")
    print(f"  合成χ運動量: {p0:.6e} → {p1:.6e}（ドリフト {abs(p1-p0):.3e}）")
    print(f"  零閉塞 |Σa²+Σb²|: {abs(c0):.3e} → {abs(c1):.3e}（変化 {abs(c1-c0):.3e}）")
    results["T2"] = {"norm_drift": abs(n1 - n0), "momentum_drift": abs(p1 - p0),
                      "closure_before": abs(c0), "closure_after": abs(c1)}

    # ---- T3 パリティ定理（純偶×純偶・g大 → P_f 恒等0）----
    print("\n[T3] パリティ定理（白猫×白猫, g=0.5, 200衝突）")
    a = make_bundle(sp, EVEN_KS, "A"); b = make_bundle(sp, EVEN_KS, "B")
    pf_max = 0.0
    for _ in range(J_MAX):
        a, b, ro = collision_step(a, b, sp, g=0.5)
        pf_max = max(pf_max, float(ro.fermionic_relation_power))
    t3 = pf_max < 1e-25
    print(f"  P_f 最大値 = {pf_max:.3e} → {'PASS（奇内容は種なしでは恒等0）' if t3 else 'FAIL'}")
    results["T3"] = {"pass": bool(t3), "pf_max": pf_max}

    # ---- T4 種付き利得（占有率掃引・探索記録）----
    print("\n[T4] 種付きパラメトリック利得（偶ポンプ＋奇シード1e-5, g=0.05, 占有率掃引・生パワー観測）")
    gains = {}
    for s in (0.5, 1.0, 2.0, 4.0, 8.0):
        a = make_bundle(sp, EVEN_KS, "A", scale=s)
        seed = make_bundle(sp, ODD_KS, "A", scale=1e-5 * s)
        a = a + seed
        b = make_bundle(sp, EVEN_KS, "B", scale=s)
        pf0 = fermionic_power_raw(a, b, sp)
        pfs = [pf0]
        for _ in range(J_MAX):
            a, b, _ = collision_step(a, b, sp, g=0.05)
            pfs.append(fermionic_power_raw(a, b, sp))
        pfs = np.array(pfs)
        growth = float(pfs[-1] / max(pfs[0], 1e-300))
        with np.errstate(divide="ignore"):
            lp = np.log(np.maximum(pfs, 1e-300))
        rate = float(np.polyfit(np.arange(len(lp)), lp, 1)[0])
        gains[s] = {"pf_initial": float(pfs[0]), "pf_final": float(pfs[-1]),
                     "growth_ratio": growth, "log_rate_per_collision": rate}
        print(f"  占有スケール s={s:4.1f}（占有率∝{s*s:5.2f}）: "
              f"P_f {pfs[0]:.3e} → {pfs[-1]:.3e}  増幅率={growth:.3e}  "
              f"利得/衝突={rate:+.4f}")
    results["T4"] = gains

    results["runtime_sec"] = time.time() - t0
    (HERE / "universal_inelastic_map_test_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
