#!/usr/bin/env python3
"""GENESIS v3: レジスタ空間の局所点火——局所性の正しい軸での稀事象実現

v2の反証の理論的整理（実行前に固定する設計論拠）:
    K_N 完全グラフは直径1でグラフ的局所性を持たない（スターのハブは
    全種を集約する——v2の濃縮盲目の構造的理由）。一方、媒介頂点は
    レジスタ点 n ごとの点ごと積であり、**この理論の局所性は最初から
    レジスタ軸（関係波形の位置軸）に住んでいる**。局所凝縮点火
    （物質=稀事象）の正しい実装は、レジスタ空間の波束濃縮である。

設定: N=12（M=66）、Nreg=16（許容k: 2k≢0 mod 16、奇数k=8本）。
    ポンプ = control親 × 単一倍音 k=2（レジスタ上一様強度）。
    種（同一総パワー δ=3e-2、同一辺プロファイル=零閉塞状態）:
    (a) 非局在: 単一奇数倍音 k=1（|s(n)|²=一様、参加比 PR_n=16）
    (b) 局在: 奇数倍音全8本の等振幅同位相パケット（n=0に局在、PR_n≈2-3）
    辺空間の零閉塞は積構造により全レジスタ点で厳密（s(n)²·Σz²=0）。

事前登録予言（実行前固定）:
    V3-P1 レジスタ局所増強: rate_packet / rate_single ≥ 2。
        濃縮率 = PR_n(single)/PR_n(packet) に対する指数 p_reg を実測
        （点ごと三次頂点は峰強度で駆動されるため p_reg > 0 を予言。
        v2のグラフ濃縮では p≈0.1 だった——軸の違いが本質かの判別）。
    V3-P2 レジスタ・ドメイン形成: パケット走行の生成後の奇数内容の
        レジスタ参加比 PR_n が 0.5·Nreg 未満に留まる（物質は生誕位置に
        留まる——レジスタ上のドメイン）。
    探索記録: crossing の種形状応答、ドメイン位置の追跡（n=0近傍か）。

使い方: python3 run_genesis_v3_register_local_v1.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec3 = importlib.util.spec_from_file_location(
    "s3g3", HERE / "run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = s3
spec3.loader.exec_module(s3)
abl = s3.abl
V2 = s3.VertexEngineV2

N_GRAPH = 12
NREG = 16
T_LONG = 800
DELTA = 3e-2


def zero_closure_state(m, rng):
    A = rng.normal(size=(m, 2))
    Q, _ = np.linalg.qr(A)
    z = (Q[:, 0] + 1j * Q[:, 1]) / np.sqrt(2)
    return z / np.linalg.norm(z)


def pr_n(x):
    s = x.sum()
    return float(s ** 2 / np.sum(x ** 2)) if s > 0 else 0.0


def odd_content_by_n(C):
    ks = np.arange(C.shape[1])
    Codd = C * ((ks % 2 == 1)[None, :])
    Wodd = np.fft.ifft(Codd, axis=1) * C.shape[1]
    return np.sum(np.abs(Wodd) ** 2, axis=0)          # P(n)


def run_case(Z0c, wp0, seed_edge, reg_profile, label):
    m = N_GRAPH * (N_GRAPH - 1) // 2
    C0 = np.zeros((m, NREG), complex)
    C0[:, 2] = Z0c
    for k in range(NREG):
        if abs(reg_profile[k]) > 0:
            C0[:, k] += DELTA * reg_profile[k] * seed_edge
    eng = V2(N_GRAPH, C0, wp0, vertex_on=True)
    p2 = Z0c.real / np.linalg.norm(Z0c.real)
    q2 = Z0c.imag - (Z0c.imag @ p2) * p2
    q2 = q2 / np.linalg.norm(q2)
    Pn0 = odd_content_by_n(C0)
    pr0 = pr_n(Pn0)
    fseeds = np.zeros(T_LONG)
    crossing = None
    for t in range(T_LONG):
        eng.step()
        fseeds[t] = eng.diagnostics()["f_seed"]
        if crossing is None:
            Z2 = eng.C[:, 2]
            Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
            if float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z2) @ Z2)) > 0.05:
                crossing = t + 1
    Pn1 = odd_content_by_n(eng.C)
    hi = min((crossing or T_LONG) - 50, 400)
    tt = np.arange(20, hi, dtype=float)
    A = np.vstack([tt, np.ones_like(tt)]).T
    coef, _, _, _ = np.linalg.lstsq(A, np.log(fseeds[20:hi]), rcond=None)
    return {"label": label, "crossing": crossing, "rate": float(coef[0]),
            "f_seed0": float(fseeds[0]), "f_seed_final": float(fseeds[-1]),
            "PR_n_init": pr0, "PR_n_final": pr_n(Pn1),
            "Pn_final_argmax": int(np.argmax(Pn1))}


def main() -> None:
    t0 = time.time()
    m = N_GRAPH * (N_GRAPH - 1) // 2
    print(f"GENESIS v3 レジスタ局所点火 N={N_GRAPH} M={m} Nreg={NREG} δ={DELTA}")
    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    seed_edge = zero_closure_state(m, np.random.default_rng(98000))

    odd_ks = [k for k in range(NREG) if k % 2 == 1]
    prof_single = np.zeros(NREG, complex)
    prof_single[1] = 1.0
    prof_packet = np.zeros(NREG, complex)
    for k in odd_ks:
        prof_packet[k] = 1.0 / np.sqrt(len(odd_ks))     # 等振幅同位相 → n=0 局在

    ru = run_case(Z0c, wp0, seed_edge, prof_single, "single_k1")
    rl = run_case(Z0c, wp0, seed_edge, prof_packet, "odd_packet")
    conc = ru["PR_n_init"] / rl["PR_n_init"]
    ratio = rl["rate"] / ru["rate"] if ru["rate"] > 0 else np.inf
    p_reg = float(np.log(ratio) / np.log(conc)) if (np.isfinite(ratio) and ratio > 0 and conc > 1) else None
    print(f"  single k=1 : PR_n={ru['PR_n_init']:.1f} rate={ru['rate']:.3e} "
          f"crossing={ru['crossing']} f_seed {ru['f_seed0']:.1e}→{ru['f_seed_final']:.1e}")
    print(f"  odd packet : PR_n={rl['PR_n_init']:.1f} rate={rl['rate']:.3e} "
          f"crossing={rl['crossing']} f_seed {rl['f_seed0']:.1e}→{rl['f_seed_final']:.1e}")
    print(f"  濃縮率={conc:.2f} 増強比={ratio:.2f} 指数 p_reg={p_reg}")
    print(f"  V3-P2: パケット走行 PR_n {rl['PR_n_init']:.1f}→{rl['PR_n_final']:.1f} "
          f"（閾 {0.5*NREG:.0f}） ドメイン位置 argmax n={rl['Pn_final_argmax']}")

    p1 = bool(np.isfinite(ratio) and ratio >= 2.0)
    p2ok = bool(rl["PR_n_final"] < 0.5 * NREG)
    out = {"N": N_GRAPH, "NREG": NREG, "DELTA": DELTA,
           "criteria": {"V3_P1": "rate_packet/rate_single >= 2",
                         "V3_P2": "PR_n(packet run, final) < 0.5*NREG"},
           "single": ru, "packet": rl,
           "concentration": float(conc), "rate_ratio": float(ratio) if np.isfinite(ratio) else None,
           "p_reg": p_reg, "V3_P1_pass": p1, "V3_P2_pass": p2ok,
           "all_pass": bool(p1 and p2ok), "runtime_sec": time.time() - t0}
    print(f"\n判定: V3-P1={p1} V3-P2={p2ok} → "
          f"{'ALL PASS——レジスタ軸が局所性の正しい住所' if out['all_pass'] else '不成立あり——反証記録'}")
    (HERE / "genesis_v3_register_local_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
