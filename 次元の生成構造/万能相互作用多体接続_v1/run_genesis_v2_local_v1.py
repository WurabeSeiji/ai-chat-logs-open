#!/usr/bin/env python3
"""GENESIS v2: 局所種濃縮実験——「物質＝局所凝縮点火の稀事象」の直接実現

GENESIS v1 の帰結: 一様原始種からの大域生成は rate∝δ⁴ で不可能。
    物質創成には局所濃縮チャネルが必然（二体fate・第9論文と整合）。
    本実験はそれを直接実現する: 同一総パワーの種を
    (a) 一様（全M辺に薄く） vs (b) 局所（頂点0のスター M_loc 辺に濃く）
    で置き、生成率と空間分布を比較する。

設定: N=12（M=66、スター=11辺、濃縮率 M/M_loc=6）、Nreg=5、
    ポンプ=control親（k=2）、種=厳密零閉塞状態（QR構成、支持を
    全辺/スター辺に限定）×δ、δ∈{1e-2, 3e-2}、T=2000。
    v2エンジン（共有O＋媒介頂点）。

事前登録予言（実行前固定）:
    V2-P1 局所増強: 同一総種パワーで rate_local/rate_uniform ≥ 3
        （頂点局所の媒介頂点は局所種密度の二乗で駆動されるため。
        濃縮率6に対し指数を実測: 完全二乗なら36、線形なら6）。
    V2-P2 ドメイン形成: 局所走行の奇数内容の辺参加比 PR_odd が
        成長窓の間 0.5·M 未満に留まる（物質は生まれた場所に留まる）。
    探索記録: crossing の種幾何応答（局所 vs 一様、同パワー）、
        増強指数 p_conc = ln(rate比)/ln(濃縮率)。

使い方: python3 run_genesis_v2_local_v1.py
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
    "s3g2", HERE / "run_stage3_sharedO_v2_and_hair_v1.py")
s3 = importlib.util.module_from_spec(spec3)
sys.modules[spec3.name] = s3
spec3.loader.exec_module(s3)
abl = s3.abl
gen3 = s3.gen3
V2 = s3.VertexEngineV2

N_GRAPH = 12
NREG = 5
T_LONG = 2000


def zero_closure_state(m, support_idx, rng):
    """支持を support_idx に限定した厳密零閉塞状態（QR構成）。"""
    ml = len(support_idx)
    A = rng.normal(size=(ml, 2))
    Q, _ = np.linalg.qr(A)
    z_loc = (Q[:, 0] + 1j * Q[:, 1]) / np.sqrt(2)
    z = np.zeros(m, complex)
    z[support_idx] = z_loc
    return z / np.linalg.norm(z)


def run_case(Z0c, wp0, seed_vec, delta, label):
    m = N_GRAPH * (N_GRAPH - 1) // 2
    C0 = np.zeros((m, NREG), complex)
    C0[:, 2] = Z0c
    C0[:, 1] = delta * seed_vec
    eng = V2(N_GRAPH, C0, wp0, vertex_on=True)
    p2 = C0[:, 2].real / np.linalg.norm(C0[:, 2].real)
    q2 = C0[:, 2].imag - (C0[:, 2].imag @ p2) * p2
    q2 = q2 / np.linalg.norm(q2)
    fseeds = np.zeros(T_LONG)
    prs = np.zeros(T_LONG)
    crossing = None
    for t in range(T_LONG):
        eng.step()
        d = eng.diagnostics()
        fseeds[t] = d["f_seed"]
        podd = np.abs(eng.C[:, 1]) ** 2 + np.abs(eng.C[:, 3]) ** 2
        s1 = podd.sum()
        prs[t] = (s1 ** 2 / np.sum(podd ** 2)) if s1 > 0 else 0.0
        if crossing is None:
            Z2 = eng.C[:, 2]
            Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
            if float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z2) @ Z2)) > 0.05:
                crossing = t + 1
    hi = min((crossing or T_LONG) - 100, 800)
    tt = np.arange(50, hi, dtype=float)
    lnf = np.log(fseeds[50:hi])
    A = np.vstack([tt, np.ones_like(tt)]).T
    coef, _, _, _ = np.linalg.lstsq(A, lnf, rcond=None)
    return {"label": label, "delta": delta, "crossing": crossing,
            "rate": float(coef[0]), "f_seed0": float(fseeds[0]),
            "PR_odd_init": float(prs[0]), "PR_odd_at_800": float(prs[799]),
            "PR_odd_final": float(prs[-1]), "f_seed_final": float(fseeds[-1])}


def main() -> None:
    t0 = time.time()
    m = N_GRAPH * (N_GRAPH - 1) // 2
    ia, ib = np.triu_indices(N_GRAPH, k=1)
    star = [e for e in range(m) if ia[e] == 0 or ib[e] == 0]
    conc = m / len(star)
    print(f"GENESIS v2 局所濃縮 N={N_GRAPH} M={m} スター={len(star)}辺 濃縮率={conc:.1f}")

    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(N_GRAPH, False)
    Z0c = Z0c / np.linalg.norm(Z0c)
    rng_u = np.random.default_rng(98000)
    rng_l = np.random.default_rng(98001)
    seed_uniform = zero_closure_state(m, list(range(m)), rng_u)
    seed_local = zero_closure_state(m, star, rng_l)

    out = {"N": N_GRAPH, "M": m, "star_edges": len(star), "concentration": conc,
           "criteria": {"V2_P1": "rate_local/rate_uniform >= 3 at equal power",
                         "V2_P2": "PR_odd(local) < 0.5*M through growth window"}}
    rows = []
    for delta in (1e-2, 3e-2):
        ru = run_case(Z0c, wp0, seed_uniform, delta, "uniform")
        rl = run_case(Z0c, wp0, seed_local, delta, "local")
        ratio = rl["rate"] / ru["rate"] if ru["rate"] > 0 else np.inf
        p_conc = float(np.log(ratio) / np.log(conc)) if ratio > 0 and np.isfinite(ratio) else None
        rows.append({"delta": delta, "uniform": ru, "local": rl,
                      "rate_ratio": float(ratio) if np.isfinite(ratio) else None,
                      "p_conc": p_conc})
        print(f"  δ={delta:.0e}: uniform rate={ru['rate']:.3e} (crossing={ru['crossing']}) "
              f"local rate={rl['rate']:.3e} (crossing={rl['crossing']})")
        print(f"       増強比={ratio:.2f} 指数p_conc={p_conc if p_conc else '—'} "
              f"PR_odd(local): {rl['PR_odd_init']:.1f}→{rl['PR_odd_at_800']:.1f} "
              f"(uniform: {ru['PR_odd_init']:.1f}→{ru['PR_odd_at_800']:.1f})")
    p1 = bool(all(r_["rate_ratio"] and r_["rate_ratio"] >= 3 for r_ in rows))
    p2ok = bool(all(r_["local"]["PR_odd_at_800"] < 0.5 * m for r_ in rows))
    out["rows"] = rows
    out["V2_P1_pass"] = p1
    out["V2_P2_pass"] = p2ok
    out["all_pass"] = bool(p1 and p2ok)
    out["runtime_sec"] = time.time() - t0
    print(f"\n判定: V2-P1 局所増強={p1}  V2-P2 ドメイン形成={p2ok} "
          f"→ {'ALL PASS' if out['all_pass'] else '不成立あり——反証記録'}")
    (HERE / "genesis_v2_local_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
