#!/usr/bin/env python3
"""GATE-0f: 単一平面（一様振幅）族の存在域プローブ——「最後の生き残り」説の検定

背景（2026-08-05）:
    第5論文の早期記録: (a) N相境界は4/5の間 (b) N=5で親不動点の収束異常
    （6試行中4不収束、存在問題として未解決登録）(c) 単一平面（rank2）親は
    N=4に固有らしく、大きいNの親は複数平面（σ₂/σ₁≈0.35〜0.5）。
    GATE-0c/0d/0e: 安定平衡はN=5のみ・安定解=完全一様振幅・完全縮退。
    統合仮説（最後の生き残り説）: 安定海解=N≤4型単一平面族がN=5まで生存した
    最後の一員であり、N≥6で族ごと消滅する。

予言（実行前固定・事後変更禁止）:
    P-1: 生成子K(arg v)の σ₂/σ₁ は、N=5安定解で < 1e-6（単一平面）、
         burst解（N=5および3..12）で 0.3〜0.55（記録3の範囲）。
    Q-2（判別実験・両仮説を事前登録）:
         H-survivor: 一様振幅の自己無撞着不動点は N≤5 にのみ存在
                     （無拘束固有モード残差 < 1e-8 が達成できるのは N≤5）
         H-everywhere: 全Nに存在するが、安定なのは N=5 のみ
         判別: 拘束付き不動点反復（振幅を毎回一様化）で収束候補を作り、
         無拘束の固有モード残差と接線安定性を測る。

方法: gen3.LowRankSystem の operator_matrix/sigma_spectrum/eigenmode_residual
    を read-only 利用。一様族プローブは N=3..12、位相初期値5シード
    （rng 97000+j）、反復2000回。接線安定性はONS-2と同じ中心差分。

使い方: python3 run_stage0f_uniform_family_probe_v1.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SERIES = HERE.parent
ABL = SERIES / "第8論文_二段階seed除去による準安定相の因果分離" / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = SERIES / "make_parent_white_managed_v1" / "make_parent_white_harmonics_n_only_v3.py"
H_FD = 1e-7
ITERS = 2000
PHASE_SEEDS = 5

spec = importlib.util.spec_from_file_location("abl_f", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen_f", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sigma_ratio(n, v):
    s = gen3.LowRankSystem(n)
    s.set_theta(np.angle(v))
    sig = np.sort(np.abs(s.sigma_spectrum()))[::-1]
    sig = sig[sig > 1e-14 * sig[0]] if sig[0] > 0 else sig
    return float(sig[1] / sig[0]) if len(sig) > 1 else 0.0


def build_K(n, theta):
    """K_ee' = sin(θ_{e'}−θ_e)（頂点共有辺対のみ・実反対称）。第8論文§1の定義。"""
    ia, ib = np.triu_indices(n, k=1)
    m = len(ia)
    K = np.zeros((m, m))
    for e in range(m):
        for f in range(m):
            if e == f:
                continue
            if ia[e] in (ia[f], ib[f]) or ib[e] in (ia[f], ib[f]):
                K[e, f] = np.sin(theta[f] - theta[e])
    return K


def eigvec_max(K):
    w, V = np.linalg.eig(K)
    idx = int(np.argmin(w.imag))          # 最小虚部 = −iσ_max
    v = V[:, idx]
    return v / np.linalg.norm(v)


def uniform_probe(n, phase_seed):
    m = n * (n - 1) // 2
    rng = np.random.default_rng(97000 + phase_seed)
    v = np.exp(1j * rng.uniform(0, 2 * np.pi, m)) / np.sqrt(m)
    last = None
    for it in range(ITERS):
        K = build_K(n, np.angle(v))
        u = eigvec_max(K)
        ip = np.vdot(v, u)
        if abs(ip) > 0:
            u = u * (np.conj(ip) / abs(ip))
        u = np.exp(1j * np.angle(u)) / np.sqrt(m)   # 一様振幅拘束
        step = float(np.linalg.norm(u - v))
        v = 0.5 * v + 0.5 * u
        v = np.exp(1j * np.angle(v)) / np.sqrt(m)
        last = step
        if step < 1e-13:
            break
    K = build_K(n, np.angle(v))
    Kv = K @ v
    lam = np.vdot(v, Kv)
    resid = float(np.linalg.norm(Kv - lam * v))
    return v, float(last), resid


def make_step(n, wp_base):
    sys_lr = abl.LowRankSystem(n)

    def F(Z):
        wp = wp_base.copy()
        sys_lr.set_theta(np.angle(Z))
        se, _ = sys_lr.sigma_max_power(wp)
        return sys_lr.cayley_step(Z.copy(), se)

    return F


def lambda_max(n, v):
    m = n * (n - 1) // 2
    wp = np.random.default_rng(91000).normal(size=m)
    F = make_step(n, wp)
    G0 = F(v)
    ip = np.conj(v) @ G0
    ph = ip / abs(ip)
    cols = [np.eye(m)[i] + 0j for i in range(m)] + [1j * np.eye(m)[i] for i in range(m)]
    J = np.zeros((2 * m, 2 * m))
    for idx, d in enumerate(cols):
        Gp = np.conj(ph) * F(v + H_FD * d)
        Gm = np.conj(ph) * F(v - H_FD * d)
        g = (Gp - Gm) / (2 * H_FD)
        J[:, idx] = np.concatenate([g.real, g.imag])
    ev = np.linalg.eigvals(J)
    return float(np.log(np.abs(ev).max()))


def main() -> None:
    t0 = time.time()
    d0d = json.load(open(HERE / "gate0d_parent_seed_sweep_v1.json"))
    out = {"P1": {}, "Q2": {}}

    print("=== P-1: 生成子 σ₂/σ₁ ===")
    for seed, kind in [(2, "stable"), (12, "stable"), (20, "stable"),
                        (7, "burst"), (10, "burst")]:
        r = gen3.make_parent(5, seed=seed)
        v = r.parent_vector / np.linalg.norm(r.parent_vector)
        sr = sigma_ratio(5, v)
        s5 = gen3.LowRankSystem(5)
        s5.set_theta(np.angle(v))
        spec_full = [float(x) for x in s5.sigma_spectrum()]
        rel = [round(x / spec_full[0], 6) for x in spec_full] if spec_full else []
        out["P1"][f"N5_seed{seed}"] = {"kind": kind, "sigma_ratio": sr,
                                        "sigma_rel": rel}
        print(f"  N=5 seed={seed:2d} [{kind}]: σ₂/σ₁ = {sr:.6f}  σ/σ₁ = {rel}")
    for n in range(3, 13):
        rows = d0d[str(n)]
        seed = next((s for s, c, _ in rows if c == "burst"), None)
        if seed is None:
            continue
        r = gen3.make_parent(n, seed=seed)
        v = r.parent_vector / np.linalg.norm(r.parent_vector)
        sr = sigma_ratio(n, v)
        out["P1"][f"N{n}_burst"] = {"kind": "burst", "sigma_ratio": sr}
        print(f"  N={n:2d} burst(seed={seed}): σ₂/σ₁ = {sr:.6f}")

    print("=== Q-2: 一様族の存在域 ===")
    for n in range(3, 13):
        best = None
        for j in range(PHASE_SEEDS):
            v, step, resid = uniform_probe(n, j)
            if best is None or resid < best[2]:
                best = (v, step, resid, j)
        v, step, resid, j = best
        exists = resid < 1e-8
        lam = lambda_max(n, v) if exists else None
        out["Q2"][str(n)] = {"best_seed": j, "constrained_step": step,
                              "eigenmode_residual": resid, "exists": bool(exists),
                              "lambda_max": lam,
                              "sigma_ratio": sigma_ratio(n, v) if exists else None}
        stat = f"存在 λ_max={lam:.5f} σ₂/σ₁={out['Q2'][str(n)]['sigma_ratio']:.4f}" if exists else "不存在"
        print(f"  N={n:2d}: 残差={resid:.2e} → {stat}")

    exist_ns = [int(n) for n, d in out["Q2"].items() if d["exists"]]
    print(f"\n一様族が存在するN: {exist_ns}")
    verdict = {"P1_stable_single_plane": all(
                   d["sigma_ratio"] < 1e-6 for k, d in out["P1"].items() if d["kind"] == "stable"),
               "P1_burst_range": all(
                   0.3 <= d["sigma_ratio"] <= 0.55 for k, d in out["P1"].items() if d["kind"] == "burst"),
               "Q2_exist_ns": exist_ns}
    print(f"判定: P-1 stable単一平面={verdict['P1_stable_single_plane']} "
          f"burst範囲={verdict['P1_burst_range']}")

    out["verdict"] = verdict
    out["criteria"] = {"P1": "stable σ2/σ1<1e-6, burst in [0.3,0.55]",
                        "Q2": "uniform fixed point exists iff residual<1e-8; H-survivor: N<=5 only"}
    out["runtime_sec"] = time.time() - t0
    (HERE / "gate0f_uniform_family_probe_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
