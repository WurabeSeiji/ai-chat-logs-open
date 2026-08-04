#!/usr/bin/env python3
"""ONS-2: 一段残差と接線写像スペクトル——「増幅⟺不安定な相対平衡」の直接検定

目的（2026-08-05 木原氏合意・開始様式判別論文の機構検証その2）:
    (1) 一段残差 r(Z0)=min_φ‖F(Z0)−e^{iφ}Z0‖ を全系列で測り、
        「自己無撞着閉包＝相対平衡」を直接観測量にする。
    (2) 4つの親状態（control N=5/40、白色v3親 N=5/40）で位相整合写像
        G̃(Z)=e^{-iφ0}F(Z) の実 2M 次元ヤコビアンを中心差分で構成し、
        固有値の最大絶対値から振幅成長率 λ_max=ln max|μ| を測る。

予言（実行前固定・事後変更禁止）:
    P1 残差: 4親は r ≤ 1e-10、白色セクター42本は r=O(1e-2〜1)。中間なし。
    P2 成長率: f は振幅二乗量なので λ_max = rate_f/2。
        control N=5: λ_max ≈ 0.04935/2 = 0.02468
        control N=40: λ_max ≈ 0.03499/2 = 0.01750
        白色v3親 N=40: λ_max ≈ 0.03468/2 = 0.01734
        白色v3親 N=5（不活性）: λ_max ≤ 1e-3（不安定固有値なし）
    P2 が通れば、必要条件は「増幅⟺不安定な相対平衡」の同値条件に昇格し、
    外れれば §5 の機構的読み（線形無生成仮説の解釈部）が反証される。

方法:
    F は第7論文力学の一段（evolve。wp は基点の値を毎評価 copy して固定）。
    ヤコビアン: h=1e-7 の中心差分、実基底 2M 本。位相整合 φ0 は基点で固定。
    固有値は実 2M×2M 行列の固有値（複素対）。

再現性: abl・生成器v3 とも read-only import（SHA-256 記録）。生成器シードは
    確定値（N=5:2 / N=40:1）。全結果 JSON 保存、固有値分布の図出力。

使い方: python3 run_onset_equilibrium_residual_v1.py <N>
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = (REPO / "次元の生成構造" / "make_parent_white_managed_v1"
        / "make_parent_white_harmonics_n_only_v3.py")
SEEDS = {5: 2, 40: 1}
H_FD = 1e-7
RATE_REF = {"control_5": 0.04935428519548658, "control_40": 0.034990401406084726,
            "white_parent_40": 0.03468168424655808}

spec = importlib.util.spec_from_file_location("abl_ons2", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen3_ons2", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def make_step(n, wp_base):
    sys_lr = abl.LowRankSystem(n)

    def F(Z):
        Zc = Z.copy()
        wp = wp_base.copy()
        sys_lr.set_theta(np.angle(Zc))
        se, _ = sys_lr.sigma_max_power(wp)
        return sys_lr.cayley_step(Zc, se)

    return F


def residual(F, Z0):
    G = F(Z0)
    ip = np.conj(Z0) @ G
    phase = ip / abs(ip) if abs(ip) > 0 else 1.0
    r = np.linalg.norm(G - phase * Z0) / np.linalg.norm(Z0)
    return float(r), complex(phase)


def jacobian_spectrum(F, Z0, phase0):
    m = Z0.shape[0]
    dim = 2 * m

    def to_real(Z):
        return np.concatenate([Z.real, Z.imag])

    J = np.zeros((dim, dim))
    basis = np.eye(m)
    cols = []
    for i in range(m):
        cols.append(basis[i] + 0j)
    for i in range(m):
        cols.append(1j * basis[i])
    for idx, d in enumerate(cols):
        Gp = np.conj(phase0) * F(Z0 + H_FD * d)
        Gm = np.conj(phase0) * F(Z0 - H_FD * d)
        J[:, idx] = to_real((Gp - Gm) / (2 * H_FD))
    ev = np.linalg.eigvals(J)
    mods = np.abs(ev)
    order = np.argsort(mods)[::-1]
    top = [{"abs": float(mods[i]), "ln_abs": float(np.log(mods[i])),
            "re": float(ev[i].real), "im": float(ev[i].imag)} for i in order[:10]]
    return {"dim": dim, "max_abs": float(mods.max()),
            "lambda_max": float(np.log(mods.max())),
            "n_unstable_gt_1e-3": int(np.sum(np.log(mods) > 1e-3)),
            "top10": top}, mods


def main() -> None:
    t0 = time.time()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    seed = SEEDS[n]
    print(f"ONS-2 平衡残差＋接線スペクトル N={n}（v3 seed={seed}）")
    print(f"  import: ABL {sha256(ABL)[:16]}…  GEN3 {sha256(GEN3)[:16]}…")

    _, v, _, _, _, p, q, Z0c, wp0 = abl.build_init(n, False)
    m = v.shape[0]
    F = make_step(n, wp0)

    r = gen3.make_parent(n, seed=seed)
    W = r.relation_waves
    C = np.fft.fft(W, axis=1) / n
    vp = r.parent_vector / np.linalg.norm(r.parent_vector)

    results = {"residuals": {}, "jacobians": {}}

    r0, ph0 = residual(F, v / np.linalg.norm(v))
    results["residuals"]["control_parent"] = r0
    rp, php = residual(F, vp)
    results["residuals"]["white_v3_parent"] = rp
    print(f"  r(control 親)={r0:.3e}  r(白色v3親)={rp:.3e}")

    sec_res = {}
    for k in range(n):
        if (2 * k) % n == 0:
            continue
        Zk = C[:, k] / np.linalg.norm(C[:, k])
        rk, _ = residual(F, Zk)
        sec_res[f"k{k}"] = rk
    results["residuals"]["sectors"] = sec_res
    vals = list(sec_res.values())
    print(f"  セクター{len(vals)}本: r ∈ [{min(vals):.3e}, {max(vals):.3e}]")

    spec_c, mods_c = jacobian_spectrum(F, v / np.linalg.norm(v), ph0)
    results["jacobians"]["control_parent"] = spec_c
    print(f"  control 親: λ_max={spec_c['lambda_max']:.5f} "
          f"(予言 {RATE_REF[f'control_{n}']/2:.5f})")
    spec_p, mods_p = jacobian_spectrum(F, vp, php)
    results["jacobians"]["white_v3_parent"] = spec_p
    pred_p = RATE_REF.get(f"white_parent_{n}")
    print(f"  白色v3親: λ_max={spec_p['lambda_max']:.5f} "
          f"(予言 {'%.5f' % (pred_p/2) if pred_p else '≤1e-3（不活性）'})")

    out = {"N": n, "seed": seed,
           "imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)},
           "criteria": {"H_FD": H_FD,
                         "prediction": {
                             "P1": "parents r<=1e-10, sectors r=O(1e-2..1)",
                             "P2_lambda_max": {k: rv / 2 for k, rv in RATE_REF.items()},
                             "P2_white_parent_5": "lambda_max <= 1e-3"}},
           **results, "runtime_sec": time.time() - t0}
    (HERE / f"onset_equilibrium_residual_N{n:05d}_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(7, 5))
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(th), np.sin(th), "k:", lw=0.5)
    evc = np.array([[e["re"], e["im"]] for e in spec_c["top10"]])
    evp = np.array([[e["re"], e["im"]] for e in spec_p["top10"]])
    ax.plot(evc[:, 0], evc[:, 1], "o", ms=6, color="#7F7F7F", label="control parent top10")
    ax.plot(evp[:, 0], evp[:, 1], "s", ms=5, color="black", label="white v3 parent top10")
    ax.set_aspect("equal")
    ax.set_xlabel("Re μ")
    ax.set_ylabel("Im μ")
    ax.set_title(f"N={n} ONS-2: tangent-map eigenvalues (top 10 by |μ|)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_onset_tangent_spectrum_N{n:05d}.png", dpi=130)
    plt.close(fig)
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
