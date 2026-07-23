#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N3 統計の梯子 測定拡張（対照テスト済み phase5 の関数を import）

【再現性の規約】
  独自再実装をしない。control_test_v1.py で公開結果の再現を確認した
  phase5_eigenphase_resonance_v2.py の解析関数（theta_from_R / scattering_matrix /
  exact_eigenvalues / resonant_R / matrix_return_error）をそのまま import して用いる。

【測定内容（N3 の主張 = 統計の梯子と有限位数根）】
  2体交換散乱行列 U(R) の固有値構造を R∈[0,1] で調べる:
    (L1) 対称チャネル固有値 λ_sym は R によらず恒等的に +1（共有＝ボゾン的チャネルは常に存在）。
    (L2) 反対称チャネル固有値 λ_a = −e^{2iθ}, θ=asin(√R) は常に単位円上（|λ_a|=1、中性）。
    (L3) 端点: R=1 → λ_a=+1（交換なし・ボゾン的自明共有）、R=0 → λ_a=−1（完全交換・フェルミオン的）。
    (L4) 中間の有限位数根 R=cos²(πm/n)（n≥3）では λ_a が1のn乗根＝有理位相（エニオン的）で、
         U^n=I が機械精度で成立（|U^n−I|≈0）。
  すなわち統計性は R=1（ボゾン）→ 有理根（エニオン）→ R=0（フェルミオン）の梯子として並ぶ。

【スコープ】根構造と統計のみ。α⁻¹同定は扱わない。

実行: python3 statistics_ladder_extension_v1.py
出力: statistics_ladder_result_v1.json
"""

import importlib.util
import json
import math
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "finite_order_resonance_v1", "src",
                   "phase5_eigenphase_resonance_v2.py")
spec = importlib.util.spec_from_file_location("phase5_verified", SRC)
p5 = importlib.util.module_from_spec(spec)
sys.modules["phase5_verified"] = p5  # dataclass 解決のため登録
spec.loader.exec_module(p5)


def eigenphase_cycles(lam):
    """λ の偏角を 1周 (=2π) 単位で返す。"""
    return (math.atan2(lam.imag, lam.real) / (2.0 * math.pi)) % 1.0


def main():
    # --- L1/L2: R 全域スイープで λ_sym=1, |λ_a|=1 を確認 ---
    Rs = np.linspace(0.0, 1.0, 2001)
    max_sym_dev = 0.0
    max_mod_dev = 0.0
    for R in Rs:
        lam_s, lam_a = p5.exact_eigenvalues(float(R))
        max_sym_dev = max(max_sym_dev, abs(lam_s - 1.0))
        max_mod_dev = max(max_mod_dev, abs(abs(lam_a) - 1.0))
    print("=== L1/L2: R∈[0,1] 全域（2001点） ===")
    print(f"  対称固有値 max|λ_sym−1|      = {max_sym_dev:.2e}（0＝ボゾン的チャネル恒存）")
    print(f"  反対称固有値 max||λ_a|−1|    = {max_mod_dev:.2e}（0＝中性・単位円上）")

    # --- L3/L4: 統計の梯子（端点＋低位数根） ---
    print("\n=== L3/L4: 統計の梯子 R=cos²(πm/n) ===")
    print(f"{'統計':>10} {'n':>3} {'m':>3} {'R':>12} {'λ_a位相/2π':>11} "
          f"{'|U^n-I|':>10}")
    ladder = []
    # 端点ボゾン R=1（n=1,m=0 相当）
    lam_s, lam_a = p5.exact_eigenvalues(1.0)
    row = {"stat": "ボゾン端", "n": 1, "m": 0, "R": 1.0,
           "lambda_a_re": lam_a.real, "lambda_a_im": lam_a.imag,
           "phase_cycles": eigenphase_cycles(lam_a), "return_err": None}
    ladder.append(row)
    print(f"{'ボゾン端':>10} {1:>3} {0:>3} {1.0:>12.6f} "
          f"{eigenphase_cycles(lam_a):>11.5f} {'—':>10}")
    # 中間エニオン根 n=3..6
    for n in range(3, 7):
        for m in range(1, (n) // 2 + 1):
            if math.gcd(m, n) != 1:
                continue
            R = p5.resonant_R(n, m) if m < n / 2 else math.cos(math.pi * m / n) ** 2
            lam_s, lam_a = p5.exact_eigenvalues(R)
            ret = p5.matrix_return_error(R, n)
            ph = eigenphase_cycles(lam_a)
            stat = "エニオン"
            ladder.append({"stat": stat, "n": n, "m": m, "R": R,
                           "lambda_a_re": lam_a.real, "lambda_a_im": lam_a.imag,
                           "phase_cycles": ph, "return_err": ret})
            print(f"{stat:>10} {n:>3} {m:>3} {R:>12.6f} {ph:>11.5f} {ret:>10.2e}")
    # 端点フェルミオン R=0（n=2,m=1 の境界）
    lam_s, lam_a = p5.exact_eigenvalues(0.0)
    ret2 = p5.matrix_return_error(0.0, 2)
    ladder.append({"stat": "フェルミオン端", "n": 2, "m": 1, "R": 0.0,
                   "lambda_a_re": lam_a.real, "lambda_a_im": lam_a.imag,
                   "phase_cycles": eigenphase_cycles(lam_a), "return_err": ret2})
    print(f"{'フェルミオン端':>10} {2:>3} {1:>3} {0.0:>12.6f} "
          f"{eigenphase_cycles(lam_a):>11.5f} {ret2:>10.2e}")

    # 検定: 端点の固有値が厳密に ±1 か、根で U^n=I か
    _, lam_boson = p5.exact_eigenvalues(1.0)
    _, lam_fermi = p5.exact_eigenvalues(0.0)
    boson_err = abs(lam_boson - 1.0)
    fermi_err = abs(lam_fermi + 1.0)
    max_root_return = max(r["return_err"] for r in ladder if r["return_err"] is not None)
    print(f"\nボゾン端 |λ_a−1|      = {boson_err:.2e}（0 が予言）")
    print(f"フェルミオン端 |λ_a+1| = {fermi_err:.2e}（0 が予言）")
    print(f"有限位数根 max|U^n−I|  = {max_root_return:.2e}（機械精度で0）")

    ok = (max_sym_dev < 1e-14 and max_mod_dev < 1e-12
          and boson_err < 1e-14 and fermi_err < 1e-14 and max_root_return < 1e-9)
    print("判定:", "PASS（統計の梯子が機械精度で成立）" if ok else "FAIL")

    with open(os.path.join(HERE, "statistics_ladder_result_v1.json"), "w") as f:
        json.dump({"max_sym_dev": max_sym_dev, "max_mod_dev": max_mod_dev,
                   "boson_endpoint_err": boson_err, "fermion_endpoint_err": fermi_err,
                   "max_root_return_err": max_root_return, "ladder": ladder},
                  f, indent=2, ensure_ascii=False)

    make_figures(ladder, os.path.join(HERE, "figures"))
    return 0 if ok else 1


def make_figures(ladder, outdir):
    """図1（単位円上の固有値軌跡）と図2（R vs 交換位相の梯子）を SVG+PNG で生成。
    図内ラベルは英数字（フォント非依存）。キャプションは論文本文（日本語）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(outdir, exist_ok=True)

    # --- 図1: 複素単位円上の λ_anti = −e^{2iθ} 軌跡（R:0→1）と根 ---
    fig, ax = plt.subplots(figsize=(5.6, 5.6))
    circ = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(circ), np.sin(circ), color="0.7", lw=1.0, zorder=1)
    Rs = np.linspace(0.0, 1.0, 600)
    xs, ys = [], []
    for R in Rs:
        _, la = p5.exact_eigenvalues(float(R))
        xs.append(la.real); ys.append(la.imag)
    ax.plot(xs, ys, color="tab:purple", lw=2.0, zorder=2,
            label=r"$\lambda_{anti}=-e^{2i\theta}$ (R:0$\to$1)")
    # λ_sym=1 固定点
    ax.plot([1], [0], "s", color="tab:gray", ms=9, zorder=3,
            label=r"$\lambda_{sym}=1$ (fixed)")
    # 端点と根
    for r in ladder:
        la = complex(r["lambda_a_re"], r["lambda_a_im"])
        if r["stat"].startswith("ボゾン"):
            ax.plot(la.real, la.imag, "o", color="tab:blue", ms=13, zorder=4,
                    label="boson end R=1")
        elif r["stat"].startswith("フェルミオン"):
            ax.plot(la.real, la.imag, "o", color="tab:red", ms=13, zorder=4,
                    label="fermion end R=0")
        else:
            ax.plot(la.real, la.imag, "^", color="tab:green", ms=9, zorder=4)
            ax.annotate(f"{r['n']}/{r['m']}", (la.real, la.imag),
                        textcoords="offset points", xytext=(6, 4), fontsize=8)
    ax.plot([], [], "^", color="tab:green", label="anyonic roots (n/m)")
    ax.set_xlabel(r"Re $\lambda$"); ax.set_ylabel(r"Im $\lambda$")
    ax.set_aspect("equal"); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
    ax.grid(True, alpha=0.2)
    ax.set_title("Statistics ladder on the unit circle")
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "N3_fig1_eigenvalue_unit_circle.svg"))
    fig.savefig(os.path.join(outdir, "N3_fig1_eigenvalue_unit_circle.png"), dpi=150)
    plt.close(fig)

    # --- 図2: R vs 交換位相（λ_anti の偏角/2π）と有理量子化点 ---
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ph = []
    for R in Rs:
        _, la = p5.exact_eigenvalues(float(R))
        ph.append((math.atan2(la.imag, la.real) / (2 * math.pi)) % 1.0)
    ax.plot(Rs, ph, color="tab:purple", lw=2.0,
            label=r"exchange phase = arg$(\lambda_{anti})/2\pi$")
    for r in ladder:
        Rr = r["R"]; pr = r["phase_cycles"]
        if r["stat"].startswith("ボゾン"):
            ax.plot(Rr, pr, "o", color="tab:blue", ms=11, zorder=5)
        elif r["stat"].startswith("フェルミオン"):
            ax.plot(Rr, pr, "o", color="tab:red", ms=11, zorder=5)
        else:
            ax.plot(Rr, pr, "^", color="tab:green", ms=9, zorder=5)
            ax.annotate(f"{r['n']-r['m']}/{r['n']}", (Rr, pr),
                        textcoords="offset points", xytext=(5, -12), fontsize=8)
    ax.set_xlabel("exchange coefficient  R")
    ax.set_ylabel(r"exchange phase / 2$\pi$")
    ax.set_title("Finite-order roots quantize the exchange phase (rational)")
    ax.set_ylim(0.45, 1.03)
    ax.grid(True, alpha=0.3)
    ax.text(0.03, 0.53, "fermion R=0", color="tab:red", fontsize=8,
            ha="left", va="center")
    ax.text(0.97, 0.97, "boson R=1", color="tab:blue", fontsize=8,
            ha="right", va="center")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "N3_fig2_statistics_ladder.svg"))
    fig.savefig(os.path.join(outdir, "N3_fig2_statistics_ladder.png"), dpi=150)
    plt.close(fig)
    print(f"図を生成: {outdir}/N3_fig1_eigenvalue_unit_circle.(svg|png), "
          f"N3_fig2_statistics_ladder.(svg|png)")


if __name__ == "__main__":
    raise SystemExit(main())
