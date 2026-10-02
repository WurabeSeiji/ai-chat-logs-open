#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
checks.py — 動作確認走行（run_all.sh）の結果を仕様の数値と照合する。
  1. test1a（δ=0, f=0, スピン 0）: 仕様 §3 の表（γ−1、a サイクル/周、θ_b/θ_a、ψ1 位相/周）との一致
  2. test2 （δ=0.065, f=0）        : 保存（同じ呼吸位相の ε̄・A・⟨m²⟩、呼吸一周平均の ε̄）
  3. test3s（δ=0.065, f=1, 全節点）: エネルギー収支（呼吸一周平均 ε̄ の変化 vs ∫P_rad dθ）と、
                                      ψ への放射反作用ステップの仕事 vs Larmor（差分の遅れの因子）
  4. 倍音の 1 節点の位相（節点時計で位相が折り返す |m| の下限）
数値だけを出す。合否判定は書かない（判定は木原氏）。
"""
import csv
import json
import math
import warnings

import numpy as np

warnings.simplefilter("ignore")
import em_two_body_node_map as E  # noqa: E402

SPEC = dict(gamma_m1=4.66e-5, a_cycles_per_orbit=37567.17, thb_over_tha=1836.2015, psi1_phase_over_pi=1 - 4.0e-5,
            xi1=274.1811, omX=2.661704e-5, F_a=0.99997339)


def load(path):
    return list(csv.DictReader(open(path)))


def col(rows, name):
    return np.array([float(r[name]) for r in rows])


def main():
    print("=" * 78)
    print("1. test1a (δ=0, f=0, spins 0) vs 仕様 §3")
    s = json.load(open("results/test1a_delta0_f0_spin0_summary.json"))
    print(f"   ξ1              program {s['init']['xi1']:.6f}   spec {SPEC['xi1']}")
    print(f"   ω_X             program {s['init']['omX0']:.6e}   spec {SPEC['omX']:.6e}")
    print(f"   γ−1 (mean)      program {s['gamma_mean'] - 1:.4e}   spec {SPEC['gamma_m1']:.3e}")
    print(f"   a cycles/orbit  program {s['a_cycles_per_orbit']:.3f}   spec {SPEC['a_cycles_per_orbit']}")
    print(f"   θ_b/θ_a         program {s['theta_b_over_theta_a']:.7f}   spec {SPEC['thb_over_tha']}")
    print(f"   ψ1 phase/orbit/π program {s['psi1_phase_per_orbit_over_pi']:.7f}   spec {SPEC['psi1_phase_over_pi']:.7f}")
    print(f"   F_a (first node) program {s['first']['Fa']:.8f}   spec {SPEC['F_a']}")
    print(f"   |S1|, |S2| last  {s['last']['S1_abs']:.2e}, {s['last']['S2_abs']:.2e}   norm−1 {s['last']['norm'] - 1:.1e}")
    print(f"   phase/node: m=1 {s['init']['phase_per_node_m1']:.4f} rad, m=31 {s['init']['phase_per_node_mM']:.1f} rad, "
          f"first |m| with phase>π: {s['init']['smallest_m_with_phase_over_pi']}")

    print("=" * 78)
    print("2. test2 (δ=0.065, f=0) 保存")
    rows = [r for r in load("results/test2_delta065_f0.csv") if int(r["k"]) % 62 == 0]   # 同じ呼吸位相（最終行 k=6199 を除く）
    eb, A, m2 = col(rows, "ebar"), col(rows, "A"), col(rows, "m2")
    print(f"   ε̄ same phase: first {eb[0]:.12e} last {eb[-1]:.12e}  max|Δ|/|ε̄| = {np.max(np.abs(eb - eb[0])) / abs(eb[0]):.2e}")
    print(f"   A: first {A[0]:.9f} last {A[-1]:.9f}   ⟨m²⟩: first {m2[0]:.12f} last {m2[-1]:.12f}")
    s2 = json.load(open("results/test2_delta065_f0_summary.json"))
    print(f"   γ_mean−1 = {s2['gamma_mean'] - 1:.4e}  （倍音の混合 ⟨m²⟩=(1+δ)² による。1/√⟨m²⟩−1 = {1 / math.sqrt(m2[0]) - 1:.4e}）")
    print(f"   e_K0 = {s2['init']['e_K0']:.6f}（Kepler: ⟨m²⟩−1 = {m2[0] - 1:.6f}）, A0 = {s2['init']['A0']:.4f}, q0 = {s2['init']['q0']:.5f}, ξ̄_c0 = {s2['init']['xi_c0']:.4f}")
    print(f"   呼吸級数: 当てはめ残差 ξ/ξ̄_c = {s2['init']['breathing_fit_resid_xi']:.1e}, 級数に沿った ε̄ の揺らぎ = {s2['init']['breathing_energy_wobble']:.1e}, "
          f"q0^(N+1) = {s2['init']['q0_pow_Np1']:.1e}（基準 {s2['init']['readout_precision']:.1e}）")
    rows_all = load("results/test2_delta065_f0.csv")
    eb_all = col(rows_all, "ebar")
    print(f"   走行中の ε̄ の全記録点の幅/|ε̄|（呼吸位相が違う点も含む）= {(eb_all.max() - eb_all.min()) / abs(eb_all[0]):.2e}")

    print("=" * 78)
    print("3. test3s (δ=0.065, f=1, 全節点) エネルギー収支")
    rows = load("results/test3s_delta065_f1_every_node.csv")
    k, eb, dth, el = col(rows, "k"), col(rows, "ebar"), col(rows, "dth"), col(rows, "eloss_em")
    m2, A, S1 = col(rows, "m2"), col(rows, "A"), col(rows, "S1_abs")
    n = len(rows)
    ncyc = n // 62
    cyc_mean = np.array([eb[i * 62:(i + 1) * 62].mean() for i in range(ncyc)])
    E_rad = np.concatenate([[0.0], np.cumsum(0.5 * (el[1:] + el[:-1]) * dth[1:])])
    E_rad_cyc = np.array([E_rad[(i + 1) * 62 - 1] for i in range(ncyc)])
    print(f"   呼吸一周平均 ε̄: 第1周 {cyc_mean[0]:.10e}  最終周 {cyc_mean[-1]:.10e}  Δ = {cyc_mean[-1] - cyc_mean[0]:.4e}")
    print(f"   ∫P_rad dθ（Larmor、差分 D''）最終周まで = {E_rad_cyc[-1] - E_rad_cyc[0]:.4e}")
    print(f"   比 Δ⟨ε̄⟩/∫P_rad = {(cyc_mean[-1] - cyc_mean[0]) / (E_rad_cyc[-1] - E_rad_cyc[0]):.4f}")
    print(f"   ε̄ の呼吸一周内の変動幅/|ε̄|（第 1 周。級数の精度の指標）= {(eb[:62].max() - eb[:62].min()) / abs(eb[0]):.2e}")
    print(f"   ⟨m²⟩ first/last {m2[0]:.7f} {m2[-1]:.7f}   A first/last {A[0]:.5f} {A[-1]:.5f}   mean|S1| {S1.mean():.4f}")
    # ψ への V ステップの仕事 vs Larmor（走行を再現して計算）
    C = E.build_constants()
    S = E.init_state(C, delta=0.065, f=1.0, sa=(0, 0, 1), sb=(0, 0, -1))
    W_V = W_rad = W_rad_x = 0.0
    Dh = []
    for _ in range(62 * 100):
        psi = S.psi.copy()
        sa_z, sb_z = S.sa[2], S.sb[2]
        xi_k, pi_k, omX, _, _ = E.breathing_eval(C, S.T, S.B, S.xi_c, S.zeta_c, S.omX0)
        dt = C.dthX / omX
        eps_m, _ = E.eps_modes(C, 1.0 / xi_k, pi_k, sa_z, sb_z)
        S1_, _ = E.quad_relations(psi)
        Dh.append(xi_k * S1_)
        Dh = Dh[-4:]
        if len(Dh) == 4:
            D3 = (Dh[3] - 3 * Dh[2] + 3 * Dh[1] - Dh[0]) / dt ** 3
            D2 = (Dh[3] - 2 * Dh[2] + Dh[1]) / dt ** 2
            Em = C.mu_over_2ma * eps_m
            rates = -1j * (Em[:-1] - Em[1:])
            comps = psi[:-1] * np.conj(psi[1:])
            D2x = xi_k * np.sum(comps * rates ** 2)
            Vpsi = -C.rr_em * xi_k * (np.conj(D3) * E.shift(psi, -1) + D3 * E.shift(psi, +1))
            ph = np.exp(-1j * Em * dt)
            W_V += float(np.sum((np.abs((psi - 1j * C.mu_over_2ma * dt * Vpsi) * ph) ** 2 - np.abs(psi * ph) ** 2) * eps_m))
            W_rad += -C.eloss_em * abs(D2) ** 2 * dt
            W_rad_x += -C.eloss_em * abs(D2x) ** 2 * dt
        E.node_map(C, S)
    print(f"   ψ への反作用ステップの仕事 W_V = {W_V:.4e}")
    print(f"   Larmor（差分 D''）= {W_rad:.4e}  比 W_V/W_rad = {W_V / W_rad:.4f}")
    print(f"   Larmor（解析 D''）= {W_rad_x:.4e}  比 W_V/W_rad_x = {W_V / W_rad_x:.4f}   （4 点後退差分の遅れ 1.5 節点の因子）")
    Em = C.mu_over_2ma * E.eps_modes(C, 1 / 311.0, 0.0, 1, -1)[0]
    dt0 = C.dthX / S.omX0
    print(f"   双極子成分の 1 節点の位相: (E0−E1)Δθ = {(Em[C.M] - Em[C.M + 1]) * dt0:.4f} rad, (E1−E2)Δθ = {(Em[C.M + 1] - Em[C.M + 2]) * dt0:.4f} rad")

    print("=" * 78)
    print("4. test3 (δ=0.065, f=1, 呼吸 1000 周) 方向")
    s3 = json.load(open("results/test3_delta065_f1_summary.json"))
    print(f"   ⟨m²⟩ first {s3['init']['psi_w0'][C.M + 2] * 4 + s3['init']['psi_w0'][C.M + 1]:.7f} → last {s3['last']['m2']:.7f}")
    print(f"   |S2| last {s3['last']['S2_abs']:.5f}（初期 0.0671）  q = A/ξ̄_c last {s3['last']['q']:.6f}（初期 {s3['init']['q0']:.6f}）")
    print(f"   ε̄ last {s3['last']['ebar']:.9e}（初期 {s3['first']['ebar']:.9e}）  norm−1 {s3['last']['norm'] - 1:.1e}  速度 {s3['nodes_per_s']:.0f} nodes/s")


if __name__ == "__main__":
    main()
