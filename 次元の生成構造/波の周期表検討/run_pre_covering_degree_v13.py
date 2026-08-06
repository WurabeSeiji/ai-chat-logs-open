#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v13: 時計被覆度の忠実度計——状態周期/観測周期の比＝スピン二値性の測定

原理（事前記録）: C論文の忠実度計（F_124 vs F_248: 状態周期=観測周期×2の実測）
を本系の種に適用。観測量場 s=Im(b̄a) の自己相関から観測周期 P_obs を検出し、
振幅重なり A(J)=⟨ψ(0)|ψ(J)⟩/‖ψ‖² を P_obs / 2P_obs で比較。
  被覆度2（フェルミオン的候補）: |A(P_obs)| 小 or Re A(P_obs)<0、|A(2P_obs)|≈1
  被覆度1（ボゾン的）: |A(P_obs)|≈1
対象: 孤立純粋+1／孤立純粋+2／中性m=0束／帯電census(D)。J≤3000。
判定 H_deg: 種ごとに被覆度を記録（差が出れば周期表スピン列の第二成分が実測で埋まる）。

使い方: python3 run_pre_covering_degree_v13.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v13", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_MAX = 3000

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    eta = 2 * np.pi * np.arange(ne) / ne

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)
    def shift_eta(v, dm):
        return (v.reshape(shape) * np.exp(1j * dm * eta)[None, :]).reshape(v.shape)

    a0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a0s = a0 + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
    pow0 = float(np.sum(np.abs(a0s) ** 2) + np.sum(np.abs(b0) ** 2))
    a1 = project_eta(a0s, {1}); b1 = project_eta(b0, {1})
    pw = float(np.sum(np.abs(a1) ** 2) + np.sum(np.abs(b1) ** 2))
    sc = np.sqrt(pow0 / pw); a1 *= sc; b1 *= sc
    neu_a = project_eta(a0s, {0}); neu_b = project_eta(b0, {0})
    pwn = float(np.sum(np.abs(neu_a) ** 2) + np.sum(np.abs(neu_b) ** 2))
    if pwn < 1e-9:
        neu_a = project_eta(v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S, {0})
        neu_b = project_eta(v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S, {0})
        pwn = float(np.sum(np.abs(neu_a) ** 2) + np.sum(np.abs(neu_b) ** 2))
    scn = np.sqrt(pow0 / pwn); neu_a *= scn; neu_b *= scn

    cases = {"pure+1孤立": (a1.copy(), b1.copy()),
             "pure+2孤立": (shift_eta(a1, 1), shift_eta(b1, 1)),
             "中性m=0束": (neu_a, neu_b),
             "帯電census(D)": (a0s.copy(), b0.copy())}
    out = {"J_MAX": J_MAX, "cases": {}}
    for name, (a, b) in cases.items():
        a0_, b0_ = a.copy(), b.copy()
        s0 = np.imag(np.conj(b0_) * a0_)
        s0c = s0 - s0.mean()
        nrm = float(np.sum(np.abs(a0_) ** 2) + np.sum(np.abs(b0_) ** 2))
        A_series, C_series = [], []
        for j in range(J_MAX):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            ov = (np.vdot(a0_, a) + np.vdot(b0_, b)) / nrm
            A_series.append(complex(ov))
            sj = np.imag(np.conj(b) * a)
            sjc = sj - sj.mean()
            C_series.append(float(np.dot(s0c, sjc) /
                                   max(np.linalg.norm(s0c) * np.linalg.norm(sjc), 1e-300)))
        A_series = np.array(A_series); C = np.array(C_series)
        # 観測周期: C の最初の強ピーク（J>10）
        Js = np.arange(1, J_MAX + 1)
        mask = Js > 10
        p_obs = int(Js[mask][np.argmax(C[mask])])
        Cp = float(C[p_obs - 1])
        def A_at(j):
            j = min(j, J_MAX) - 1
            return A_series[j]
        A1 = A_at(p_obs); A2 = A_at(2 * p_obs)
        deg = None
        if abs(A1) > 0.8 and A1.real > 0:
            deg = 1
        elif (abs(A1) < 0.5 or A1.real < 0) and abs(A2) > 0.8 and A2.real > 0:
            deg = 2
        print(f"{name}: P_obs={p_obs}（C={Cp:.3f}） "
              f"A(P)={A1.real:+.3f}{A1.imag:+.3f}i(|{abs(A1):.3f}|) "
              f"A(2P)={A2.real:+.3f}{A2.imag:+.3f}i(|{abs(A2):.3f}|) → 被覆度={deg}")
        out["cases"][name] = {"P_obs": p_obs, "C_at_P": Cp,
            "A_P": [A1.real, A1.imag], "A_2P": [A2.real, A2.imag],
            "degree": deg,
            "C_top5": sorted([[int(j), float(c)] for j, c in zip(Js[mask], C[mask])],
                              key=lambda x: -x[1])[:5]}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_covering_degree_result_v13.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_covering_degree_result_v13.json")

if __name__ == "__main__":
    main()
