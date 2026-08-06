#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v13b 被覆度計の較正版: 連続回帰列での符号交代検定

較正の要点（事前記録）:
  v13は単一P_obsの1点判定で、位相の偶然のπ近傍と構造的ロックを区別できない。
  本版は観測回帰ピーク列 J_1<J_2<... を検出し、各ピークでの振幅位相 Φ_k=arg A(J_k)
  の系列で判定する:
    被覆度2: Re A が交代（−,+,−,+…）＝Φ_k が kπ にロック
    被覆度1: Re A が全て正（Φ_k≈0 mod 2π）
    被覆なし: Φ_k が連続ドリフト（有理ロックなし＝エネルギー巻きのみ）
  符号のみ使うため減衰に頑健。J=6000・ピーク条件 C>0.85（局所最大）。

使い方: python3 run_pre_covering_degree_v13b.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_v13b", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

S = 8.0; J_MAX = 6000; C_TH = 0.85

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    shape = (n, ne)
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)

    def single_winding(v):
        f = np.fft.fft(v.reshape(shape), axis=0, norm="ortho")
        f[0, :] = 0.0; f[n // 2:, :] = 0.0
        return np.fft.ifft(f, axis=0, norm="ortho").reshape(v.shape)
    def project_eta(v, m_set):
        f = np.fft.fft(v.reshape(shape), axis=1, norm="ortho")
        keep = np.isin(mm, list(m_set)); f[:, ~keep] = 0.0
        return np.fft.ifft(f, axis=1, norm="ortho").reshape(v.shape)

    a0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "A", scale=1.0)) * S
    b0 = single_winding(v1.make_bundle(sp, (30, 32, 34), "B", scale=1.0)) * S
    a0s = a0 + single_winding(v1.make_bundle(sp, (21,), "A", scale=1.0)) * (0.2 * S)
    pow0 = float(np.sum(np.abs(a0s) ** 2) + np.sum(np.abs(b0) ** 2))
    a1 = project_eta(a0s, {1}); b1 = project_eta(b0, {1})
    pw = float(np.sum(np.abs(a1) ** 2) + np.sum(np.abs(b1) ** 2))
    sc = np.sqrt(pow0 / pw); a1 *= sc; b1 *= sc
    neu_a = project_eta(v1.make_bundle(sp, (29, 31, 33), "A", scale=1.0) * S, {0})
    neu_b = project_eta(v1.make_bundle(sp, (29, 31, 33), "B", scale=1.0) * S, {0})
    pwn = float(np.sum(np.abs(neu_a) ** 2) + np.sum(np.abs(neu_b) ** 2))
    scn = np.sqrt(pow0 / pwn); neu_a *= scn; neu_b *= scn

    cases = {"中性m=0束": (neu_a, neu_b),
             "帯電census(D)": (a0s.copy(), b0.copy()),
             "pure+1孤立": (a1.copy(), b1.copy())}
    out = {"J_MAX": J_MAX, "C_TH": C_TH, "cases": {}}
    for name, (a, b) in cases.items():
        a0_, b0_ = a.copy(), b.copy()
        s0 = np.imag(np.conj(b0_) * a0_); s0c = s0 - s0.mean()
        nrm = float(np.sum(np.abs(a0_) ** 2) + np.sum(np.abs(b0_) ** 2))
        A, C = [], []
        for j in range(J_MAX):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            A.append(complex((np.vdot(a0_, a) + np.vdot(b0_, b)) / nrm))
            sj = np.imag(np.conj(b) * a); sjc = sj - sj.mean()
            C.append(float(np.dot(s0c, sjc) /
                            max(np.linalg.norm(s0c) * np.linalg.norm(sjc), 1e-300)))
        A = np.array(A); C = np.array(C)
        # ピーク列: 局所最大かつ C>閾（種ごとに実効閾=0.85×max後半）
        th = min(C_TH, 0.85 * float(C[10:].max()))
        peaks = []
        for j in range(11, J_MAX - 1):
            if C[j] > th and C[j] >= C[j - 1] and C[j] >= C[j + 1]:
                if not peaks or j - peaks[-1] > 5:
                    peaks.append(j)
            if len(peaks) >= 12:
                break
        rows = []
        signs = []
        print(f"{name}（実効閾C>{th:.2f}, ピーク{len(peaks)}個）:")
        for k, j in enumerate(peaks[:10]):
            ph = float(np.angle(A[j])); re = float(A[j].real)
            rows.append({"k": k + 1, "J": j + 1, "C": float(C[j]),
                          "ReA": re, "absA": float(abs(A[j])), "phase": ph})
            signs.append(np.sign(re))
            print(f"   k={k+1:2d} J={j+1:4d} C={C[j]:.3f} ReA={re:+.3f} "
                  f"|A|={abs(A[j]):.3f} Φ/π={ph/np.pi:+.3f}")
        # 判定
        verdict = "判定不能"
        if len(signs) >= 4:
            alt = all(signs[i] != signs[i + 1] for i in range(len(signs) - 1))
            allpos = all(s > 0 for s in signs)
            if alt and signs[0] < 0:
                verdict = "被覆度2（符号交代・πロック）"
            elif allpos:
                verdict = "被覆度1"
            else:
                phs = np.unwrap([r["phase"] for r in rows])
                drift = float(np.polyfit(range(len(phs)), phs, 1)[0])
                verdict = f"被覆なし/ドリフト（Φ増分={drift/np.pi:.3f}π/回帰）"
        print(f"   → {verdict}")
        out["cases"][name] = {"threshold": th, "peaks": rows, "verdict": verdict}
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_covering_degree_result_v13b.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_covering_degree_result_v13b.json")

if __name__ == "__main__":
    main()
