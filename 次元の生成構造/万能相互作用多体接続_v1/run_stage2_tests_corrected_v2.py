#!/usr/bin/env python3
"""段階2 単体検証の修正再実験 v2——判定基準の誤設定を修正した再検定

修正の記録（仮説→反証→修正の規律。v1 の誤りを隠さず記載）:
    v1-U2 の誤り: 「点ごと場閉塞の合計」を全力学の不変量として判定したが、
        これは v1 アーキテクチャ（スライス独立線形部）の不変量ではない
        （各スライスが別の直交行列で回るため、スライス間交差モーメントは
        正当に変化する）。計器の誤りであり、力学の欠陥ではない。
        正しい不変量: (a) 全ノルム（線形=直交・頂点=定理、両部で厳密保存）
        (b) 頂点流単独での点ごと場閉塞（一意化定理＋RK4積分器の検定）
        (c) 線形部単独でのスライス毎閉塞 c_k^T c_k（Cayley直交性の検定）。
    v1-U3 の誤り: 閾値 1e-6 に根拠がなかった。修正: 閾値を同一実験系で
        独立に測った点火則（E1a: rate=C·f^p）から導出する。

判定基準 v2（実行前固定）:
    U2a 全ノルム: 結合力学 T=200 で相対ドリフト ≤1e-10。
    U2b 頂点流単独: T=200（頂点のみ適用）で点ごと場閉塞 max|Σ_e w_e[n]²|
        の変化 ≤1e-10（初期値は非零でよい——変化量を測る）。
    U2c 線形部単独: T=200 で各占有スライスの |c_k^T c_k| の変化 ≤1e-12。
    U3′ 生成の存在と法則整合:
        OFF: |Δf_seed| ≤1e-14（線形無生成）
        ON:  Δf_seed(ON)/max(Δf_seed(OFF),1e-20) ≥1e3（分離）かつ
             Δf_seed(ON) / Δf_pred ∈ [0.1, 10]（点火則整合）。
             Δf_pred = f₀·(exp(C·f₀^p·T)−1)、C,p は本実験の E1a 掃引から取得。

使い方: python3 run_stage2_tests_corrected_v2.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

spec_e = importlib.util.spec_from_file_location(
    "s2v1", HERE / "run_stage2_vertex_engine_v1.py")
s2 = importlib.util.module_from_spec(spec_e)
sys.modules[spec_e.name] = s2
spec_e.loader.exec_module.__self__ if False else None
# 実行部を走らせずにクラスだけ使うため、__main__ ガード付きモジュールとしてロード
spec_e.loader.exec_module(s2)

abl = s2.abl
gen3 = s2.gen3
VertexEngine = s2.VertexEngine


def main() -> None:
    t0 = time.time()
    n, m, nreg, T = 5, 10, 5, 200
    out = {"correction_record": "v1-U2は非不変量を判定（計器誤り）、v1-U3閾値1e-6は無根拠。"
                                  "v2は不変量の分離検定＋点火則導出閾値に修正。"}

    _, v, _, _, _, _, _, Z0c, wp0 = abl.build_init(n, False)
    r = gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])

    def mixed_C0(delta):
        C0 = np.zeros((m, nreg), complex)
        C0[:, 2] = Z0c
        C0[:, 1] = delta * seed_state
        return C0

    wps = {1: np.random.default_rng(92001).normal(size=m), 2: wp0.copy()}

    # ---- U2a: 結合力学の全ノルム ----
    eng = VertexEngine(n, mixed_C0(0.01), wps, vertex_on=True)
    n0 = eng.diagnostics()["norm"]
    for _ in range(T):
        eng.step()
    n1 = eng.diagnostics()["norm"]
    u2a = abs(n1 - n0) / n0
    out["U2a_norm_drift_rel"] = float(u2a)
    out["U2a_pass"] = bool(u2a <= 1e-10)
    print(f"  U2a 結合・全ノルム: 相対ドリフト={u2a:.2e} → {out['U2a_pass']}")

    # ---- U2b: 頂点流単独の点ごと場閉塞 ----
    eng = VertexEngine(n, mixed_C0(0.1), wps, vertex_on=True)   # 大きめの種で頂点を強く駆動
    W0 = np.fft.ifft(eng.C, axis=1) * nreg
    cl0 = np.abs(np.sum(W0 ** 2, axis=0))
    for _ in range(T):
        eng._nonlinear()                                          # 頂点のみ
    W1 = np.fft.ifft(eng.C, axis=1) * nreg
    cl1 = np.abs(np.sum(W1 ** 2, axis=0))
    u2b = float(np.max(np.abs(cl1 - cl0)))
    out["U2b_vertex_only_closure_drift"] = u2b
    out["U2b_pass"] = bool(u2b <= 1e-10)
    print(f"  U2b 頂点流単独・点ごと閉塞: max変化={u2b:.2e}（初期値max={cl0.max():.2e}） → {out['U2b_pass']}")

    # ---- U2c: 線形部単独のスライス毎閉塞 ----
    eng = VertexEngine(n, mixed_C0(0.01), wps, vertex_on=False)
    sl0 = {k: complex(eng.C[:, k] @ eng.C[:, k]) for k in (1, 2)}
    for _ in range(T):
        eng.step()
    sl1 = {k: complex(eng.C[:, k] @ eng.C[:, k]) for k in (1, 2)}
    u2c = max(abs(sl1[k] - sl0[k]) for k in (1, 2))
    out["U2c_linear_only_slice_closure_drift"] = float(u2c)
    out["U2c_pass"] = bool(u2c <= 1e-12)
    print(f"  U2c 線形単独・スライス毎閉塞: max変化={u2c:.2e} → {out['U2c_pass']}")

    # ---- E1a 再測定（U3′ の閾値導出のため同一実験系で） ----
    rows = []
    for delta in (1e-3, 3e-3, 1e-2, 3e-2):
        eng = VertexEngine(n, mixed_C0(delta), wps, vertex_on=True)
        po = []
        for _ in range(60):
            eng.step()
            po.append(eng.diagnostics()["f_seed"])
        po = np.array(po)
        tt = np.arange(5, 40, dtype=float)
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, np.log(po[5:40]), rcond=None)
        rows.append({"delta": delta, "f0": float(po[0]), "rate": float(coef[0])})
    ln_f = np.log([q["f0"] for q in rows])
    ln_r = np.log([q["rate"] for q in rows])
    A = np.vstack([ln_f, np.ones_like(ln_f)]).T
    coef, _, _, _ = np.linalg.lstsq(A, ln_r, rcond=None)
    p_exp = float(coef[0])
    Cm = float(np.exp(coef[1]))
    out["E1a_refit"] = {"p": p_exp, "C": Cm}
    print(f"  E1a 再フィット: p={p_exp:.3f} C={Cm:.3f}")

    # ---- U3′: 生成の存在と法則整合 ----
    delta = 0.01
    eng_on = VertexEngine(n, mixed_C0(delta), wps, vertex_on=True)
    f0 = eng_on.diagnostics()["f_seed"]
    for _ in range(T):
        eng_on.step()
    df_on = abs(eng_on.diagnostics()["f_seed"] - f0)
    eng_off = VertexEngine(n, mixed_C0(delta), wps, vertex_on=False)
    g0 = eng_off.diagnostics()["f_seed"]
    for _ in range(T):
        eng_off.step()
    df_off = abs(eng_off.diagnostics()["f_seed"] - g0)
    rate_pred = Cm * f0 ** p_exp
    df_pred = f0 * (np.exp(rate_pred * T) - 1.0)
    sep = df_on / max(df_off, 1e-20)
    law_ratio = df_on / df_pred
    out["U3"] = {"df_on": float(df_on), "df_off": float(df_off),
                  "df_pred": float(df_pred), "separation": float(sep),
                  "law_ratio": float(law_ratio)}
    out["U3_pass"] = bool(df_off <= 1e-14 and sep >= 1e3 and 0.1 <= law_ratio <= 10.0)
    print(f"  U3′ 生成: Δf ON={df_on:.2e} OFF={df_off:.2e} 予測={df_pred:.2e} "
          f"分離={sep:.1e}倍 法則比={law_ratio:.2f} → {out['U3_pass']}")

    ok = all(out[k] for k in ("U2a_pass", "U2b_pass", "U2c_pass", "U3_pass"))
    out["all_pass"] = bool(ok)
    out["runtime_sec"] = time.time() - t0
    print(f"\n修正再実験 判定: {'ALL PASS' if ok else 'FAIL あり'}")
    (HERE / "stage2_tests_corrected_result_v2.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
