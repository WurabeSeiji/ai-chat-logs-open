# -*- coding: utf-8 -*-
"""
円周住民：ゲージ間隔に対する測地線長の不等間隔（二次元・円弧）
設定（余分な仮定なし）
  成分      x_n(θ) = iA cos(n(θ-θ0)/2)          （本当の二乗を使う）
  殻        Σ x_n^2 = R^2  ->  R(θ)^2 = -A^2 Σ cos^2(n(θ-θ0)/2)  （負：R は虚数＝振幅）
  ゲージ    最高次成分 n=N の節 θ_k = θ0 + (2k+1)π/N  （位相で等間隔 2π/N）
  測地線長  ds = R(θ) dθ  ->  Δs_k = ∫_{θ_k}^{θ_{k+1}} R dθ = i ∫ |R| dθ
  読出し    Δs_k / <Δs>   （i は共通因子なので比は実数）
曲率は使わない。外は持ち込まない。
"""
import numpy as np, json

A = 1.0
THETA0 = 0.0
M = 200000            # 一周あたりの分割数

SYSTEMS = {
    "all_N2":   list(range(1, 3)),
    "all_N4":   list(range(1, 5)),
    "all_N18":  list(range(1, 19)),
    "odd_N17":  list(range(1, 18, 2)),
    "even_N18": list(range(2, 19, 2)),
}

def R2(theta, ns):
    phi = theta - THETA0
    return -A**2 * sum(np.cos(n*phi/2)**2 for n in ns)   # 本当の二乗の和

def run(ns):
    N = max(ns)
    # 一周 2π を最高次の節で N 区間に切る（R^2 は周期 2π）
    k = np.arange(N+1) - (N+1)//2
    marks = THETA0 + (2*k + 1)*np.pi/N        # 最高次成分の節そのもの（一周分 N 区間）
    ds_im = []
    for a, b in zip(marks[:-1], marks[1:]):
        t = np.linspace(a, b, M//N + 1)
        absR = np.sqrt(-R2(t, ns))                  # |R| = R/i
        ds_im.append(np.trapezoid(absR, t))             # Δs_k / i
    ds_im = np.array(ds_im)
    ratio = ds_im/ds_im.mean()
    centers = 0.5*(marks[:-1]+marks[1:])
    return dict(N=N, ns=ns, marks=marks.tolist(), centers=centers.tolist(),
                ds_over_i=ds_im.tolist(), L_over_i=float(ds_im.sum()),
                ratio=ratio.tolist(), max_ratio=float(ratio.max()),
                min_ratio=float(ratio.min()), max_over_min=float(ratio.max()/ratio.min()))

if __name__ == "__main__":
    out = {"settings": dict(A=A, theta0=THETA0, M=M,
            definition="R^2 = Σx_n^2 = -A^2 Σcos^2(n(θ-θ0)/2); Δs_k = i∫|R|dθ over gauge interval 2π/N"),
           "results": {k: run(v) for k, v in SYSTEMS.items()}}
    # 曲線データ
    th = np.linspace(-np.pi, np.pi, 4001)
    np.savez("data_gauge_geodesic_s1.npz", theta=th,
             **{k: np.sqrt(-R2(th, v)) for k, v in SYSTEMS.items()})
    with open("results_gauge_geodesic_s1.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    for k, r in out["results"].items():
        print(k, "L/i=%.4f" % r["L_over_i"], "max=%.3f min=%.3f max/min=%.3f" %
              (r["max_ratio"], r["min_ratio"], r["max_over_min"]))
        print("  ratio:", np.round(r["ratio"], 3).tolist())
