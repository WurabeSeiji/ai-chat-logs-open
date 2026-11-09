# -*- coding: utf-8 -*-
"""核 F(d) の形を二つの包絡で比べる（確認計算。実験ではない）。
sinc 型：現行の c_m = (2/π) sin(mπφ/360)/m、φ = 1°、M = 100（c_m はほぼ一定なので自己相関は Dirichlet 核、副峰は 1/d で減る）。
cos² 型：c_m ∝ C(2k, k−m)、k = M = 100（自己相関は cos^{2k}(d/2) 型で副峰がない、幅 ≈ 2/√k）。"""
import json, math
M = 100; phi = 1.0
cs_sinc = [(2 / math.pi) * math.sin(m * math.pi * phi / 360) / m for m in range(1, M + 1)]
lg = [0.0] * (2 * M + 1)
for n in range(1, 2 * M + 1): lg[n] = lg[n - 1] + math.log(n)
cs_cos = [math.exp(lg[2 * M] - lg[M - m] - lg[M + m] - 2 * M * math.log(2) + math.log(2)) for m in range(1, M + 1)]   # 2^{1−2k} C(2k, k−m)
def F(cs, d):
    r = math.radians(d); num = sum(c * c * math.cos((m + 1) * r) for m, c in enumerate(cs)); den = sum(c * c for c in cs); return num / den
ds = [0, 1, 2, 3, 5, 8, 10, 15, 20, 30, 45, 90, 180]
out = {"M": M, "phi_deg": phi, "d_deg": ds, "F_sinc": [F(cs_sinc, d) for d in ds], "F_cos2": [F(cs_cos, d) for d in ds],
       "sinc_max_abs_F_beyond_30deg": max(abs(F(cs_sinc, d)) for d in range(30, 181)), "cos2_max_abs_F_beyond_30deg": max(abs(F(cs_cos, d)) for d in range(30, 181)),
       "half_width_deg_sinc": next(d / 10 for d in range(0, 3600) if F(cs_sinc, d / 10) < 0.5), "half_width_deg_cos2": next(d / 10 for d in range(0, 3600) if F(cs_cos, d / 10) < 0.5),
       "two_over_sqrtk_deg": math.degrees(2 / math.sqrt(M))}
json.dump(out, open("results/kernel_profiles.json", "w"), ensure_ascii=False, indent=2)
print("d°      F_sinc      F_cos2"); [print(f"{d:4d}  {a:+.5f}  {b:+.5f}") for d, a, b in zip(ds, out["F_sinc"], out["F_cos2"])]
print("半値幅（度）sinc", out["half_width_deg_sinc"], "cos2", out["half_width_deg_cos2"], "2/√k =", round(out["two_over_sqrtk_deg"], 2))
print("30° 以遠の |F| 最大：sinc", f"{out['sinc_max_abs_F_beyond_30deg']:.4f}", " cos2", f"{out['cos2_max_abs_F_beyond_30deg']:.2e}")
