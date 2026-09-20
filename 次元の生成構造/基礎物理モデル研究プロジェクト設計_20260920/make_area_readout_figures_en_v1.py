#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_area_readout_figures_en_v1.py

English-label PNG versions of fig01-fig03 of the first thought experiment
("Tracing the Equivalence Principle Down to Two Nameless Degrees of Freedom").

It reuses run_area_readout_experiments_v1.py in the same folder and does NOT rerun
the experiments: fig01 and fig02 are deterministic, fig03 is drawn from the stored
area_readout_experiments_results_v1.json, so the figures show the same data as the
Japanese ones. No existing file is overwritten.

Output (same folder):
  fig01_area_readout_three_classes_en_v1.png
  fig02_area_readout_closure_staircase_en_v1.png
  fig03_area_readout_metric_from_orbit_en_v1.png

Requires: numpy, mpmath, matplotlib.
Run: python3 make_area_readout_figures_en_v1.py
"""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("area_readout", HERE / "run_area_readout_experiments_v1.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)                      # main() is not executed
if not m.HAVE_MPL:
    raise SystemExit("matplotlib is required to draw the PNG figures.")

# Japanese label in the main script -> English label
EN = {
    "(a) 楕円型 0<κ<4：一般の記述では楕円": "(a) Elliptic type, 0<κ<4: an ellipse in a generic description",
    "(b) 双曲型 κ<0：a²−b²=C²（破線は零錐 a=±b）": "(b) Hyperbolic type, κ<0: a²−b²=C² (dashed: null cone a=±b)",
    "(c) 放物型 κ=0：a が一定、b が k に比例して増える": "(c) Parabolic type, κ=0: a constant, b grows linearly in k",
    "(d) (a) の軌道に沿った二つの読み出し（初期値で規格化）": "(d) Two readouts along orbit (a), normalized to the first value",
    "素朴な a²+b²（保存しない）": "naive a²+b² (not conserved)",
    "面積読み出し ω(X_k, X_{k+1})（保存）": "area readout ω(X_k, X_{k+1}) (conserved)",
    "上界 1/δ ≈ 2πC/ε": "upper bound 1/δ ≈ 2πC/ε",
    "ρ = C/ε（曲率半径 ÷ 分解能）": "ρ = C/ε (curvature radius ÷ resolution)",
    "初回帰還数 n（＝有効な閉路数）": "first-return number n (= effective closure number)",
    "有限分解能のもとでの閉路数：n は θ/2π の連分数近似分母だけを取る":
        "Closure number under finite resolution: n takes only convergent denominators of θ/2π",
    "n≥3：相対誤差の最大値（200 試行）": "n ≥ 3: maximum relative error (200 trials)",
    "n≧3：相対誤差の最大値（200 試行）": "n ≥ 3: maximum relative error (200 trials)",
    "事前に宣言した許容値 1e-10": "pre-declared tolerance 1e-10",
    "n=2：S=−I で面積の自己読み出しが\n恒等的にゼロ（計量が定義できない）":
        "n=2: S=−I, the area self-readout vanishes\nidentically (no metric can be defined)",
    "閉路数 n": "closure number n",
    "計量読み出しの相対誤差": "relative error of the metric readout",
    "閉軌道の二次モーメントから計量を読む": "Reading the metric from the second moment of a closed orbit",
}

m.JA, m.FP = False, None                        # use matplotlib's default font
m.L = lambda ja, en: EN.get(ja, en)
m.FIG1 = "fig01_area_readout_three_classes_en_v1"
m.FIG2 = "fig02_area_readout_closure_staircase_en_v1"
m.FIG3 = "fig03_area_readout_metric_from_orbit_en_v1"

m.make_fig1(m.figure1_data())

_, curves, rhos = m.closure_tests()
names = {"golden (sqrt5-1)/2": "α = (√5−1)/2 (golden ratio)", "sqrt2-1": "α = √2−1", "e-2": "α = e−2", "pi-3": "α = π−3"}
m.make_fig2(rhos, {names.get(k, k): v for k, v in curves.items()})

with open(HERE / "area_readout_experiments_results_v1.json", encoding="utf-8") as f:
    mr = json.load(f)["metric_readout"]
mr["per_n"] = {int(k): v for k, v in mr["per_n"].items()}     # JSON keys are strings
m.make_fig3(mr)

for name in (m.FIG1, m.FIG2, m.FIG3):
    print("written:", name + ".png")
