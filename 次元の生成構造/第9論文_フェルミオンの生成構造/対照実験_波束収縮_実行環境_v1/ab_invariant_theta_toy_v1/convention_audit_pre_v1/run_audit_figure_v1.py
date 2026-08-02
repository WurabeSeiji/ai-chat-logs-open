#!/usr/bin/env python3
"""規約監査の図 v1（論文B用）: マスク×状態の R 行列と有理数注記"""
import json, importlib.util, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("toy_fig", HERE.parent / "run_ab_invariant_theta_toy_v1.py")
toy = importlib.util.module_from_spec(spec); sys.modules[spec.name] = toy; spec.loader.exec_module(toy)
plt = toy.base.plt
data = json.loads((HERE / "convention_audit_pre_result_v1.json").read_text(encoding="utf-8"))
rows = [r for r in data["readout_rows"] if "spread" not in str(r["state"])]
masks = sorted({r["mask"] for r in rows})
states = ["B63_equal", "B63_inverse_k", "odds_5_to_63_equal", "evens_equal"]
M = np.array([[next(r["R_measured"] for r in rows if r["mask"] == m and r["state"] == s)
               for s in states] for m in masks])
fig, ax = plt.subplots(figsize=(9.5, 5), constrained_layout=True)
im = ax.imshow(M, cmap="viridis", vmin=0, vmax=1)
ax.set_xticks(range(len(states))); ax.set_xticklabels(states, rotation=20, ha="right", fontsize=8)
ax.set_yticks(range(len(masks))); ax.set_yticklabels(masks, fontsize=8)
for i, m in enumerate(masks):
    for j, s in enumerate(states):
        r = next(r for r in rows if r["mask"] == m and r["state"] == s)
        label = r["nearest_rational"] if "PASS" in str(r["C2_rational_at_equal_weight"]) else f"{r['R_measured']:.4f}"
        ax.text(j, i, label, ha="center", va="center", fontsize=8,
                color="white" if M[i, j] < 0.6 else "black")
ax.set_title("Convention audit: equal-weight R degenerates to exact rationals under every mask\n(amplitude-deformed column stays irrational)")
fig.colorbar(im, ax=ax, shrink=0.8, label="R")
for ext in ("png", "svg"):
    fig.savefig(HERE / f"audit_figure_v1.{ext}", dpi=160)
print("audit figure generated")
