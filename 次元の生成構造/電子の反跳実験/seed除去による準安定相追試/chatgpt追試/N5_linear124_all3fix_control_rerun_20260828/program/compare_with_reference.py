# -*- coding: utf-8 -*-
"""Compare this control rerun (data/) with the reference package (../N5_linear124_all3fix_seedless_parentnorm_removed_20260828/data/).
Output: results/compare_with_reference.json"""
import os, json, csv, numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(os.path.dirname(HERE), "N5_linear124_all3fix_seedless_parentnorm_removed_20260828", "data")
MINE = os.path.join(HERE, "data"); R = os.path.join(HERE, "results"); os.makedirs(R, exist_ok=True)
def load(p):
    rows = list(csv.reader(open(p))); return rows[0], np.array([[float(x) for x in r] for r in rows[1:]])
out = {}
for name in ("baseline_linear124_phase_only_timeseries.csv", "treatment_linear124_amplitude_aware_timeseries.csv", "key_steps.csv"):
    h, a = load(os.path.join(REF, name)); h2, b = load(os.path.join(MINE, name))
    assert h == h2 and a.shape == b.shape
    d = np.abs(a - b); rel = d / np.maximum(np.abs(a), 1e-300)
    cols = {h[j]: {"max_abs_diff": float(np.nanmax(d[:, j])), "max_rel_diff": float(np.nanmax(np.where(np.abs(a[:, j]) > 1e-12, rel[:, j], 0)))} for j in range(len(h))}
    out[name] = {"shape": list(a.shape), "max_abs_diff_all": float(np.nanmax(d)), "columns": cols}
    print(f"{name}: shape {a.shape}, max|diff| = {np.nanmax(d):.3e}")
for name in ("states_baseline.npz", "states_treatment.npz"):
    A = np.load(os.path.join(REF, name))["Z"]; B = np.load(os.path.join(MINE, name))["Z"]
    d = np.abs(A - B); out[name] = {"shape": list(A.shape), "max_abs_diff": float(d.max()), "max_abs_diff_step": int(np.unravel_index(d.argmax(), d.shape)[0]),
                                    "max_abs_diff_first_1000": float(d[:1000].max()), "max_abs_diff_last_1000": float(d[-1000:].max())}
    print(f"{name}: max|ΔZ| = {d.max():.3e} (at step {out[name]['max_abs_diff_step']}), first 1000 steps {d[:1000].max():.1e}, last 1000 {d[-1000:].max():.1e}")
sa = json.load(open(os.path.join(REF, "summary.json"))); sb = json.load(open(os.path.join(MINE, "summary.json")))
keys = ["parent_residual", "parent_sigma"]; cmp = {k: [sa[k], sb[k]] for k in keys}
for br in ("baseline", "treatment"):
    for k, v in sa[br].items():
        if k == "growth_fit": continue
        cmp[f"{br}.{k}"] = [v, sb[br][k]]
    cmp[f"{br}.growth_fit.slope"] = [sa[br]["growth_fit"]["slope_ln_Hperp_per_step"], sb[br]["growth_fit"]["slope_ln_Hperp_per_step"]]
out["summary_compare"] = cmp
mism = {k: v for k, v in cmp.items() if isinstance(v[0], (int, float)) and v[0] is not None and abs(v[0] - v[1]) > 1e-9 * max(1, abs(v[0]))}
out["summary_mismatch_gt_1e-9"] = mism; print("summary keys differing by >1e-9 relative:", mism)
json.dump(out, open(os.path.join(R, "compare_with_reference.json"), "w"), indent=1)
