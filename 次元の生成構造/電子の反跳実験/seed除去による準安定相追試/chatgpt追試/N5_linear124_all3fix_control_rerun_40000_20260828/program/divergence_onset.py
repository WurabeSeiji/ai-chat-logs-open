# -*- coding: utf-8 -*-
"""Where does this rerun diverge from the reference (round-off amplification)? Output: results/divergence_onset.json"""
import csv, json, os, numpy as np
HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); REF=os.path.join(os.path.dirname(HERE),"N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828","data")
out={}
for name in ("treatment_linear124_amplitude_aware_timeseries.csv","baseline_linear124_phase_only_timeseries.csv"):
    a=np.array([[float(x) for x in r] for r in list(csv.reader(open(os.path.join(REF,name))))[1:]]); b=np.array([[float(x) for x in r] for r in list(csv.reader(open(os.path.join(HERE,"data",name))))[1:]])
    d=np.abs(a[:,2]-b[:,2]); rec={}
    for thr in (1e-8,1e-6,1e-4,1e-2):
        ix=np.where(d>thr)[0]; rec[f"first_step_|dHperp|>{thr:.0e}"]=int(ix[0]) if len(ix) else None
    out[name.split("_")[0]]=rec; print(name.split("_")[0], rec)
A=np.load(os.path.join(REF,"states_treatment.npz"))["Z"]; B=np.load(os.path.join(HERE,"data","states_treatment.npz"))["Z"]
g=np.vdot(A[0],B[0]); out["treatment_step0_gauge_removed_state_distance"]=float(np.linalg.norm(A[0]-B[0]*np.exp(-1j*np.angle(g)))/np.linalg.norm(A[0]))
out["treatment_step0_abs_amplitude_maxdiff_sorted"]=float(np.abs(np.sort(np.abs(A[0]))-np.sort(np.abs(B[0]))).max())
print("treatment step 0: gauge-removed state distance %.3f (parents differ by a symmetry: sorted-amplitude max diff %.1e)"%(out["treatment_step0_gauge_removed_state_distance"],out["treatment_step0_abs_amplitude_maxdiff_sorted"]))
json.dump(out,open(os.path.join(HERE,"results","divergence_onset.json"),"w"),indent=1)
