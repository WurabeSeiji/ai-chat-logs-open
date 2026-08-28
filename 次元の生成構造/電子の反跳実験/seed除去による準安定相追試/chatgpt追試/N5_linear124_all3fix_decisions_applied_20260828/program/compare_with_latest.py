# -*- coding: utf-8 -*-
"""判断1〜6 適用版と、最新版（無変更）をこの計算機で走らせた対照 (`../N5_linear124_all3fix_control_rerun_40000_20260828`) を突合。
 力学（H∥,H⊥,H_total,PR,|ZᵀZ|,振幅,K ノルム）と全状態 npz は一致すべき（読出しだけが変わる）。σ₁ 列は読出し器が違うので差が出る。
 出力: results/compare_with_latest.json"""
import os, json, numpy as np, pandas as pd
H=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); C=os.path.join(os.path.dirname(H),"N5_linear124_all3fix_control_rerun_40000_20260828","data"); R=os.path.join(os.path.dirname(H),"N5_linear124_all3fix_seedless_parentnorm_removed_40000_20260828","data")
out={}
for br in ("baseline_linear124_phase_only","treatment_linear124_amplitude_aware"):
    a=pd.read_csv(os.path.join(C,br+"_timeseries.csv")); b=pd.read_csv(os.path.join(H,"data",br+"_timeseries.csv")); r={}
    same=[c for c in a.columns if c in b.columns and c!="sigma1_phase_reader"]
    r["dynamics_columns_max_abs_diff"]={c:float(np.nanmax(np.abs(a[c].to_numpy()-b[c].to_numpy()))) for c in same}
    r["dynamics_bit_identical"]=all(v==0.0 for v in r["dynamics_columns_max_abs_diff"].values())
    za=np.load(os.path.join(C,"states_"+br.split("_")[0]+".npz"))["Z"]; zb=np.load(os.path.join(H,"data","states_"+br.split("_")[0]+".npz"))["Z"]
    r["states_max_abs_diff"]=float(np.abs(za-zb).max())
    # 読出しの比較：旧 σ₁（位相のみ読出し器）と新 σ₁（実際の生成子）
    s_old=a["sigma1_phase_reader"].to_numpy(); s_new=b["sigma1_of_step_generator"].to_numpy()
    r["sigma1_old_phase_reader"]={"step1":float(s_old[1]),"final":float(s_old[-1]),"min":float(np.nanmin(s_old)),"max":float(np.nanmax(s_old))}
    r["sigma1_new_actual_generator"]={"step1":float(s_new[1]),"final":float(s_new[-1]),"min":float(np.nanmin(s_new)),"max":float(np.nanmax(s_new))}
    # 位相進みの検算：実測 arg<Z_t,Z_{t+1}> と ANGLE·σ₁(K_t)。K_t は t→t+1 を回した生成子（csv では t+1 行に格納）
    mp=b["measured_phase_advance"].to_numpy(); asig=b["angle_times_sigma1_next"].to_numpy(); m=np.isfinite(mp)&np.isfinite(asig)
    r["phase_advance_check"]={"step0_measured":float(mp[0]),"step0_angle_sigma1":float(asig[0]),"step0_rel_diff":float(abs(mp[0]-asig[0])/asig[0]),
        "max_abs_diff_first_50_steps":float(np.abs(mp[:50]-asig[:50]).max()),"mean_ratio_measured_over_angle_sigma_all":float(np.mean(mp[m]/asig[m]))}
    out[br]=r; print(br); print(json.dumps(r,indent=1))
# 最新版の保存データ（ChatGPT 計算機）との差は既知の丸め発散（対照実験で確認済）。参考として H⊥ の一致区間だけ記録
a=pd.read_csv(os.path.join(R,"treatment_linear124_amplitude_aware_timeseries.csv")); b=pd.read_csv(os.path.join(H,"data","treatment_linear124_amplitude_aware_timeseries.csv"))
d=np.abs(a.H_perp.to_numpy()-b.H_perp.to_numpy()); ix=np.where(d>1e-8)[0]; out["vs_saved_latest_data_first_step_|dHperp|>1e-8"]=int(ix[0]) if len(ix) else None
print("vs saved latest data (other machine): first |dH_perp|>1e-8 at step", out["vs_saved_latest_data_first_step_|dHperp|>1e-8"])
json.dump(out,open(os.path.join(H,"results","compare_with_latest.json"),"w"),indent=1)
