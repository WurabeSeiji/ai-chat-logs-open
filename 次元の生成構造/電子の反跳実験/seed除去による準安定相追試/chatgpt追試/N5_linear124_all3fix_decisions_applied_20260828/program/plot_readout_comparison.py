# -*- coding: utf-8 -*-
"""σ₁ の旧読出し（位相のみ K）と新読出し（実際の生成子）、および実測位相進みの比較図。出力 figures/N5_sigma1_readout_comparison.png"""
import os, math, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Hiragino Sans","font.size":11})
H=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); C=os.path.join(os.path.dirname(H),"N5_linear124_all3fix_control_rerun_40000_20260828","data")
fig,ax=plt.subplots(1,3,figsize=(16,4.6))
for i,(br,lab) in enumerate((("baseline_linear124_phase_only","baseline（位相のみ K）"),("treatment_linear124_amplitude_aware","treatment（振幅込み K）"))):
    a=pd.read_csv(os.path.join(C,br+"_timeseries.csv")); b=pd.read_csv(os.path.join(H,"data",br+"_timeseries.csv"))
    ax[i].plot(a.step,a.sigma1_phase_reader,color="#c0392b",lw=1,label="旧読出し：位相のみ K の σ₁")
    ax[i].plot(b.step,b.sigma1_of_step_generator,color="#1f4e79",lw=1,label="新読出し：実際に回した K の σ₁")
    ax[i].set_title(lab); ax[i].set_xlabel("step"); ax[i].set_xlim(0,5000); ax[i].legend(fontsize=9)
b=pd.read_csv(os.path.join(H,"data","treatment_linear124_amplitude_aware_timeseries.csv")); bb=pd.read_csv(os.path.join(H,"data","baseline_linear124_phase_only_timeseries.csv"))
ax[2].plot(bb.step,bb.measured_phase_advance,color="#c0392b",lw=1,label="baseline 実測 arg⟨Z_t,Z_{t+1}⟩"); ax[2].plot(bb.step,bb.angle_times_sigma1_next,color="#c0392b",lw=1,ls="--",label="baseline (2π/L)·σ₁")
ax[2].plot(b.step,b.measured_phase_advance,color="#1f4e79",lw=1,label="treatment 実測"); ax[2].plot(b.step,b.angle_times_sigma1_next,color="#1f4e79",lw=1,ls="--",label="treatment (2π/L)·σ₁")
ax[2].set_title("1 step の位相進み：実測 vs (2π/L)·σ₁"); ax[2].set_xlabel("step"); ax[2].set_xlim(0,2000); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(H,"figures","N5_sigma1_readout_comparison.png"),dpi=160); print("ok")
