# -*- coding: utf-8 -*-
"""figures/N5_amplitude_std_compare.png（std(|Z_m|) の 10⁻¹² スケールの振動）の原因調査。入力: ../data/states_treatment.npz（力学は再実行しない）。
 1) std(|Z_m|) 系列の極値から周期を実測
 2) 親 v = Z[0] のまわりの偏差 δ(t) = Z(t) − e^{iθ_t} v（大域位相を除去）を、K_amp(v) の固有ベクトル基底で分解し、各成分の大きさと位相回転率を実測
 3) 線形回転の予言：成分 k は回転系で角速度 ANGLE·(w_k − w_0) → 周期 L/|w_k − w_0| step。実測の周期・各成分の位相回転率と比較
 4) 保存量（Σ|z|², H_int）の推移
出力: results.json, analysis.md, figures"""
import os, json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__)); Z=np.load(os.path.join(H,"..","data","states_treatment.npz"))["Z"]; T,M=Z.shape; N=16; L=1000; ANGLE=2*np.pi/L
ea,eb=np.triu_indices(N,k=1); A=np.zeros((M,M))
for i in range(M):
    sh=(ea==ea[i])|(ea==eb[i])|(eb==ea[i])|(eb==eb[i]); sh[i]=False; A[i,sh]=1
def Kamp(z): K=A*np.imag(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(K,0); return K
v=Z[0]; K0=Kamp(v); w,U=np.linalg.eigh(1j*K0); out={}
# 親は w[0]=−σ_max のモードか
ov0=abs(np.vdot(U[:,0],v))/np.linalg.norm(v); out["parent_overlap_with_w0_mode"]=float(ov0); out["w0"]=float(w[0]); out["sigma_max"]=float(-w[0])
# 1) std 系列と周期
amp=np.abs(Z); s=amp.std(axis=1); out["std_t0"]=float(s[0]); out["std_range"]=[float(s.min()),float(s.max())]
from scipy.signal import argrelextrema
pk=argrelextrema(s,np.greater,order=200)[0]; tr=argrelextrema(s,np.less,order=200)[0]
out["std_peaks_steps"]=pk.tolist(); out["std_troughs_steps"]=tr.tolist()
ext=np.sort(np.r_[pk,tr]); out["half_periods_steps"]=np.diff(ext).tolist(); out["period_estimate_steps"]=float(2*np.median(np.diff(ext))) if len(ext)>1 else None
# 2) 偏差の固有基底分解
th=np.angle(np.einsum("j,tj->t",np.conj(v),Z)); delta=Z-np.exp(1j*th)[:,None]*v[None,:]
c=delta@np.conj(U)   # c[t,k] = <U_k, δ(t)>
mag0=np.abs(c[0]); order=np.argsort(mag0)[::-1]
out["delta_norm_t0"]=float(np.linalg.norm(delta[0])); out["delta_norm_final"]=float(np.linalg.norm(delta[-1]))
top=[]
for k in order[:8]:
    ck=c[:,k]; ph=np.unwrap(np.angle(ck)); rate=np.polyfit(np.arange(T),ph,1)[0] if np.abs(ck).min()>0 else np.nan
    top.append({"k":int(k),"w_k":float(w[k]),"|c_k(0)|":float(mag0[k]),"|c_k(final)|":float(abs(ck[-1])),"measured_phase_rate_per_step":float(rate),"predicted_rate_ANGLE*(w_k-w0)":float(ANGLE*(w[k]-w[0])),"predicted_period_steps_L/|w_k-w0|":float(L/abs(w[k]-w[0])) if abs(w[k]-w[0])>1e-12 else None})
out["top_components_of_delta"]=top
# 非等モジュラー成分（振幅ずれ）が乗っているモード：Re(v̄_m δ_m) の寄与
amp_dev0=amp[0]-amp[0].mean(); out["amp_dev_t0_max"]=float(np.abs(amp_dev0).max())
# 3) 振幅分散の周波数：各辺の |Z_m| 偏差系列を FFT し、支配周波数（step⁻¹）を出す
dev=amp-amp.mean(axis=1,keepdims=True); dev=dev-dev.mean(axis=0); F=np.abs(np.fft.rfft(dev,axis=0)).sum(axis=1); freqs=np.fft.rfftfreq(T,1.0); j=np.argsort(F[1:])[::-1][:5]+1
out["amp_dev_dominant_periods_steps"]=[float(1/freqs[i]) for i in j]; out["amp_dev_dominant_power"]=[float(F[i]) for i in j]
# 予言周期の一覧（w_k − w0 の全組）
pred=sorted(set(round(float(L/abs(wk-w[0])),1) for wk in w if abs(wk-w[0])>1e-9)); out["predicted_periods_all_modes_steps"]=pred[:15]
# 4) 保存量
Hint=np.array([0.5*np.sum(np.triu(Kamp(Z[t]))**2) for t in range(0,T,50)]); out["Htotal_drift"]=float(np.ptp(np.sum(amp**2,axis=1))); out["Hint_rel_drift"]=float((Hint.max()-Hint.min())/Hint[0])
json.dump(out,open(os.path.join(H,"results.json"),"w"),indent=1)
# 図
fig,ax=plt.subplots(2,2,figsize=(14,9))
ax[0,0].plot(s); ax[0,0].set_title("std(|Z_m|) (input figure)"); ax[0,0].set_xlabel("step")
for kk in top[:4]: ax[0,1].semilogy(np.abs(c[:,kk["k"]]),label=f"k={kk['k']} w={kk['w_k']:.3f}")
ax[0,1].set_title("|<U_k, δ(t)>| of dominant components"); ax[0,1].legend(); ax[0,1].set_xlabel("step")
ax[1,0].semilogy(np.linalg.norm(delta,axis=1)); ax[1,0].set_title("‖δ(t)‖ (global phase removed)"); ax[1,0].set_xlabel("step")
ax[1,1].plot(1/freqs[1:200],F[1:200]); ax[1,1].set_xscale("log"); ax[1,1].set_title("power of |Z_m| deviations vs period (steps)"); ax[1,1].set_xlabel("period [steps]")
fig.tight_layout(); fig.savefig(os.path.join(H,"oscillation_analysis.png"),dpi=140)
print(json.dumps(out,indent=1,ensure_ascii=False))
# 5) 偏差の成長源の検算：親の固有モード残差 r（summary.json の parent_residual）が毎 step ANGLE·r·‖v‖ の偏差を生み、線形に蓄積するという仮説
s_json=json.load(open(os.path.join(H,"..","data","summary.json"))); r=s_json["parent_residual"]; nv=float(np.linalg.norm(v))
dn=np.linalg.norm(delta,axis=1); t=np.arange(T); m=t>100; slope=np.polyfit(t[m],dn[m],1)[0]
extra={"parent_residual_r":r,"norm_v":nv,"predicted_drift_per_step_ANGLE*r*|v|":ANGLE*r*nv,"measured_d|delta|/dt_per_step":float(slope),"ratio_measured/predicted":float(slope/(ANGLE*r*nv)),
       "kernel_dim_of_Kamp(v)":int(np.sum(np.abs(w)<1e-9)),"fraction_of_|delta(final)|^2_in_kernel":float(np.sum(np.abs(c[-1][np.abs(w)<1e-9])**2)/np.sum(np.abs(c[-1])**2)),
       "predicted_period_kernel_modes_L/sigma_max":float(L/(-w[0]))}
out.update(extra); json.dump(out,open(os.path.join(H,"results.json"),"w"),indent=1); print(json.dumps(extra,indent=1))
# 6) 訂正：|c_k(0)| は丸め誤差（2e-19）の分解で意味が薄い。δ(final) の分解で支配モードを特定し、その w_k と位相回転率・周期を出す
magF=np.abs(c[-1]); orderF=np.argsort(magF)[::-1]; topF=[]
for k in orderF[:10]:
    ck=c[:,k]; seg=slice(1000,T); ph=np.unwrap(np.angle(ck[seg])); rate=np.polyfit(np.arange(1000,T),ph,1)[0]
    topF.append({"k":int(k),"w_k":float(w[k]),"|c_k(final)|":float(magF[k]),"share_of_|delta(final)|^2":float(magF[k]**2/np.sum(magF**2)),"measured_phase_rate_per_step(t>1000)":float(rate),"predicted_rate_ANGLE*(w_k-w0)":float(ANGLE*(w[k]-w[0])),"predicted_period_steps":float(L/abs(w[k]-w[0])) if abs(w[k]-w[0])>1e-12 else None})
# 固有値の群ごとの |δ(final)|² の割合
groups={}
for k in range(M):
    key=round(float(w[k]),3); groups[key]=groups.get(key,0.0)+float(magF[k]**2)
tot=float(np.sum(magF**2)); share={str(k):v/tot for k,v in sorted(groups.items(),key=lambda kv:-kv[1])[:6]}
out["delta_final_top_components"]=topF; out["delta_final_share_by_eigenvalue"]=share; out["eigenvalue_spectrum_distinct"]=sorted(set(round(float(x),4) for x in w))
json.dump(out,open(os.path.join(H,"results.json"),"w"),indent=1)
print("δ(final) の固有値群別割合:",share); print("スペクトル:",out["eigenvalue_spectrum_distinct"])
for t_ in topF[:5]: print(t_)
