#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import json, hashlib
from pathlib import Path
import numpy as np
import unified_engine as eng
import unified_init as init

HERE=Path(__file__).resolve().parent
OUT=HERE/'results_three_odd_harmonics_v1'
PATTERN={'pattern_id':'three_odd_harmonics_v1','charge_multiple_n':1,'sign_a':1,'sign_b':-1}
MACROSTEPS=4000

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def indices():
    L=init.load_supplement(); size=L.NST; ab=size
    return L,size,ab+L.UR,ab+L.UI,ab+L.QB

def build_ext(base,kind0,arg0):
    n0=base.size; n=n0+6
    psi=np.empty(n); psi[:n0]=base; psi[n0:]=[1.,0.,1.,0.,1.,0.]
    kind=np.full((n,n),eng.K_ZERO,dtype=np.int64); arg=np.zeros((n,n,4),dtype=np.int64)
    kind[:n0,:n0]=kind0; arg[:n0,:n0]=arg0
    _,_,_,_,ab_qb=indices(); v1r,v1i,v3r,v3i,v5r,v5i=range(n0,n0+6)
    # q0..q9: keep all harmonic readouts unchanged.
    for j in range(10):
        col=ab_qb+j
        for row in (v1r,v1i,v3r,v3i,v5r,v5i):
            kind[row,col]=eng.K_IDENT; arg[row,col,0]=row
    # q10: update V1 by half angle, derive V3 and V5 from the same updated V1 by multiplication only.
    col=ab_qb+10
    kind[v1r,col]=eng.K_ROT_V1_RE; arg[v1r,col,0]=v1r; arg[v1r,col,1]=v1i
    kind[v1i,col]=eng.K_ROT_V1_IM; arg[v1i,col,0]=v1r; arg[v1i,col,1]=v1i
    for row,k in ((v3r,eng.K_DERIVE_V3_RE),(v3i,eng.K_DERIVE_V3_IM),(v5r,eng.K_DERIVE_V5_RE),(v5i,eng.K_DERIVE_V5_IM)):
        kind[row,col]=k; arg[row,col,0]=v1r; arg[row,col,1]=v1i
    return psi,kind,arg

def cpow(v,p): return v**p

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    b,names,metas,norm=init.initial_state(PATTERN); k0,a0=init.action_definition(b.size)
    e,k,a=build_ext(b,k0,a0); n0=b.size
    _,_,ur,ui,_=indices()
    bn=np.empty_like(b); en=np.empty_like(e); bs=np.empty((n0,n0)); es=np.empty((e.size,e.size))
    # compile
    bt=b.copy(); et=e.copy(); eng.microstep(bt,bn,bs,k0,a0); eng.microstep(et,en,es,k,a)
    b,_,_,_=init.initial_state(PATTERN); e,k,a=build_ext(b,k0,a0)
    maxdiff=0.; first=None; maxv1sq=0.; maxv3=0.; maxv5=0.; maxnorm=[0.,0.,0.]; checks=[]
    def audit(mac,mic):
        nonlocal maxdiff,first,maxv1sq,maxv3,maxv5,maxnorm
        d=np.abs(b-e[:n0]); maxdiff=max(maxdiff,float(d.max()))
        if first is None and not np.array_equal(b,e[:n0]): first=(mac,mic,int(np.flatnonzero(b!=e[:n0])[0]))
        v1=complex(e[n0],e[n0+1]); v3=complex(e[n0+2],e[n0+3]); v5=complex(e[n0+4],e[n0+5]); u=complex(e[ur],e[ui])
        maxv1sq=max(maxv1sq,abs(v1*v1-u)); maxv3=max(maxv3,abs(v3-v1**3)); maxv5=max(maxv5,abs(v5-v1**5))
        for i,v in enumerate((v1,v3,v5)): maxnorm[i]=max(maxnorm[i],abs(abs(v)**2-1.0))
    audit(0,0)
    for mac in range(1,MACROSTEPS+1):
        for mic in range(1,eng.MICROSTEPS_PER_MACRO+1):
            eng.microstep(b,bn,bs,k0,a0); eng.microstep(e,en,es,k,a); audit(mac,mic)
        if mac in (1,10,100,1000,2000,4000):
            checks.append({'macrostep':mac,'P_ab':float(e[init.load_supplement().NST+init.load_supplement().P]),
                           'U_ab':[float(e[ur]),float(e[ui])],
                           'V1':[float(e[n0]),float(e[n0+1])],
                           'V3':[float(e[n0+2]),float(e[n0+3])],
                           'V5':[float(e[n0+4]),float(e[n0+5])]})
    res={'experiment':'Paper 8 numerical experiment 3: three odd harmonics V,V^3,V^5, no added coupling',
         'principle':'original 123-state S unchanged; harmonics derived from V by multiplication only; no phase-to-orbit coupling added because original S has none',
         'macrosteps':MACROSTEPS,'microsteps_per_macrostep':eng.MICROSTEPS_PER_MACRO,'original_dim':n0,'extended_dim':e.size,
         'bitwise_equal_original_123':first is None,'first_mismatch':first,'max_abs_diff_original_123':maxdiff,
         'max_abs_V1_squared_minus_Uab':maxv1sq,'max_abs_V3_minus_V1_cubed':maxv3,'max_abs_V5_minus_V1_fifth':maxv5,
         'max_abs_norm2_minus_1':maxnorm,'checkpoints':checks,
         'final_P_ab':checks[-1]['P_ab'],'final_U_ab':checks[-1]['U_ab'],'final_V1':checks[-1]['V1'],'final_V3':checks[-1]['V3'],'final_V5':checks[-1]['V5']}
    (OUT/'result.json').write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    np.save(OUT/'final_original_123.npy',b); np.save(OUT/'final_extended_129.npy',e)
    md=f'''# 第八思考実験 数値実験3：3奇数倍音 V, V^3, V^5（作り込みなし）\n\n## 原則\n\n現行123状態の相互作用 S は変更しない。監査すると軌道更新 dP は P,N,C,D のみを読み、U,H等の位相を読まない。従って、倍音を軌道へ効かせる新項を追加することは行わない。\n\n追加したのは二重被覆 V と、その同一時点の V から複素乗算だけで生成する\n\n- V\n- V^3\n- V^5\n\nの読み出し状態だけである。等振幅で初期値はすべて1。q10で V を半角回転し、同じ更新後 V から V^3,V^5 を生成する。\n\n## 結果\n\n- 既存123成分 bitwise完全一致: **{first is None}**\n- 最大差: `{maxdiff:.17g}`\n- max |V^2-U_ab|: `{maxv1sq:.17g}`\n- max |V^3-(V)^3|: `{maxv3:.17g}`\n- max |V^5-(V)^5|: `{maxv5:.17g}`\n- max ||V_h|^2-1|: `{maxnorm}`\n- 4000 macrosteps後 P_ab: `{checks[-1]['P_ab']:.17g}`\n\n## 分析\n\n3奇数倍音を追加しても既存軌道は変化しなかった。これは失敗ではなく、現行 S の構造から必然である。現行の dP は位相 U を入力に持たないため、位相の二重被覆や奇数倍音を追加しても、作り込みなしでは軌道へ戻る経路が存在しない。\n\nしたがって本実験は、現在のモデルが「4π位相構造を保持できる」一方で、「その位相構造を力学へフィードバックする規則はまだ含んでいない」ことを切り分けた。次に軌道変化を検討する場合、恣意的な新しい力を加えるのではなく、既存の積和 S のどの交差項が位相から P への写像として匿名性・ステートレス性から要求されるかを導出する必要がある。\n'''
    (OUT/'EXPERIMENT_CHANGE_AND_RESULT_ANALYSIS_ja.md').write_text(md,encoding='utf-8')
    sums={p.name:sha(p) for p in [Path(__file__),HERE/'unified_engine.py',HERE/'unified_init.py',OUT/'result.json',OUT/'final_original_123.npy',OUT/'final_extended_129.npy',OUT/'EXPERIMENT_CHANGE_AND_RESULT_ANALYSIS_ja.md']}
    (OUT/'SHA256SUMS.json').write_text(json.dumps(sums,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(res,ensure_ascii=False,indent=2))
if __name__=='__main__': run()
