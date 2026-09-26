#!/usr/bin/env python3
import math, json, importlib.util, sys
import numpy as np
spec=importlib.util.spec_from_file_location('base','/mnt/data/audit_two_paths_rk4.py')
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
H,PH,PH2,TAU0=base.H,base.PH,base.PH2,base.TAU0
rr,gN,ref_step,step_A,init_B,coreB=base.rr,base.gN,base.ref_step,base.step_A,base.init_B,base.coreB

def micro_fast(s,nu):
    j=int(np.argmax(s['q']))
    # copy old -> new first (synchronous preserve unless this phase updates field)
    o={k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in s.items()}
    n={k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in s.items()}
    U,P,E=o['U'],o['P'],o['E']
    if j==0:
        n['k1']=rr(U,P,E,nu)
    elif j==1:
        n['Us']=U*PH2; n['Ps']=P+0.5*H*o['k1'][0]; n['Es']=E+0.5*H*o['k1'][1]
    elif j==2:
        n['k2']=rr(o['Us'],o['Ps'],o['Es'],nu)
    elif j==3:
        n['Us']=U*PH2; n['Ps']=P+0.5*H*o['k2'][0]; n['Es']=E+0.5*H*o['k2'][1]
    elif j==4:
        n['k3']=rr(o['Us'],o['Ps'],o['Es'],nu)
    elif j==5:
        n['Us']=U*PH; n['Ps']=P+H*o['k3'][0]; n['Es']=E+H*o['k3'][1]
    elif j==6:
        n['k4']=rr(o['Us'],o['Ps'],o['Es'],nu)
    elif j==7:
        n['d']=(H/6)*(o['k1']+2*o['k2']+2*o['k3']+o['k4'])
    elif j==8:
        n['Um']=U*PH2; n['Pm']=P+0.5*o['d'][0]; n['Em']=E+0.5*o['d'][1]
    elif j==9:
        n['dphi']=H*gN(o['Um'],o['Pm'],o['Em'])
    elif j==10:
        n['U']=U*PH; n['P']=P+o['d'][0]; n['E']=E+o['d'][1]
        n['H']=o['H']*complex(math.cos(o['dphi']),math.sin(o['dphi']))
        n['Q']=o['Q']*math.exp(o['d'][2]/TAU0)
        # reset work to explicit safe state of the just-committed macro state
        n['k1']=np.zeros(3); n['k2']=np.zeros(3); n['k3']=np.zeros(3); n['k4']=np.zeros(3)
        n['Us']=U*PH; n['Ps']=P+o['d'][0]; n['Es']=E+o['d'][1]; n['d']=np.zeros(3)
        n['Um']=U*PH; n['Pm']=P+o['d'][0]; n['Em']=E+o['d'][1]; n['dphi']=0.0
    n['q']=np.roll(o['q'],1)
    return n

def macro_fast(s,nu):
    for _ in range(11): s=micro_fast(s,nu)
    return s

def compare_case(p0,e0,qmass,nsteps):
    nu=qmass/(1+qmass)**2
    x0=np.array([1+0j,p0,e0,1+0j,1+0j],complex)
    r=x0.copy(); a=x0.copy(); b=init_B(*x0)
    ma=np.zeros(5); mb=np.zeros(5)
    for _ in range(nsteps):
        r=ref_step(r,nu); a=step_A(a,nu); b=macro_fast(b,nu)
        ma=np.maximum(ma,np.abs(a-r)); mb=np.maximum(mb,np.abs(coreB(b)-r))
    return ma,mb

def csv_check():
    data=np.genfromtxt('/mnt/data/experiment02_sn_2p5pn_rr_C1_N12.csv',delimiter=',',names=True)
    p0=float(data['p'][0]); e0=float(data['e'][0]); nu=.25
    a=np.array([1+0j,p0,e0,1+0j,1+0j],complex); b=init_B(*a)
    err={k:0.0 for k in ['A_p','A_e','A_t','A_xy','B_p','B_e','B_t','B_xy']}
    for i,row in enumerate(data):
        for label,z in [('A',a),('B',coreB(b))]:
            U,P,E,Hh,Q=z; ce=base.trigU(U)[0]; r=P.real/(1+E.real*ce)
            xx=r*Hh.real; yy=r*Hh.imag; tt=TAU0*math.log(Q.real)
            err[label+'_p']=max(err[label+'_p'],abs(P.real-row['p']))
            err[label+'_e']=max(err[label+'_e'],abs(E.real-row['e']))
            err[label+'_t']=max(err[label+'_t'],abs(tt-row['t_osc']))
            err[label+'_xy']=max(err[label+'_xy'],abs(xx-row['x']),abs(yy-row['y']))
        if i+1<len(data): a=step_A(a,nu); b=macro_fast(b,nu)
    return err

if __name__=='__main__':
    tests=[(24.3,0.4358898943540673,1.0,12000),(60,.45,4.0,12000),(120,.55,1.0,12000)]
    rows=[]; ga=np.zeros(5); gb=np.zeros(5)
    for t in tests:
        ma,mb=compare_case(*t); ga=np.maximum(ga,ma); gb=np.maximum(gb,mb); print(t,'A',ma,'B',mb); rows.append([t,ma.tolist(),mb.tolist()])
    cc=csv_check(); print('GLOBAL_A',ga); print('GLOBAL_B',gb); print('CSV',cc)
    with open('/mnt/data/two_path_rk4_fast_audit.json','w') as f: json.dump({'tests':rows,'global_A':ga.tolist(),'global_B':gb.tolist(),'csv':cc},f,indent=2)
