#!/usr/bin/env python3
import math, json, sys
import numpy as np
sys.path.insert(0,'/mnt/data')
import experiment03_core as e3

LP,LE,LT=0.01,1.0,1.0e-5

def enc3(p,e,t=0.0):
    return np.array([math.exp(LP*p), math.exp(LE*e), math.exp(LT*t)],float)

def dec3(Z):
    return np.array([math.log(Z[0])/LP, math.log(Z[1])/LE, math.log(Z[2])/LT],float)

def exp_inc(v, scales):
    return np.exp(np.asarray(scales)*np.asarray(v))

def f(chi, Z, nu):
    # factor-generator oracle only; decoder is precisely the remaining item to eliminate later.
    p,e,_=dec3(Z)
    dp,de,dt,_=e3.rr_rates_chi(chi,float(p),float(e),nu)
    return np.array([dp,de,dt],float)

def rk4_mul(chi,p,e,hchi,nu):
    scales=np.array([LP,LE,LT],float)
    Z=enc3(p,e,0.0)
    k1=f(chi,Z,nu)
    Z2=Z*exp_inc(0.5*hchi*k1,scales)
    k2=f(chi+0.5*hchi,Z2,nu)
    Z3=Z*exp_inc(0.5*hchi*k2,scales)
    k3=f(chi+0.5*hchi,Z3,nu)
    Z4=Z*exp_inc(hchi*k3,scales)
    k4=f(chi+hchi,Z4,nu)
    # exp[h/6(k1+2k2+2k3+k4)] = product of exponential factors
    Zout=Z.copy()
    Zout*=exp_inc((hchi/6.0)*k1,scales)
    Zout*=exp_inc((hchi/3.0)*k2,scales)
    Zout*=exp_inc((hchi/3.0)*k3,scales)
    Zout*=exp_inc((hchi/6.0)*k4,scales)
    p1,e1,dt=dec3(Zout)
    return float(p1),float(e1),float(dt)

def main():
    tests=[]
    for p in (24.3,60.0,120.0):
      for e in (.2,.45,.55):
       for q in (1.,2.,4.,10.):
        nu=q/(1+q)**2
        h=2*math.pi/4000
        worst=[0.,0.,0.]
        chi=0.; pp=p; ee=e
        for _ in range(200):
            a=e3.rk4_elements(chi,pp,ee,h,nu)
            b=rk4_mul(chi,pp,ee,h,nu)
            for i in range(3): worst[i]=max(worst[i],abs(a[i]-b[i]))
            pp,ee=a[0],a[1]; chi+=h
        tests.append({'p':p,'e':e,'q':q,'max_dp':worst[0],'max_de':worst[1],'max_dt':worst[2]})
    summary={'tests':tests,'worst':{
      'p':max(x['max_dp'] for x in tests),'e':max(x['max_de'] for x in tests),'dt':max(x['max_dt'] for x in tests)}}
    print(json.dumps(summary,indent=2))
    with open('/mnt/data/multiplicative_rk4_shell_audit.json','w') as f: json.dump(summary,f,indent=2)
if __name__=='__main__': main()
