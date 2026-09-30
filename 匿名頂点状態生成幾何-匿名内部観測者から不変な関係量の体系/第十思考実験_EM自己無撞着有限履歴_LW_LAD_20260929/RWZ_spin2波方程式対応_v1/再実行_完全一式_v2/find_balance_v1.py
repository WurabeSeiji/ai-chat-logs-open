import json, sys
sys.path.insert(0,'/mnt/data/rwz_spin2_v1')
from rwz_spin2_flux_v1 import run

def diff(parity,x0):
    r=run(parity=parity,x0=x0,tmax=360,nx=5201,sample_stride=500)
    return r['fraction_infinity']-r['fraction_horizon'],r

out={}
for parity in ['odd','even']:
    a,b=-10.0,10.0
    fa,ra=diff(parity,a); fb,rb=diff(parity,b)
    # expand if needed
    if fa*fb>0:
        a,b=-30.0,30.0; fa,ra=diff(parity,a); fb,rb=diff(parity,b)
    for _ in range(12):
        m=(a+b)/2
        fm,rm=diff(parity,m)
        if fa*fm<=0:
            b,fb,rb=m,fm,rm
        else:
            a,fa,ra=m,fm,rm
    x=(a+b)/2
    f,r=diff(parity,x)
    out[parity]={'x0_balance':x,'difference_inf_minus_hor':f,'result':{k:r[k] for k in ['fraction_horizon','fraction_infinity','fraction_remaining','closure','E0','E_horizon','E_infinity']}}
    print(parity,out[parity])
with open('/mnt/data/rwz_spin2_v1/rwz_balance_search_v1.json','w') as fh: json.dump(out,fh,indent=2)
