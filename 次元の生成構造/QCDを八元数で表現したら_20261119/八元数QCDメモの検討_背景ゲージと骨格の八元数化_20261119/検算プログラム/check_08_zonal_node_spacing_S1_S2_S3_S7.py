import numpy as np
from scipy.special import gegenbauer, eval_legendre
from numpy.polynomial import polynomial as P
k=10
chi=np.linspace(1e-6,np.pi-1e-6,2000001)
def zeros(y):
    s=np.sign(y); idx=np.where(s[:-1]*s[1:]<0)[0]
    return np.degrees(chi[idx])
cases={'S1 (円周)        cos(k chi)':np.cos(k*chi),
       'S2 (球面)        Legendre  ':eval_legendre(k,np.cos(chi)),
       'S3 (四元数の単位) lambda=1  ':gegenbauer(k,1.0)(np.cos(chi)),
       'S7 (八元数の単位) lambda=3  ':gegenbauer(k,3.0)(np.cos(chi))}
for name,y in cases.items():
    z=zeros(y); d=np.diff(z)
    print(name,'| 節の間隔: 最小 %.3f 度, 最大 %.3f 度'%(d.min(),d.max()),'| 端から最初の節まで %.3f 度'%z[0])
# cos(k chi) を S3 と S7 の固有モードで展開したときの成分数
x=np.cos(chi)
def ncomp(lam):
    A=np.array([gegenbauer(j,lam)(x) for j in range(k+1)]).T
    w=np.sin(chi)**(2*lam)        # 球面の測度
    coef=np.linalg.lstsq(A*np.sqrt(w)[:,None],np.cos(k*chi)*np.sqrt(w),rcond=None)[0]
    return [(j,round(c,4)) for j,c in enumerate(coef) if abs(c)>1e-6]
print('cos(10 chi) の S3 固有モード展開:',ncomp(1.0))
print('cos(10 chi) の S7 固有モード展開:',ncomp(3.0))
