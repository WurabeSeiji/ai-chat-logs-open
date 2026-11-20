import numpy as np
f=np.zeros((8,8,8))
for (i,j,k) in [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]:
    for (a,b,c),s in [((i,j,k),1),((j,k,i),1),((k,i,j),1),((j,i,k),-1),((i,k,j),-1),((k,j,i),-1)]: f[a,b,c]=s
def mul(x,y):
    z=np.zeros(8); z[0]=x[0]*y[0]-x[1:]@y[1:]; z[1:]+=x[0]*y[1:]+y[0]*x[1:]
    z[1:]+=np.einsum('ijk,i,j->k',f[1:,1:,1:],x[1:],y[1:]); return z
conj=lambda z: np.concatenate([[z[0]],-z[1:]])
def ex(n,ang):  # e^{n*ang}, n: unit imaginary octonion (7-vector), ang in rad
    z=np.zeros(8); z[0]=np.cos(ang); z[1:]=n*np.sin(ang); return z
M=100; phi=20.0
psi=np.arange(0,360,0.05)
def pulse(n,theta,cm):   # 倍音部分のみ: sum_m c_m e^{n m pi (psi-theta)/180}
    Q=np.zeros((len(psi),8))
    for m in range(1,M+1):
        ang=m*np.pi*(psi-theta)/180
        Q[:,0]+=cm[m-1]*np.cos(ang); Q[:,1:]+=cm[m-1]*np.outer(np.sin(ang),n)
    return Q
def overlap(Qi,Qj):  # 交差項 Q_i conj(Q_j) のスカラー部を円周で平均
    return np.mean(np.einsum('pk,pk->p',Qi,Qj))   # Re(x conj y) = x.y
rng=np.random.default_rng(2)
def unit7():
    v=rng.normal(size=7); return v/np.linalg.norm(v)
ni,nj=unit7(),unit7()
m=np.arange(1,M+1)
box=(2/np.pi)*np.sin(m*np.pi*phi/360)/m
mu=6.0; lor=1/(m**2+mu**2)
print('軸の内積 n_i.n_j =',round(float(ni@nj),4),' (1+n_i.n_j)/2 =',round(float((1+ni@nj)/2),4))
for name,cm in [('窓(幅20度)',box),('重み 1/(m^2+mu^2), mu=6',lor)]:
    print('---',name)
    Qi=pulse(ni,0.0,cm); self_same=overlap(Qi,Qi)
    for d in [0,5,10,15,20,25,30,45,60,90]:
        same=overlap(Qi,pulse(ni,float(d),cm))/self_same
        diff=overlap(Qi,pulse(nj,float(d),cm))/self_same
        extra=''
        if name.startswith('重み') and same>0: extra='  -ln(比)/mu [度換算] = %.2f'%(-np.log(same)/mu*180/np.pi)
        print('中心の差 %3d度: 同じ軸の比 = %+.4f  違う軸の比/同じ軸の比 = %.4f%s'%(d,same,diff/same if abs(same)>1e-6 else float('nan'),extra))
