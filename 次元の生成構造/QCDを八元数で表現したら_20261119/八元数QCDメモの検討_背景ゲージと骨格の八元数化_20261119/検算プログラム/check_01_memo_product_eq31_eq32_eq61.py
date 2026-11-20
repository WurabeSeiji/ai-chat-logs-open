import numpy as np, itertools
f=np.zeros((8,8,8))
for (i,j,k) in [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]:
    for (a,b,c),s in [((i,j,k),1),((j,k,i),1),((k,i,j),1),((j,i,k),-1),((i,k,j),-1),((k,j,i),-1)]:
        f[a,b,c]=s
def mul(x,y):
    z=np.zeros(8); z[0]=x[0]*y[0]-x[1:]@y[1:]
    z[1:]+=x[0]*y[1:]+y[0]*x[1:]
    for i in range(1,8):
        for j in range(1,8):
            for k in range(1,8):
                z[k]+=f[i,j,k]*x[i]*y[j]
    return z
E=np.eye(8)
# alternativity check
rng=np.random.default_rng(0)
x,y=rng.normal(size=8),rng.normal(size=8)
print('alt',np.abs(mul(mul(x,x),y)-mul(x,mul(x,y))).max(), 'norm',abs(np.linalg.norm(mul(x,y))-np.linalg.norm(x)*np.linalg.norm(y)))
# complex structure: left mult by e7
def emb(a,psi):  # a complex, psi complex 3-vector -> octonion
    z=a.real*E[0]+a.imag*E[7]
    for n in range(3):
        z=z+psi[n].real*E[n+1]+psi[n].imag*mul(E[7],E[n+1])
    return z
def dec(z):
    a=z[0]+1j*z[7]
    B=np.array([E[n+1] for n in range(3)]+[mul(E[7],E[n+1]) for n in range(3)]).T
    full=np.column_stack([E[0],E[7],B]); c=np.linalg.solve(full,z)
    return c[0]+1j*c[1], c[2:5]+1j*c[5:8]
for n in range(3): print('e7*e%d ='%(n+1), mul(E[7],E[n+1]))
a,b=rng.normal()+1j*rng.normal(),rng.normal()+1j*rng.normal()
ps=rng.normal(size=3)+1j*rng.normal(size=3); ch=rng.normal(size=3)+1j*rng.normal(size=3)
A,P=dec(mul(emb(a,ps),emb(b,ch)))
print('C part',A,'| ab-chi^dag psi',a*b-np.vdot(ch,ps),'| ab-psi^dag chi',a*b-np.vdot(ps,ch))
cands={'a chi+bbar psi+cross(psibar,chibar)':a*ch+np.conj(b)*ps+np.cross(np.conj(ps),np.conj(ch)),
'a chi+bbar psi-cross':a*ch+np.conj(b)*ps-np.cross(np.conj(ps),np.conj(ch)),
'a chi+b psi+cross':a*ch+b*ps+np.cross(np.conj(ps),np.conj(ch)),
'abar chi+b psi+cross':np.conj(a)*ch+b*ps+np.cross(np.conj(ps),np.conj(ch)),
'abar chi+b psi-cross':np.conj(a)*ch+b*ps-np.cross(np.conj(ps),np.conj(ch))}
print('C3 part',P)
for k,v in cands.items(): print(k, np.abs(v-P).max())
# (3.2): psi^dag chi = [chi psi*]_C
conj=lambda z: np.concatenate([[z[0]],-z[1:]])
A2,_=dec(mul(emb(0,ch),conj(emb(0,ps)))); print('[chi psi*]_C',A2,'psi^dag chi',np.vdot(ps,ch))
# (6.1)
t,R=rng.normal(),rng.normal(); xv=rng.normal(size=3); Q=rng.normal(size=3)
q1=t*E[0]+xv@E[1:4]; q2=R*E[0]+Q@E[1:4]
for lab,z in [('q1+e7*q2',q1+mul(E[7],q2)),('q1+q2*e7',q1+mul(q2,E[7]))]:
    A3,P3=dec(z); print(lab,'C:',A3,'vs t+iR',t+1j*R,'| C3-(x-iQ):',np.abs(P3-(xv-1j*Q)).max(),' C3-(x+iQ):',np.abs(P3-(xv+1j*Q)).max())
