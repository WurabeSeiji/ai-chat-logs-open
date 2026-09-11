import numpy as np
from pathlib import Path

src=Path('/mnt/data/run_N3_N40_stage123_v1.py').read_text(encoding='utf-8')
g={'__name__':'canonical','__file__':'/mnt/data/run_N3_N40_stage123_v1.py'}
exec(compile(src.split('rows=[]; summaries=[]')[0],g['__file__'],'exec'),g)
adjacency=g['adjacency']

def integer_K(labels,A):
    diff=(labels[None,:]-labels[:,None])%4
    s=np.zeros_like(diff,dtype=float)
    s[diff==1]=1.0; s[diff==3]=-1.0
    return A*s

def floor(N,max_iter=500):
    ea,eb=np.triu_indices(N,1); A=adjacency(N); labels=((ea+eb)%4).astype(int); seen={}
    for _ in range(max_iter):
        K=integer_K(labels,A); w,V=np.linalg.eigh(1j*K); v=V[:,-1]
        new=(np.round(np.degrees(np.angle(v))%360/90).astype(int))%4
        cands=[((new+sh)%4) for sh in range(4)]
        cand=max(cands,key=lambda c:int(np.sum(c==labels)))
        key=tuple(cand.tolist())
        if np.array_equal(cand,labels) or key in seen: return labels,K
        seen[key]=1; labels=cand
    return labels,integer_K(labels,A)

worst_factor=0.0; worst_spec=0.0
for N in range(3,41):
    labels,K=floor(N); ea,eb=np.triu_indices(N,1); M=len(ea)
    B=np.zeros((N,M)); B[ea,np.arange(M)]=1; B[eb,np.arange(M)]=1
    th=labels*np.pi/2; c=np.cos(th); s=np.sin(th)
    X=np.vstack([B*c,B*s]); J=np.block([[np.zeros((N,N)),np.eye(N)],[-np.eye(N),np.zeros((N,N))]])
    K2=X.T@J@X
    worst_factor=max(worst_factor,float(np.max(np.abs(K-K2))))

    # predicted positive squared spectrum (N>=5 generic; small N compare only factorization)
    if N>=5:
        p=(N+1)//2 if N%2 else N//2
        q=N//2
        pred=[]
        if q-1>0 and p*(q-2)>1e-12: pred += [p*(q-2)]*(q-1)
        if p-1>0 and q*(p-2)>1e-12: pred += [q*(p-2)]*(p-1)
        glob=4*p*q-2*N
        if glob>1e-12: pred += [glob]
        obs=np.linalg.eigvalsh(1j*K)
        obs2=sorted([float(x*x) for x in obs if x>1e-8])
        pred=sorted(map(float,pred))
        if len(obs2)!=len(pred):
            raise RuntimeError((N,len(obs2),len(pred)))
        worst_spec=max(worst_spec,max(abs(a-b) for a,b in zip(obs2,pred)))

print('worst |K-X^TJX| =',worst_factor)
print('worst positive omega^2 formula error N=5..40 =',worst_spec)
