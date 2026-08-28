import pandas as pd,numpy as np,itertools,math,sys,time
from scipy.spatial import cKDTree
N=int(sys.argv[1]); out=sys.argv[2]
df=pd.read_csv(f'/mnt/data/N3_N16_partial_zero_closure_analysis_20260826/SOURCE_N{N}_step5000_final_edges.csv')
z=df.z2_re.to_numpy()+1j*df.z2_im.to_numpy(); n=len(z); absz=np.abs(z)
def comb_array(n,k):
 return np.fromiter((x for c in itertools.combinations(range(n),k) for x in c),dtype=np.int16,count=math.comb(n,k)*k).reshape(-1,k)
def best(k,neighbors=16,batch=25000):
 p=k//2;q=k-p;A=comb_array(n,p);B=comb_array(n,q);SB=z[B].sum(axis=1);tree=cKDTree(np.c_[SB.real,SB.imag]);br=1e99;bi=None
 for s in range(0,len(A),batch):
  X=A[s:s+batch];SA=z[X].sum(axis=1);_,ids=tree.query(np.c_[-SA.real,-SA.imag],k=min(neighbors,len(B)))
  if ids.ndim==1:ids=ids[:,None]
  for row,a in enumerate(X):
   aset=set(map(int,a))
   for bj in ids[row]:
    b=B[int(bj)]
    if any(int(x) in aset for x in b):continue
    idx=tuple(sorted(aset|set(map(int,b))))
    if len(idx)!=k:continue
    rr=abs(z[list(idx)].sum())/absz[list(idx)].sum()
    if rr<br:br=rr;bi=idx
 return br,bi
rows=[]
for k in [5,6]:
 t=time.time();r,idx=best(k);edges=';'.join(f'{int(df.iloc[e].i)}-{int(df.iloc[e].j)}' for e in idx);rows.append([N,len(df),k,r,edges,time.time()-t]);print(N,k,r,edges)
pd.DataFrame(rows,columns=['N','M','k','best_residual','edges','runtime_s']).to_csv(out,index=False)
