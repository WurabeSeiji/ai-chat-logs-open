import pandas as pd,numpy as np,sys,itertools
N=int(sys.argv[1]); edgestr=sys.argv[2]
df=pd.read_csv(f'/mnt/data/N3_N16_partial_zero_closure_analysis_20260826/SOURCE_N{N}_step5000_final_edges.csv')
q=df.z2_re.to_numpy()+1j*df.z2_im.to_numpy(); w=np.abs(q); mp={(int(r.i),int(r.j)):i for i,r in df.iterrows()}
S=set(mp[tuple(map(int,e.split('-')))] for e in edgestr.split(';') if e)
def rr(S): idx=np.array(sorted(S)); return abs(q[idx].sum())/w[idx].sum()
print('start',rr(S),len(S))
for it in range(20):
 best=rr(S); bestS=None; sel=list(S); uns=[i for i in range(len(df)) if i not in S]
 # single swaps
 for o in sel:
  for x in uns:
   T=(S-{o})|{x}; r=rr(T)
   if r<best: best=r;bestS=T
 # double swaps
 for o1,o2 in itertools.combinations(sel,2):
  baseS=S-{o1,o2}
  for x1,x2 in itertools.combinations(uns,2):
   T=baseS|{x1,x2};r=rr(T)
   if r<best:best=r;bestS=T
 if bestS is None:break
 S=bestS;print('iter',it,best)
print('final',rr(S),';'.join(f'{int(df.iloc[i].i)}-{int(df.iloc[i].j)}' for i in sorted(S)))
