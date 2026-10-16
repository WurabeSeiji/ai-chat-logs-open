import itertools, math, sys
import numpy as np
from collections import Counter
def shell(R):
    C=np.array((0.5,)*4); R2=R*R
    inside=lambda p: float(np.sum((np.array(p)-C)**2))<=R2+1e-12
    A=[];B=[]; rng=range(-R-1,R+1)
    for a in itertools.product(rng,repeat=4):
        c=tuple(x+0.5 for x in a); v=[tuple(ci+s for ci,s in zip(c,sg)) for sg in itertools.product((-0.5,0.5),repeat=4) if round(sum(ci+s for ci,s in zip(c,sg)))%2==0]
        if all(inside(p) for p in v): A.append(v)
    for c in itertools.product(range(-R-1,R+2),repeat=4):
        if sum(c)%2==1:
            v=[tuple(x+(s if k==i else 0) for k,x in enumerate(c)) for i in range(4) for s in (1,-1)]
            if all(inside(p) for p in v): B.append(v)
    cubes=sum(1 for a in itertools.product(rng,repeat=4) if all(inside(tuple(x+s for x,s in zip(a,sg))) for sg in itertools.product((0,1),repeat=4)))
    fc=Counter(); owner={}
    for v in A+B:
        vv=sorted(v); cen=tuple(sum(p[k] for p in vv)/8 for k in range(4)); pairs=[];used=set()
        for p in vv:
            if p in used: continue
            q=tuple(2*c-x for c,x in zip(cen,p)); pairs.append((p,q)); used|={p,q}
        for ch in itertools.product(*pairs):
            f=frozenset(ch); fc[f]+=1; owner[f]=cen
    bnd=[f for f,n in fc.items() if n==1]
    normals=set()
    for f in bnd:
        P=np.array(sorted(f)); n=np.linalg.svd(P[1:]-P[0])[2][-1]
        if np.dot(n,P.mean(0)-np.array(owner[f]))<0: n=-n
        normals.add(tuple(np.round(n,6)))
    bv={p for f in bnd for p in f}; allv={p for v in A+B for p in v}
    ec=Counter(frozenset(e) for f in bnd for e in itertools.combinations(sorted(f),2))
    print(f'R={R}: cubes {cubes} | inscribed 16-cells {len(A)} + odd-centred {len(B)} = {len(A)+len(B)} | boundary tetra {len(bnd)}, triangles {len(bnd)*2}, vertices {len(bv)}, interior vertices {len(allv-bv)} | normal directions {len(normals)} | tets per boundary edge {sorted(Counter(ec.values()).items())}', flush=True)
for R in map(int, sys.argv[1:]): shell(R)
