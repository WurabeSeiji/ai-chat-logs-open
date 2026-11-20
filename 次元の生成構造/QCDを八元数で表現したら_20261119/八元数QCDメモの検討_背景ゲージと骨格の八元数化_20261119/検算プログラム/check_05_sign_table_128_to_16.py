import numpy as np, itertools
lines=[(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]
def table(s):
    f=np.zeros((8,8,8))
    for (i,j,k),sg in zip(lines,s):
        for (a,b,c),p in [((i,j,k),1),((j,k,i),1),((k,i,j),1),((j,i,k),-1),((i,k,j),-1),((k,j,i),-1)]: f[a,b,c]=p*sg
    return f
def mulf(f,x,y):
    z=np.zeros(8); z[0]=x[0]*y[0]-x[1:]@y[1:]; z[1:]+=x[0]*y[1:]+y[0]*x[1:]
    z[1:]+=np.einsum('ijk,i,j->k',f[1:,1:,1:],x[1:],y[1:]); return z
rng=np.random.default_rng(0); X=[rng.normal(size=8) for _ in range(4)]
valid=[]
for s in itertools.product([1,-1],repeat=7):
    f=table(s); ok=True
    for a in range(0,4,2):
        x,y=X[a],X[a+1]
        if abs(np.linalg.norm(mulf(f,x,y))-np.linalg.norm(x)*np.linalg.norm(y))>1e-9: ok=False
    if ok: valid.append(s)
print('符号の付け方 全128通りのうち、八元数になるもの:',len(valid))
# 単位の符号反転 e_i -> -e_i は、iを通る3本の線の符号を反転
def flip(s,units):
    s=list(s)
    for n,(i,j,k) in enumerate(lines):
        c=sum(u in (i,j,k) for u in units)
        if c%2: s[n]=-s[n]
    return tuple(s)
base=(1,)*7; orbit=set()
for r in range(8):
    for units in itertools.combinations(range(1,8),r): orbit.add(flip(base,units))
print('単位の符号反転だけで移れる表の数:',len(orbit),' 八元数になる表と一致:',orbit==set(valid))
triv=[units for r in range(8) for units in itertools.combinations(range(1,8),r) if flip(base,units)==base]
print('表を変えない反転の組:',len(triv),triv)
print('e7を通る線:',[l for l in lines if 7 in l]); print('e7を通らない線:',[l for l in lines if 7 not in l])
