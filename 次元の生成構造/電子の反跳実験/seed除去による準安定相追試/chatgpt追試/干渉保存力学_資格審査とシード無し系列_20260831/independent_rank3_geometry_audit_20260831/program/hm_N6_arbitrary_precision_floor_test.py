import mpmath as mp
import time

N=6
M=N*(N-1)//2
L=124

# edge ordering identical to common.py
E=[(i,j) for i in range(N) for j in range(i+1,N)]
# 1-factor classes identical to common.py/state_provider.py
n=N-1
col={}
for rr in range(n):
    col[tuple(sorted((rr,N-1)))]=rr
    for k in range(1,N//2):
        col[tuple(sorted(((rr-k)%n,(rr+k)%n)))]=rr
cls=[col[e] for e in E]

# adjacency on edges
A=[[mp.mpf('0')]*M for _ in range(M)]
for a,(i,j) in enumerate(E):
    for b,(k,l) in enumerate(E):
        if a!=b and (i==k or i==l or j==k or j==l):
            A[a][b]=mp.mpf('1')

def dot_real_basis(p,z):
    return sum(p[i]*z[i] for i in range(M))

def norm2(z):
    return mp.fsum([abs(x)**2 for x in z])

def make_hm():
    # hm_N6 observed mean_amp2 = 0.08; use exact 2/25 to avoid importing double noise
    r=mp.sqrt(mp.mpf(2)/25)
    q=N-1
    return [r*mp.e**(1j*mp.pi*mp.mpf(c)/q) for c in cls]

def make_pq(v):
    vr=[mp.re(x) for x in v]
    vi=[mp.im(x) for x in v]
    nr=mp.sqrt(mp.fsum([x*x for x in vr])); p=[x/nr for x in vr]
    proj=mp.fsum([vi[i]*p[i] for i in range(M)])
    q0=[vi[i]-proj*p[i] for i in range(M)]
    nq=mp.sqrt(mp.fsum([x*x for x in q0])); q=[x/nq for x in q0]
    return p,q

def hperp_frac(z,p,q):
    ap=dot_real_basis(p,z); aq=dot_real_basis(q,z)
    zp=[z[i]-p[i]*ap-q[i]*aq for i in range(M)]
    return norm2(zp)/norm2(z)

def Hmat(z):
    H=mp.matrix(M,M)
    for i in range(M):
        zi=mp.conj(z[i])
        for j in range(M):
            H[i,j]=A[i][j]*zi*z[j]
    return H

def step(z,delta):
    H=Hmat(z)
    vals,V=mp.eighe(H)
    # c=V^H z
    zz=mp.matrix(z)
    c=V.H*zz
    for i in range(M):
        c[i] *= mp.e**(-1j*delta*vals[i])
    out=V*c
    return [out[i] for i in range(M)]

def residual_equilibrium(v):
    H=Hmat(v)
    hv=H*mp.matrix(v)
    den=norm2(v)
    mu=mp.re(sum(mp.conj(v[i])*hv[i] for i in range(M))/den)
    rr=[hv[i]-mu*v[i] for i in range(M)]
    return mu, mp.sqrt(norm2(rr)/den)

def run(dps,subdiv,macro_steps=2):
    mp.mp.dps=dps
    v=make_hm(); p,q=make_pq(v)
    mu,res=residual_equilibrium(v)
    z=list(v)
    delta=2*mp.pi/(mp.mpf(L)*subdiv)
    vals=[hperp_frac(z,p,q)]
    for t in range(macro_steps):
        for _ in range(subdiv):
            z=step(z,delta)
        vals.append(hperp_frac(z,p,q))
    # exact phase solution error after macro steps
    phase=mp.e**(-1j*(2*mp.pi/L)*mu*macro_steps)
    zex=[phase*x for x in v]
    err=mp.sqrt(norm2([z[i]-zex[i] for i in range(M)])/norm2(v))
    return mu,res,vals,err

for dps in [30,50,80]:
    for subdiv in [1,2,4,8]:
        t=time.time()
        mu,res,vals,err=run(dps,subdiv,2)
        print('dps',dps,'sub',subdiv,'mu',mp.nstr(mu,10),'res',mp.nstr(res,6),'f0',mp.nstr(vals[0],8),'f1',mp.nstr(vals[1],8),'f2',mp.nstr(vals[2],8),'err2',mp.nstr(err,8),'sec',round(time.time()-t,2))
