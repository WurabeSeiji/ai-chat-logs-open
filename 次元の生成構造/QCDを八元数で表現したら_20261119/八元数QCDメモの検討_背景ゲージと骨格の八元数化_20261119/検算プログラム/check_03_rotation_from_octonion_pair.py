import numpy as np
f=np.zeros((8,8,8))
for (i,j,k) in [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]:
    for (a,b,c),s in [((i,j,k),1),((j,k,i),1),((k,i,j),1),((j,i,k),-1),((i,k,j),-1),((k,j,i),-1)]: f[a,b,c]=s
def mul(x,y):
    z=np.zeros(8); z[0]=x[0]*y[0]-x[1:]@y[1:]; z[1:]+=x[0]*y[1:]+y[0]*x[1:]
    z[1:]+=np.einsum('ijk,i,j->k',f[1:,1:,1:],x[1:],y[1:]); return z
E=np.eye(8)
def D(a,b):
    def act(x):
        c=mul(a,b)-mul(b,a)
        return mul(c,x)-mul(x,c)-3*(mul(mul(a,b),x)-mul(a,mul(b,x)))
    return np.array([act(E[i]) for i in range(8)]).T
rng=np.random.default_rng(3)
a,b,x,y=[rng.normal(size=8) for _ in range(4)]
M=D(a,b)
print('derivation err',np.abs(M@mul(x,y)-mul(M@x,y)-mul(x,M@y)).max())
print('antisym',np.abs(M+M.T).max(),' rule f_ijk M_jk',np.abs(np.einsum('ijk,jk->i',f[1:,1:,1:],M[1:,1:])).max())
print('M e7 (generic pair):',np.abs(M@E[7]).max())
# span of D(e_i,e_j)
S=np.array([D(E[i],E[j]).ravel() for i in range(1,8) for j in range(i+1,8)])
print('rank all pairs',np.linalg.matrix_rank(S))
# pairs that kill e7: combos; find dimension of subspace of span killing e7
U,s,Vt=np.linalg.svd(S); basis=Vt[:14]
K=np.array([ (v.reshape(8,8)@E[7]) for v in basis])
print('dim killing e7',14-np.linalg.matrix_rank(K))
# commutator closes
a2,b2=rng.normal(size=8),rng.normal(size=8); N=D(a2,b2); C=M@N-N@M
print('commutator nonzero',np.abs(C).max(),' still obeys rule',np.abs(np.einsum('ijk,jk->i',f[1:,1:,1:],C[1:,1:])).max())
