import numpy as np, itertools, sympy as sp
f=np.zeros((8,8,8),dtype=int)
for (i,j,k) in [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]:
    for (a,b,c),s in [((i,j,k),1),((j,k,i),1),((k,i,j),1),((j,i,k),-1),((i,k,j),-1),((k,j,i),-1)]:
        f[a,b,c]=s
def mul(x,y):
    z=np.zeros(8); z[0]=x[0]*y[0]-x[1:]@y[1:]; z[1:]+=x[0]*y[1:]+y[0]*x[1:]
    z[1:]+=np.einsum('ijk,i,j->k',f[1:,1:,1:],x[1:],y[1:]); return z
E=np.eye(8)
# complex coords: psi_a = alpha_a + i beta_a  <->  alpha_a e_a + beta_a (e7 e_a)
B=np.zeros((8,6))   # columns: Re psi1..3, Im psi1..3 -> octonion
for n in range(3):
    B[:,n]=E[n+1]; B[:,3+n]=mul(E[7],E[n+1])
print('e7*e_a:',[ (n+1, int(np.argmax(abs(B[:,3+n]))), int(B[:,3+n][np.argmax(abs(B[:,3+n]))])) for n in range(3)])
# Gell-Mann
l=np.zeros((8,3,3),dtype=complex)
l[0][0,1]=l[0][1,0]=1; l[1][0,1]=-1j; l[1][1,0]=1j; l[2][0,0]=1; l[2][1,1]=-1
l[3][0,2]=l[3][2,0]=1; l[4][0,2]=-1j; l[4][2,0]=1j; l[5][1,2]=l[5][2,1]=1
l[6][1,2]=-1j; l[6][2,1]=1j; l[7]=np.diag([1,1,-2])/np.sqrt(3)
T=[]
for a in range(8):
    X=-0.5j*l[a]   # acts on C^3
    R=np.block([[X.real,-X.imag],[X.imag,X.real]])  # on (Re,Im)
    M=B@R@B.T   # B has orthonormal columns
    T.append(M)
# checks: derivation property, kills e7, trace
rng=np.random.default_rng(1)
Fa=rng.normal(size=8); F=sum(Fa[a]*T[a] for a in range(8))
x,y=rng.normal(size=8),rng.normal(size=8)
print('derivation err',np.abs(F@mul(x,y)-mul(F@x,y)-mul(x,F@y)).max(),' F e7:',np.abs(F@E[7]).max(),' F 1:',np.abs(F@E[0]).max())
print('TrO(FF)',np.trace(F@F),' -F^aF^a',-Fa@Fa, ' -sum Fij^2',-(F**2).sum())
print('g2 condition f_ijk F_jk:',np.abs(np.einsum('ijk,jk->i',f[1:,1:,1:],F[1:,1:])).max())
# dimension of {A in so(7): f_ijk A_jk=0, A_7j=0}
idx=[(i,j) for i in range(1,8) for j in range(i+1,8)]
rows=[]
for i in range(1,8):
    r=[2*f[i,p,q] for (p,q) in idx]; rows.append(r)
for j in range(1,7):
    r=[1 if (p,q)==(j,7) else 0 for (p,q) in idx]; rows.append(r)
print('dim solution space', len(idx)-np.linalg.matrix_rank(np.array(rows)))
# symbolic table F_ij in terms of F^a
Fs=sp.symbols('F1:9'); Ms=sp.zeros(8,8)
for a in range(8): Ms+=Fs[a]*sp.Matrix(np.round(T[a]*2*np.sqrt(3)).astype(int))/ (2*sp.sqrt(3)) if a==7 else Fs[a]*sp.Matrix(np.round(T[a]*2).astype(int))/2
for a in range(8):
    chk=np.abs(T[a]-np.array((Ms.subs({Fs[b]:(1 if b==a else 0) for b in range(8)})).evalf(),dtype=float)).max()
    assert chk<1e-12, (a,chk)
for i in range(1,7):
    for j in range(i+1,7):
        print('F_%d%d ='%(i,j), sp.simplify(Ms[i,j]))
print('row0/col7 zero:', Ms[0,:], Ms[7,:])
# constraints written out
for i in range(1,8):
    terms=[('+' if f[i,p,q]>0 else '-')+'F_%d%d'%(p,q) for (p,q) in idx if f[i,p,q]!=0]
    print('i=%d:'%i,' '.join(terms),'= 0')
