import numpy as np
f=np.zeros((8,8,8))
for (i,j,k) in [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]:
    for (a,b,c),s in [((i,j,k),1),((j,k,i),1),((k,i,j),1),((j,i,k),-1),((i,k,j),-1),((k,j,i),-1)]: f[a,b,c]=s
def mul(x,y):
    z=np.zeros(8); z[0]=x[0]*y[0]-x[1:]@y[1:]; z[1:]+=x[0]*y[1:]+y[0]*x[1:]
    z[1:]+=np.einsum('ijk,i,j->k',f[1:,1:,1:],x[1:],y[1:]); return z
E=np.eye(8); conj=lambda z: np.concatenate([[z[0]],-z[1:]])
L=lambda a: np.array([mul(a,E[i]) for i in range(8)]).T
def unit(v): return v/np.linalg.norm(v)
rng=np.random.default_rng(5)
def test(dim,label):
    q=[]
    for _ in range(3):
        v=np.zeros(8); v[:dim]=rng.normal(size=dim); q.append(unit(v))
    a,b,c=q
    # (A) 交差項を先に取り、それを掛ける: U_ij = L(q_i conj(q_j))
    H=L(mul(a,conj(b)))@L(mul(b,conj(c)))@L(mul(c,conj(a)))
    # (B) 値を順に掛ける: U_ij = L(q_i) L(conj q_j)
    H2=L(a)@L(conj(b))@L(b)@L(conj(c))@L(c)@L(conj(a))
    print(label,'| (A) 一周のずれ |H-1| =',round(np.abs(H-np.eye(8)).max(),6),' Tr(1-H) =',round(np.trace(np.eye(8)-H),6),'| (B) |H-1| =',round(np.abs(H2-np.eye(8)).max(),6))
    return H
test(2,'複素数 ')
test(4,'四元数 ')
H=test(8,'八元数 ')
print('直交行列か',np.abs(H@H.T-np.eye(8)).max(),' det',round(np.linalg.det(H),6))
print('1を動かすか |H e0 - e0| =',round(np.abs(H@E[0]-E[0]).max(),6))
x,y=rng.normal(size=8),rng.normal(size=8)
print('積を保つか（自己同型か） err =',round(np.abs(H@mul(x,y)-mul(H@x,H@y)).max(),6))
print('--- 自分の代数の中だけで見た場合 ---')
for dim,label in [(2,'複素数'),(4,'四元数'),(8,'八元数')]:
    devs=[]
    for t in range(200):
        q=[]
        for _ in range(3):
            v=np.zeros(8); v[:dim]=rng.normal(size=dim); q.append(unit(v))
        a,b,c=q
        H=L(mul(a,conj(b)))@L(mul(b,conj(c)))@L(mul(c,conj(a)))
        Hs=H[:dim,:dim]
        devs.append(np.trace(np.eye(dim)-Hs))
    print(label,' Tr(1-H) 平均 =',round(float(np.mean(devs)),6),' 最大 =',round(float(np.max(devs)),6))
