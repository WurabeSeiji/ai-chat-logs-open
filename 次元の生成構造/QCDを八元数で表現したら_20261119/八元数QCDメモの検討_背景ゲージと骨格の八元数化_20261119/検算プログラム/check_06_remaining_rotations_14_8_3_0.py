import numpy as np
f=np.zeros((8,8,8))
for (i,j,k) in [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]:
    for (a,b,c),s in [((i,j,k),1),((j,k,i),1),((k,i,j),1),((j,i,k),-1),((i,k,j),-1),((k,j,i),-1)]: f[a,b,c]=s
def mul(x,y):
    z=np.zeros(8); z[0]=x[0]*y[0]-x[1:]@y[1:]; z[1:]+=x[0]*y[1:]+y[0]*x[1:]
    z[1:]+=np.einsum('ijk,i,j->k',f[1:,1:,1:],x[1:],y[1:]); return z
E=np.eye(8); conj=lambda z: np.concatenate([[z[0]],-z[1:]])
def D(a,b):
    def act(x):
        c=mul(a,b)-mul(b,a); return mul(c,x)-mul(x,c)-3*(mul(mul(a,b),x)-mul(a,mul(b,x)))
    return np.array([act(E[i]) for i in range(8)]).T
S=np.array([D(E[i],E[j]).ravel() for i in range(1,8) for j in range(i+1,8)])
basis=np.linalg.svd(S)[2][:14].reshape(14,8,8)   # 積の表を保つ回転 14個
def remaining(vecs):
    if not vecs: return 14
    K=np.array([np.concatenate([B@v for v in vecs]) for B in basis])
    return 14-np.linalg.matrix_rank(K,tol=1e-9)
rng=np.random.default_rng(11)
q=[rng.normal(size=8) for _ in range(5)]
Y=lambda i,j: mul(q[i],conj(q[j]))          # 交差項
print('交差項なし            : 残る回転',remaining([]))
print('交差項1個 Y12         : 残る回転',remaining([Y(0,1)]))
print('交差項2個 Y12,Y23     : 残る回転',remaining([Y(0,1),Y(1,2)]))
print('交差項3個 Y12,Y23,Y34 : 残る回転',remaining([Y(0,1),Y(1,2),Y(2,3)]))
print('参考 Y12,Y23,Y13      : 残る回転',remaining([Y(0,1),Y(1,2),Y(0,2)]))
# 3個目が最初の2個の作る四元数の中にある場合
u=Y(0,1); v=Y(1,2); w=mul(u,v)
print('3個目が2個の積(同じ四元数内): 残る回転',remaining([u,v,w]))
# 3個で読める符号つき不変量 phi(u,v,w)=f_ijk u_i v_j w_k
phi=lambda a,b,c: np.einsum('ijk,i,j,k->',f[1:,1:,1:],a[1:],b[1:],c[1:])
g=np.eye(8)
import scipy.linalg as sl
R=sl.expm(sum(rng.normal()*B for B in basis))   # 積の表を保つ回転で全体を呼び替え
a,b,c=Y(0,1),Y(1,2),Y(2,3)
print('phi(Y12,Y23,Y34) =',round(phi(a,b,c),6),' 呼び替え後 =',round(phi(R@a,R@b,R@c),6))
