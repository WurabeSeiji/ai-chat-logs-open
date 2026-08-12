import numpy as np, itertools
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.spatial import ConvexHull
d2=lambda p,q: float(((p-q)**2).sum())
COL={1:"#1f77b4",2:"#d62728",3:"#2ca02c"}

A=np.array([[1.00,0,0],[0.42,0.95,0],[0.28,-0.33,0.88]]).T*1.15/2   # 列 = u_i/2
S=np.array(list(itertools.product([1,-1],repeat=3)),dtype=float)
V=S@A.T
d=3; N=len(V)
T=V.T@V                                   # 慣性テンソル
lam,U=np.linalg.eigh(T)                   # 固有値＝半軸²（スケール前）
# 全頂点を通る楕円体: x^T T^-1 x = d / 2^d
c=d/2**d
semi=np.sqrt(lam*c)                       # 半軸 = √(c·λ)
order=np.argsort(-semi); semi=semi[order]; U=U[:,order]
resid=max(abs(v@np.linalg.inv(T)@v - c) for v in V)
ABC=semi                                  # ← A, B, C（序数的な抽象名）
r2=float((ABC**2).sum())
Rbar2=float((V**2).sum())/N               # 平均二乗半径
sumd2=sum(d2(V[i],V[j]) for i in range(N) for j in range(i+1,N))

H=ConvexHull(V); Ae,be=H.equations[:,:3],H.equations[:,3]
on=lambda i,f: abs(Ae[f]@V[i]+be[f])<1e-7
cls={}
for i in range(N):
    for j in range(i+1,N):
        F=[f for f in range(len(Ae)) if on(i,f) and on(j,f)]
        k=3 if not F else 3-np.linalg.matrix_rank(Ae[F],tol=1e-8); cls.setdefault(k,[]).append((i,j))
alt=sum((-1)**(k+1)*sum(d2(V[i],V[j]) for i,j in cls[k]) for k in cls)

fig=plt.figure(figsize=(13.6,6.0))
ax=fig.add_subplot(121,projection="3d")
u=np.linspace(0,2*np.pi,80); v=np.linspace(0,np.pi,40)
X=np.outer(np.cos(u),np.sin(v))*semi[0]; Y=np.outer(np.sin(u),np.sin(v))*semi[1]; Z=np.outer(np.ones_like(u),np.cos(v))*semi[2]
P=np.stack([X,Y,Z],axis=-1)@U.T
ax.plot_wireframe(P[...,0],P[...,1],P[...,2],rstride=6,cstride=6,color="#999999",linewidth=0.5,alpha=0.55)
for k in sorted(cls):
    ax.add_collection3d(Line3DCollection([[V[i],V[j]] for i,j in cls[k]],colors=COL[k],
        linewidths=2.4 if k==1 else 1.1, alpha=0.95 if k==1 else 0.45,
        label=f"k={k} ({'+' if k%2 else '−'})  {len(cls[k])} seg"))
ax.scatter(V[:,0],V[:,1],V[:,2],c="k",s=40,depthshade=False)
lbl=["A","B","C"]; colax=["#c00000","#0060c0","#00a000"]
for a3 in range(3):
    e=U[:,a3]*semi[a3]
    ax.plot(*zip(-e,e),color=colax[a3],linewidth=3.0,zorder=20)
    ax.text(*(e*1.14),lbl[a3],color=colax[a3],fontsize=17,fontweight="bold",ha="center")
rr=semi.max()*1.25
ax.set_xlim(-rr,rr); ax.set_ylim(-rr,rr); ax.set_zlim(-rr,rr); ax.set_box_aspect((1,1,1))
ax.view_init(elev=20,azim=32); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
ax.set_title("the three semi-axes of the ellipsoid are  A, B, C",fontsize=12,pad=0)
ax.legend(loc="upper left",fontsize=8.4,framealpha=0.9)
for p in (ax.xaxis,ax.yaxis,ax.zaxis): p.pane.fill=False; p.pane.set_edgecolor("#e8e8e8")

ax2=fig.add_subplot(122); ax2.axis("off")
t_,R_,Q_=ABC
txt =f"""Skew parallelotope,  d = 3,  N = 2^d = 8,  M = 28

The ellipsoid  x^T T^(-1) x = d / 2^d  passes through
all 8 vertices.   max residual = {resid:.2e}

   semi-axis  A = {t_:.6f}
   semi-axis  B = {R_:.6f}
   semi-axis  C = {Q_:.6f}

   A^2 + B^2 + C^2 = {r2:.6f}  =  r^2      <-- fixed by the closure
   individual A, B, C are NOT fixed

Lagrange identity (Claim 5):
   sum_(i<j) d_ij^2 = N * tr(T) = {sumd2:.6f}
   N^2 * (A^2+B^2+C^2)/d        = {N*N*r2/d:.6f}

Signed closure (Claim 4):
   k=1 (+) {len(cls[1]):2d} seg    k=2 (-) {len(cls[2]):2d} seg    k=3 (+) {len(cls[3]):2d} seg
   alternating sum = {alt:+.2e}

What the closure fixes and what it leaves free (Claim 6):
   l=0   1 dof   size  A^2+B^2+C^2        FIXED
   l=1   3 dof   centre offset            zero by symmetry
   l=2   5 dof   2 ratios + 3 orientations  FREE
   l>=3          not carried by a quadratic form

The names A, B, C are ordinal only (by eigenvalue size).
Whether they correspond to physical spacetime axes
is an open question -- there is no basis for the identification.
"""
ax2.text(0.0,0.99,txt,va="top",family="monospace",fontsize=10.2)
fig.suptitle("The ellipsoid of the zero closure:  its three semi-axes are A, B, C,  and the closure fixes only  A²+B²+C² = r²",
             fontsize=12.5,y=0.995)
fig.tight_layout(); fig.savefig("fig4_semiaxes_ABC.png",dpi=190,bbox_inches="tight"); plt.close(fig)
print("semi-axes A,B,C =",ABC,"  r^2 =",r2,"  resid",resid,"  alt",alt)
print("sum d2 =",sumd2,"  N^2 r^2/d =",N*N*r2/d)
