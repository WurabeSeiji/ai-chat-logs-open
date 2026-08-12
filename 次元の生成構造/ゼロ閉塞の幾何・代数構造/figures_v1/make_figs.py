import numpy as np, itertools, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.spatial import ConvexHull

rng=np.random.default_rng(2026)
d2=lambda p,q: float(((p-q)**2).sum())
COL={1:"#1f77b4",2:"#d62728",3:"#2ca02c",4:"#9467bd"}

def face_dim_classes(V,tol=1e-7):
    d=V.shape[1]; H=ConvexHull(V); A,b=H.equations[:,:d],H.equations[:,d]
    on=lambda i,f: abs(A[f]@V[i]+b[f])<tol
    out={}
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            F=[f for f in range(len(A)) if on(i,f) and on(j,f)]
            k=d if not F else d-np.linalg.matrix_rank(A[F],tol=1e-8)
            out.setdefault(k,[]).append((i,j))
    return out

def ellipsoid_mesh(Q,n=60):
    """x^T Q x = 1 の楕円体メッシュ"""
    w,U=np.linalg.eigh(Q); ax=1/np.sqrt(w)
    u=np.linspace(0,2*np.pi,n); v=np.linspace(0,np.pi,n//2)
    x=np.outer(np.cos(u),np.sin(v)); y=np.outer(np.sin(u),np.sin(v)); z=np.outer(np.ones_like(u),np.cos(v))
    P=np.stack([x*ax[0],y*ax[1],z*ax[2]],axis=-1)@U.T
    return P[...,0],P[...,1],P[...,2]

def draw(ax,V3,cls,Q=None,title="",elev=22,azim=35):
    if Q is not None:
        X,Y,Z=ellipsoid_mesh(Q)
        ax.plot_surface(X,Y,Z,color="#888888",alpha=0.13,linewidth=0,antialiased=True,zorder=0)
    for k in sorted(cls):
        segs=[[V3[i],V3[j]] for i,j in cls[k]]
        sign="+" if k%2 else "−"
        ax.add_collection3d(Line3DCollection(segs,colors=COL.get(k,"#777"),
              linewidths=2.4 if k==1 else 1.3, alpha=0.95 if k==1 else 0.55,
              label=f"k={k} ({sign})  {len(cls[k])} seg"))
    ax.scatter(V3[:,0],V3[:,1],V3[:,2],c="k",s=34,depthshade=False,zorder=10)
    r=np.abs(V3).max()*1.12
    ax.set_xlim(-r,r); ax.set_ylim(-r,r); ax.set_zlim(-r,r)
    ax.set_box_aspect((1,1,1)); ax.view_init(elev=elev,azim=azim)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_title(title,fontsize=11,pad=2)
    for pane in (ax.xaxis,ax.yaxis,ax.zaxis): pane.pane.fill=False; pane.pane.set_edgecolor("#dddddd")

def alt_and_table(V,cls):
    rows=[]; alt=0.0
    for k in sorted(cls):
        s=sum(d2(V[i],V[j]) for i,j in cls[k]); alt+=(-1)**(k+1)*s
        rows.append((k,len(cls[k]),s))
    return rows,alt

# ---------- Fig 1 : d=3 skew parallelotope on its ellipsoid ----------
A=np.array([[1.00,0.00,0.00],[0.42,0.95,0.00],[0.28,-0.33,0.88]]).T*1.15
S=np.array(list(itertools.product([1,-1],repeat=3)),dtype=float)
V=S@A.T/2
Ai=np.linalg.inv(A/2); Q=Ai.T@Ai/3
cls=face_dim_classes(V); rows,alt=alt_and_table(V,cls)
fig=plt.figure(figsize=(13.2,5.6))
ax=fig.add_subplot(121,projection="3d")
draw(ax,V,cls,Q,"(a)  d=3 skew parallelotope,  N=8,  M=28\nall 8 vertices lie exactly on the ellipsoid")
ax.legend(loc="upper left",fontsize=8,framealpha=0.9)
ax2=fig.add_subplot(122); ax2.axis("off")
txt =("Classification: k = dimension of the smallest face\ncontaining both endpoints.  Sign = (-1)^(k+1)\n\n")
txt+=f"{'class':>7}{'segments':>11}{'sum d^2':>14}{'sign':>7}\n"+"-"*40+"\n"
for k,n,s in rows: txt+=f"{'k='+str(k):>7}{n:>11}{s:>14.6f}{'+' if k%2 else '−':>7}\n"
txt+="-"*40+f"\n{'alternating sum':>18} = {alt:+.3e}\n\n"
resid=max(abs(v@Q@v-1) for v in V)
txt+=f"max |x^T Q x - 1| over the 8 vertices = {resid:.2e}\n"
txt+="\nClosed form (any generators, orthogonal or not):\n"
txt+="   sum over class k  =  2^(d-1) C(d-1,k-1) * sum|u_i|^2\n"
txt+="   alternating sum   =  2^(d-1) * sum|u_i|^2 * (1-1)^(d-1)  =  0\n"
txt+="\nInertia ellipsoid and circumscribed ellipsoid coincide:\n"
T=V.T@V; txt+=f"   T = 2^d A A^T      max deviation {np.abs(T-2**3*(A/2)@(A/2).T).max():.1e}\n"
ax2.text(0.0,0.98,txt,va="top",family="monospace",fontsize=9.4)
fig.suptitle("Zero closure selects parallelotopes — the signed sum of squared lengths vanishes",fontsize=12.5,y=0.99)
fig.tight_layout(); fig.savefig("fig1_parallelotope_d3.png",dpi=190,bbox_inches="tight"); plt.close(fig)
print("fig1 done  alt=",alt,"resid=",resid)

# ---------- Fig 2 : d=4 parallelotope (N=16, M=120) projected into 3D ----------
A4=np.array([[1,0,0,0],[0.35,0.94,0,0],[0.20,-0.28,0.94,0],[-0.15,0.22,0.31,0.91]],dtype=float).T
S4=np.array(list(itertools.product([1,-1],repeat=4)),dtype=float)
V4=S4@A4.T/2
cls4={}
for i in range(16):
    for j in range(i+1,16):
        k=int(sum(x!=y for x,y in zip(S4[i],S4[j]))); cls4.setdefault(k,[]).append((i,j))
rows4,alt4=alt_and_table(V4,cls4)
# 4D -> 3D : 慣性テンソルの上位3主軸へ直交射影
T4=V4.T@V4; w,U=np.linalg.eigh(T4); P=U[:,::-1][:,:3]
V4p=V4@P
fig=plt.figure(figsize=(13.2,5.6))
ax=fig.add_subplot(121,projection="3d")
for k in sorted(cls4):
    segs=[[V4p[i],V4p[j]] for i,j in cls4[k]]
    ax.add_collection3d(Line3DCollection(segs,colors=COL[k],
        linewidths=2.2 if k==1 else 1.0, alpha=0.9 if k==1 else 0.4,
        label=f"k={k} ({'+' if k%2 else '−'})  {len(cls4[k])} seg"))
ax.scatter(V4p[:,0],V4p[:,1],V4p[:,2],c="k",s=30,depthshade=False)
r=np.abs(V4p).max()*1.12
ax.set_xlim(-r,r); ax.set_ylim(-r,r); ax.set_zlim(-r,r); ax.set_box_aspect((1,1,1))
ax.view_init(elev=20,azim=30); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
ax.set_title("(b)  d=4 parallelotope,  N=16,  M=120\northogonal projection onto the 3 principal axes",fontsize=11,pad=2)
ax.legend(loc="upper left",fontsize=8,framealpha=0.9)
for pane in (ax.xaxis,ax.yaxis,ax.zaxis): pane.pane.fill=False; pane.pane.set_edgecolor("#dddddd")
ax2=fig.add_subplot(122); ax2.axis("off")
Ai4=np.linalg.inv(A4/2); Q4=Ai4.T@Ai4/4
resid4=max(abs(v@Q4@v-1) for v in V4)
txt =("The 4-dimensional case.  M = 16*15/2 = 120 segments\nsplit into exactly 4 classes.\n\n")
txt+=f"{'class':>7}{'segments':>11}{'sum d^2':>14}{'sign':>7}\n"+"-"*40+"\n"
for k,n,s in rows4: txt+=f"{'k='+str(k):>7}{n:>11}{s:>14.6f}{'+' if k%2 else '−':>7}\n"
txt+="-"*40+f"\n{'total segments':>18} = {sum(n for _,n,_ in rows4)}\n"
txt+=f"{'alternating sum':>18} = {alt4:+.3e}\n\n"
txt+=f"max |x^T Q x - 1| over the 16 vertices = {resid4:.2e}\n"
txt+="   (the 4D ellipsoid is exact; the picture is its 3D shadow)\n\n"
txt+="Counts for the regular tesseract: 32 / 48 / 32 / 8\n"
txt+="   32 - 96 + 96 - 32 = 0     (weighted by d^2 = k)\n\n"
txt+="Number of pairs 2^(d-1) = 8  <=  d(d+1)/2 = 10,\nbut the identity holds for every d, including d >= 5\nwhere 2^(d-1) exceeds d(d+1)/2."
ax2.text(0.0,0.98,txt,va="top",family="monospace",fontsize=9.4)
fig.suptitle("A 4-dimensional zero closure, projected into 3 dimensions",fontsize=12.5,y=0.99)
fig.tight_layout(); fig.savefig("fig2_parallelotope_d4_projection.png",dpi=190,bbox_inches="tight"); plt.close(fig)
print("fig2 done  alt=",alt4,"resid=",resid4)
