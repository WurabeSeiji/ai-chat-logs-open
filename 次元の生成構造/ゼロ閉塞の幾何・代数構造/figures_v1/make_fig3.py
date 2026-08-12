import numpy as np, itertools
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.spatial import ConvexHull
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

def panel(ax,V,title,show_sphere=True):
    cls=face_dim_classes(V); rows=[];alt=0
    if show_sphere:
        R=np.linalg.norm(V[0]); u=np.linspace(0,2*np.pi,60); v=np.linspace(0,np.pi,30)
        x=R*np.outer(np.cos(u),np.sin(v)); y=R*np.outer(np.sin(u),np.sin(v)); z=R*np.outer(np.ones_like(u),np.cos(v))
        ax.plot_surface(x,y,z,color="#888888",alpha=0.12,linewidth=0)
    for k in sorted(cls):
        s=sum(d2(V[i],V[j]) for i,j in cls[k]); alt+=(-1)**(k+1)*s; rows.append((k,len(cls[k]),s))
        ax.add_collection3d(Line3DCollection([[V[i],V[j]] for i,j in cls[k]],colors=COL.get(k,"#777"),
            linewidths=2.2 if k==1 else 0.9, alpha=0.9 if k==1 else 0.32,
            label=f"k={k} ({'+' if k%2 else '−'})  {len(cls[k])}"))
    ax.scatter(V[:,0],V[:,1],V[:,2],c="k",s=28,depthshade=False)
    r=np.abs(V).max()*1.1
    ax.set_xlim(-r,r); ax.set_ylim(-r,r); ax.set_zlim(-r,r); ax.set_box_aspect((1,1,1))
    ax.view_init(elev=18,azim=28); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_title(title,fontsize=10.5,pad=1)
    ax.legend(loc="upper left",fontsize=7.4,framealpha=0.85)
    for p in (ax.xaxis,ax.yaxis,ax.zaxis): p.pane.fill=False; p.pane.set_edgecolor("#dddddd")
    return rows,alt

phi=(1+5**0.5)/2
ico=np.unique(np.round(np.array([p for c in [(0,1,phi),(1,phi,0),(phi,0,1)] for p in
      {(s1*c[0],s2*c[1],s3*c[2]) for s1 in(1,-1) for s2 in(1,-1) for s3 in(1,-1)}],dtype=float),9),axis=0)
octa=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],dtype=float)
A=np.array([[1.00,0,0],[0.42,0.95,0],[0.28,-0.33,0.88]]).T*1.15
par=np.array(list(itertools.product([1,-1],repeat=3)),dtype=float)@A.T/2

fig=plt.figure(figsize=(14.5,5.4))
a1=fig.add_subplot(131,projection="3d"); r1,alt1=panel(a1,ico,"(a) icosahedron  N=12, M=66\ncentrally symmetric, convex, on a sphere")
a2=fig.add_subplot(132,projection="3d"); r2,alt2=panel(a2,octa,"(b) octahedron  N=6, M=15\ncentrally symmetric, convex, on a sphere")
a3=fig.add_subplot(133,projection="3d"); r3,alt3=panel(a3,par,"(c) skew parallelotope  N=8, M=28\ncentrally symmetric, convex, on an ellipsoid",show_sphere=False)
Ai=np.linalg.inv(A/2); Q=Ai.T@Ai/3
w,U=np.linalg.eigh(Q); ax=1/np.sqrt(w)
u=np.linspace(0,2*np.pi,60); v=np.linspace(0,np.pi,30)
X=np.outer(np.cos(u),np.sin(v))*ax[0]; Y=np.outer(np.sin(u),np.sin(v))*ax[1]; Z=np.outer(np.ones_like(u),np.cos(v))*ax[2]
P=np.stack([X,Y,Z],axis=-1)@U.T
a3.plot_surface(P[...,0],P[...,1],P[...,2],color="#888888",alpha=0.12,linewidth=0)

def cap(rows,alt):
    s=" / ".join(f"k={k}: {n} seg, Σd²={v:.3f}" for k,n,v in rows)
    return s+f"\nalternating sum = {alt:+.4f}"
fig.text(0.175,0.045,cap(r1,alt1),ha="center",family="monospace",fontsize=8.4,
         color="#b00000" if abs(alt1)>1e-9 else "#006000")
fig.text(0.505,0.045,cap(r2,alt2),ha="center",family="monospace",fontsize=8.4,
         color="#b00000" if abs(alt2)>1e-9 else "#006000")
fig.text(0.835,0.045,cap(r3,alt3),ha="center",family="monospace",fontsize=8.4,
         color="#b00000" if abs(alt3)>1e-9 else "#006000")
fig.suptitle("Lying on an ellipsoid is not enough:  (a) and (b) are centrally symmetric, convex and inscribed in a sphere,\n"
             "yet only the parallelotope (c) makes the signed sum vanish",fontsize=12,y=1.02)
fig.tight_layout(rect=[0,0.10,1,1])
fig.savefig("fig3_counterexamples.png",dpi=190,bbox_inches="tight"); plt.close(fig)
print("fig3", alt1, alt2, alt3)
