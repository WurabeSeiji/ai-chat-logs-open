import json, math, argparse
import numpy as np

# Full nonlinear Einstein-Maxwell-Lambda exact solution (Kastor-Traschen)
# ds^2 = -Omega^-2 dt^2 + a(t)^2 Omega^2 d x^2
# Omega = 1 + sum M_A/(a r_A), a=exp(H t), Lambda=3 H^2
# Define W = a Omega = a + sum M_A/r_A, then gamma_ij=W^2 delta_ij,
# alpha=a/W, beta=0, K_ij=-H gamma_ij.
# Gauge potential A_t = +Omega^-1 (overall sign flips charge).
# E_i = - d_i W / W, E^i = - d_i W / W^3.


def centered_grad(f, h):
    gx = (f[2:,1:-1,1:-1]-f[:-2,1:-1,1:-1])/(2*h)
    gy = (f[1:-1,2:,1:-1]-f[1:-1,:-2,1:-1])/(2*h)
    gz = (f[1:-1,1:-1,2:]-f[1:-1,1:-1,:-2])/(2*h)
    return gx,gy,gz

def centered_lap(f,h):
    c=f[1:-1,1:-1,1:-1]
    return ((f[2:,1:-1,1:-1]+f[:-2,1:-1,1:-1]+f[1:-1,2:,1:-1]+f[1:-1,:-2,1:-1]+f[1:-1,1:-1,2:]+f[1:-1,1:-1,:-2]-6*c)/(h*h))

def run(N=49, L=12.0, H=-0.05, sep=8.0, m1=1.0, m2=1.0, times=(0.0,4.0,8.0), exclusion=1.2):
    x=np.linspace(-L,L,N); h=x[1]-x[0]
    X,Y,Z=np.meshgrid(x,x,x,indexing='ij')
    centers=[(-sep/2,0,0),(sep/2,0,0)]
    rs=[]
    for cx,cy,cz in centers:
        rs.append(np.sqrt((X-cx)**2+(Y-cy)**2+(Z-cz)**2))
    # avoid exact division at punctures; excluded anyway
    eps=1e-15
    Ustatic=m1/np.maximum(rs[0],eps)+m2/np.maximum(rs[1],eps)
    results=[]
    Lam=3*H*H
    for t in times:
        a=math.exp(H*t)
        W=a+Ustatic
        alpha=a/W
        # numerical derivatives of W
        gx,gy,gz=centered_grad(W,h)
        lap=centered_lap(W,h)
        Wi=W[1:-1,1:-1,1:-1]
        grad2=gx*gx+gy*gy+gz*gz
        # 3-Ricci scalar for gamma_ij=W^2 delta_ij
        R3=-4*lap/(Wi**3)+2*grad2/(Wi**4)
        rho_em=grad2/(8*math.pi*Wi**4)
        # K=-3H, KijKij=3H^2
        Ham=R3+6*H*H-16*math.pi*rho_em-2*Lam
        # Maxwell Gauss: D_i E^i = W^-3 div(W^3 E^i), W^3 E^i=-grad W
        Gauss=-lap/(Wi**3)
        # exact momentum residual is identically zero for K_ij=-H gamma_ij, S_i=0
        # Construct mask away from boundaries and punctures
        r1=rs[0][1:-1,1:-1,1:-1]; r2=rs[1][1:-1,1:-1,1:-1]
        mask=(r1>exclusion)&(r2>exclusion)
        # also avoid very large W near excluded edges of punctures
        def stats(v):
            q=np.abs(v[mask])
            return {"max_abs":float(q.max()),"rms":float(np.sqrt(np.mean(q*q))),"median_abs":float(np.median(q))}
        # constraint normalized scale
        source_scale=np.abs(R3)+np.abs(16*math.pi*rho_em)+2*Lam+1e-30
        normHam=np.abs(Ham)/source_scale
        # Dynamic geometric quantities at sample point origin (midpoint)
        i0=N//2
        W0=float(W[i0,i0,i0]); alpha0=float(alpha[i0,i0,i0])
        bg_sep=a*sep
        # nonlinear aa/ab/bb metric-factor components at origin
        ua=float(m1/max(rs[0][i0,i0,i0],eps)); ub=float(m2/max(rs[1][i0,i0,i0],eps))
        # W^2 = a^2 + 2a(ua+ub) + ua^2 + 2uaub + ub^2
        decomp={"a2":a*a,"linear_a":2*a*ua,"linear_b":2*a*ub,"aa":ua*ua,"ab":2*ua*ub,"bb":ub*ub,
                "sum":a*a+2*a*ua+2*a*ub+ua*ua+2*ua*ub+ub*ub,"W2":W0*W0}
        results.append({
            "t":t,"a":a,"Lambda":Lam,"background_physical_separation":bg_sep,
            "midpoint":{"W":W0,"alpha":alpha0,"decomposition":decomp},
            "constraints":{"Hamiltonian":stats(Ham),"Gauss":stats(Gauss),"Momentum":{"max_abs":0.0,"rms":0.0},
                           "Hamiltonian_normalized":{"max":float(normHam[mask].max()),"rms":float(np.sqrt(np.mean(normHam[mask]**2)))}}
        })
    return {"model":"two-center Kastor-Traschen exact Einstein-Maxwell-Lambda","N":N,"L":L,"h":h,"H":H,"Lambda":Lam,
            "masses":[m1,m2],"charges":[m1,m2],"coordinate_separation":sep,"exclusion":exclusion,"times":list(times),"results":results}

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--N',type=int,default=49); ap.add_argument('--out',default='kt_two_center.json')
    args=ap.parse_args()
    d=run(N=args.N)
    with open(args.out,'w') as f: json.dump(d,f,indent=2)
    print(json.dumps(d,indent=2))
