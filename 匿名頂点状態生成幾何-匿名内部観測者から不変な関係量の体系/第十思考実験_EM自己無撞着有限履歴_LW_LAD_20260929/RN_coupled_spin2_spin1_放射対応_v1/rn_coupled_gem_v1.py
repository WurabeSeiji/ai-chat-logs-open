import argparse, json, math
import numpy as np
from scipy.interpolate import interp1d


def rn_rstar_grid(x, M, Q):
    q = abs(Q)
    if q >= M:
        raise ValueError('black-hole scattering domain requires |Q| < M for this v1 implementation')
    rp = M + math.sqrt(M*M-q*q)
    rm = M - math.sqrt(M*M-q*q)
    # logarithmic radial grid resolves the outer horizon and extends far away
    eps = 1e-10 * max(1.0, rp)
    rmin = rp + eps
    rmax = max(600.0*M, float(np.max(x))+50.0*M)
    y = np.linspace(math.log(rmin-rp), math.log(rmax-rp), 250000)
    r = rp + np.exp(y)
    if q < 1e-14:
        rs = r + 2*M*np.log(r/(2*M)-1.0)
    else:
        d = rp-rm
        rs = r + (rp*rp/d)*np.log(r/rp-1.0) - (rm*rm/d)*np.log(r/rm-1.0)
    # shift r*=0 at r=3M (or just outside horizon if needed)
    rref = max(3.0*M, rp*1.05)
    rs_ref = float(np.interp(rref, r, rs))
    rs = rs-rs_ref
    # ensure requested x is covered
    if x[0] < rs[0] or x[-1] > rs[-1]:
        raise RuntimeError(f'r* interpolation range insufficient: have [{rs[0]}, {rs[-1]}], need [{x[0]}, {x[-1]}]')
    return np.interp(x, rs, r)


def potential_matrix(r, M, Q, ell):
    f = 1.0 - 2.0*M/r + Q*Q/r**2
    L = ell*(ell+1)
    mu2 = (ell-1)*(ell+2)
    Vgg = f*(L/r**2 - 6*M/r**3 + 4*Q*Q/r**4)
    Vee = f*(L/r**2 + 4*Q*Q/r**4)
    Vge = f*(2*Q*math.sqrt(mu2)/r**3)
    return f, Vgg, Vee, Vge


def exact_decoupled_potentials(r, M, Q, ell):
    f = 1.0 - 2.0*M/r + Q*Q/r**2
    L=ell*(ell+1)
    mu2=(ell-1)*(ell+2)
    s=np.sqrt(9*M*M+4*Q*Q*mu2)
    q1=3*M+s
    q2=3*M-s
    # q1 branch -> gravitational RW as Q->0, q2 branch -> EM as Q->0
    Vgrav=f/r**2*(L-q1/r+4*Q*Q/r**2)
    Vem=f/r**2*(L-q2/r+4*Q*Q/r**2)
    return Vgrav, Vem, q1, q2


def canonical_energy(Gt,Et,Gx,Ex,G,E,Vgg,Vee,Vge,dx):
    dens=0.5*(Gt*Gt+Et*Et+Gx*Gx+Ex*Ex+Vgg*G*G+Vee*E*E+2*Vge*G*E)
    return float(np.sum(dens)*dx)


def run(M=1.0,Q=0.0,ell=2,xmin=-40.0,xmax=220.0,nx=5201,cfl=0.42,tmax=420.0,
        x0=5.0,sigma=7.0,amp=1e-3,initial_channel='grav',initial='displacement',sample_stride=100):
    x=np.linspace(xmin,xmax,nx); dx=x[1]-x[0]; dt=cfl*dx; nsteps=int(tmax/dt)
    r=rn_rstar_grid(x,M,Q)
    f,Vgg,Vee,Vge=potential_matrix(r,M,Q,ell)
    V1,V2,q1,q2=exact_decoupled_potentials(r,M,Q,ell)
    # eigenvalue check of coupled matrix
    tr=Vgg+Vee; disc=np.sqrt((Vgg-Vee)**2+4*Vge**2)
    evlo=0.5*(tr-disc); evhi=0.5*(tr+disc)
    exactlo=np.minimum(V1,V2); exacthi=np.maximum(V1,V2)
    eigen_err=float(max(np.max(np.abs(evlo-exactlo)),np.max(np.abs(evhi-exacthi))))

    g=amp*np.exp(-0.5*((x-x0)/sigma)**2)
    gx=-(x-x0)/(sigma*sigma)*g
    G0=np.zeros_like(x); E0f=np.zeros_like(x); Gpi=np.zeros_like(x); Epi=np.zeros_like(x)
    if initial_channel=='grav': G0=g.copy(); base_x=gx
    elif initial_channel=='em': E0f=g.copy(); base_x=gx
    else: raise ValueError(initial_channel)
    if initial=='outgoing':
        if initial_channel=='grav': Gpi=-base_x
        else: Epi=-base_x
    elif initial=='ingoing':
        if initial_channel=='grav': Gpi=base_x
        else: Epi=base_x
    elif initial!='displacement': raise ValueError(initial)

    def lap(a):
        out=np.zeros_like(a); out[1:-1]=(a[2:]-2*a[1:-1]+a[:-2])/dx**2; return out
    Gacc=lap(G0)-Vgg*G0-Vge*E0f
    Eacc=lap(E0f)-Vge*G0-Vee*E0f
    Gprev=G0-dt*Gpi+0.5*dt*dt*Gacc; Eprev=E0f-dt*Epi+0.5*dt*dt*Eacc
    G=G0.copy(); E=E0f.copy(); Gprev[0]=G0[0]-dt*Gpi[0]; Gprev[-1]=G0[-1]-dt*Gpi[-1]; Eprev[0]=E0f[0]-dt*Epi[0]; Eprev[-1]=E0f[-1]-dt*Epi[-1]

    Gx=np.gradient(G0,dx); Ex=np.gradient(E0f,dx)
    Etot0=canonical_energy(Gpi,Epi,Gx,Ex,G0,E0f,Vgg,Vee,Vge,dx)
    flux={'H_grav':0.0,'H_em':0.0,'I_grav':0.0,'I_em':0.0}; samples=[]
    c=dt/dx
    for n in range(nsteps):
        Gn=np.empty_like(G); En=np.empty_like(E)
        Gn[1:-1]=(2*G[1:-1]-Gprev[1:-1]+c*c*(G[2:]-2*G[1:-1]+G[:-2])-dt*dt*(Vgg[1:-1]*G[1:-1]+Vge[1:-1]*E[1:-1]))
        En[1:-1]=(2*E[1:-1]-Eprev[1:-1]+c*c*(E[2:]-2*E[1:-1]+E[:-2])-dt*dt*(Vge[1:-1]*G[1:-1]+Vee[1:-1]*E[1:-1]))
        Gn[0]=G[0]+c*(G[1]-G[0]); Gn[-1]=G[-1]-c*(G[-1]-G[-2])
        En[0]=E[0]+c*(E[1]-E[0]); En[-1]=E[-1]-c*(E[-1]-E[-2])
        Gt=(Gn-Gprev)/(2*dt); Et=(En-Eprev)/(2*dt); Gx=np.gradient(G,dx); Ex=np.gradient(E,dx)
        Sg=-Gt*Gx; Se=-Et*Ex
        flux['H_grav']+=max(0.0,-Sg[1])*dt; flux['H_em']+=max(0.0,-Se[1])*dt
        flux['I_grav']+=max(0.0,Sg[-2])*dt; flux['I_em']+=max(0.0,Se[-2])*dt
        if n%sample_stride==0 or n==nsteps-1:
            Ed=canonical_energy(Gt,Et,Gx,Ex,G,E,Vgg,Vee,Vge,dx)
            samples.append({'t':float(n*dt),'E_domain':Ed,
                            **{k:float(v) for k,v in flux.items()},
                            'closure':float((Ed+sum(flux.values()))/Etot0 if Etot0 else 0.0)})
        Gprev,G=G,Gn; Eprev,E=E,En
    Gt=(G-Gprev)/dt; Et=(E-Eprev)/dt; Gx=np.gradient(G,dx); Ex=np.gradient(E,dx)
    Efinal=canonical_energy(Gt,Et,Gx,Ex,G,E,Vgg,Vee,Vge,dx)
    out={
      'equation':'RN odd-parity coupled gravitational-electromagnetic Moncrief system',
      'M':M,'Q':Q,'Q_over_M':Q/M,'ell':ell,'initial_channel':initial_channel,'initial_type':initial,
      'grid':{'xmin':xmin,'xmax':xmax,'nx':nx,'dx':dx,'dt':dt,'tmax':tmax,'steps':nsteps},
      'initial':{'x0':x0,'sigma':sigma,'amp':amp},
      'q1':float(q1),'q2':float(q2),'decoupling_eigenvalue_max_error':eigen_err,
      'E0':Etot0,'Efinal_domain':Efinal, **{k:float(v) for k,v in flux.items()},
      'fractions':{k:float(v/Etot0) for k,v in flux.items()},
      'fraction_remaining':float(Efinal/Etot0),'closure':float((Efinal+sum(flux.values()))/Etot0),
      'samples':samples
    }
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--Q',type=float,required=True); ap.add_argument('--channel',choices=['grav','em'],default='grav'); ap.add_argument('--out',required=True); ap.add_argument('--tmax',type=float,default=420.0)
    args=ap.parse_args(); res=run(Q=args.Q,initial_channel=args.channel,tmax=args.tmax)
    with open(args.out,'w') as f: json.dump(res,f,indent=2)
    s={k:res[k] for k in ['Q_over_M','initial_channel','decoupling_eigenvalue_max_error','E0','H_grav','H_em','I_grav','I_em','fraction_remaining','closure']}
    s['fractions']=res['fractions']; print(json.dumps(s,indent=2))
if __name__=='__main__': main()
