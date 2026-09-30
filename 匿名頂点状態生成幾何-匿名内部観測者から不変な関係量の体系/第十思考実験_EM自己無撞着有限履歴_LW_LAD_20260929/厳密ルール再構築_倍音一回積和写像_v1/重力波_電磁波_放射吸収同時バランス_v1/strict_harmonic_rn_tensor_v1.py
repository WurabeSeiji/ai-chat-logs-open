#!/usr/bin/env python3
"""Strict harmonic-state RN gravito-EM experiment.

Persistent state: ONLY complex harmonic coefficients Z[v,m].
Transition: ONE fixed sparse quadratic product-sum map
    Z' = T[Z,Z]
compiled once from the discrete RN coupled wave recurrence.

No RK stages, no previous-state arrays, no external Q in transition, no readout feedback.
Q^0..Q^4 are state channels and are preserved by the same quadratic map.
The same compiled interaction array is used for every Q/M case; only initial state differs.

Units: M=c=G=1. Odd-parity RN coupled perturbation benchmark.
The two wave channels are anonymous in generation. Their grav/EM interpretation is readout only,
fixed by the Q->0 limit of the standard coupled system.
"""
from __future__ import annotations
import argparse, json, math, hashlib
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy import sparse

TWOPI=2*np.pi

@dataclass
class Layout:
    nx:int
    K:int
    nfield:int=2
    nq:int=5
    @property
    def nvar(self): return self.nfield*self.nx+self.nq
    @property
    def qbase(self): return self.nfield*self.nx
    @property
    def nstate(self): return self.nvar*self.K
    def field_var(self,ch,i): return ch*self.nx+i
    def q_var(self,p): return self.qbase+p
    def idx(self,v,m): return v*self.K+m


def dft_encode(samples):
    return np.fft.fft(samples,axis=1)/samples.shape[1]

def dft_decode(coeff):
    return np.fft.ifft(coeff*coeff.shape[1],axis=1)


def radial_operator_bases(r,dt,ell=2):
    """Return A[p,ch_out,ch_in,i,j] such that newest wave sample is
       y_out = sum_p q^p A_p @ current + (-previous) interior.
       Boundary rows implement first-order ingoing/outgoing characteristic update.
       p=0..4. No Q is used here.
    """
    nx=len(r); dr=r[1]-r[0]; L=ell*(ell+1); mu2=(ell-1)*(ell+2)
    A=np.zeros((5,2,2,nx,nx),float)
    # polynomial pieces f=f0+q^2 f2; f'=fp0+q^2 fp2
    f0=1-2/r; f2=1/r**2
    fp0=2/r**2; fp2=-2/r**3
    f2b=[f0*f0, np.zeros(nx), 2*f0*f2, np.zeros(nx), f2*f2]
    ffp=[f0*fp0, np.zeros(nx), f0*fp2+f2*fp0, np.zeros(nx), f2*fp2]
    # V diagonal polynomial pieces
    gg0=L/r**2-6/r**3; gg2=4/r**4
    ee0=L/r**2; ee2=4/r**4
    Vgg=[f0*gg0,np.zeros(nx),f0*gg2+f2*gg0,np.zeros(nx),f2*gg2]
    Vee=[f0*ee0,np.zeros(nx),f0*ee2+f2*ee0,np.zeros(nx),f2*ee2]
    b1=2*math.sqrt(mu2)/r**3
    Vge=[np.zeros(nx),f0*b1,np.zeros(nx),f2*b1,np.zeros(nx)]
    # interior: direct recurrence y=2cur-prev+dt^2(...)
    for i in range(1,nx-1):
        for p in (0,2,4):
            # f^2 D2 + f f' D1
            c2=dt*dt*f2b[p][i]/(dr*dr)
            c1=dt*dt*ffp[p][i]/(2*dr)
            for ch in (0,1):
                A[p,ch,ch,i,i-1]+=c2-c1
                A[p,ch,ch,i,i]+=-2*c2
                A[p,ch,ch,i,i+1]+=c2+c1
        # direct 2*current belongs p=0
        A[0,0,0,i,i]+=2.0; A[0,1,1,i,i]+=2.0
        for p in range(5):
            A[p,0,0,i,i]+=-dt*dt*Vgg[p][i]
            A[p,1,1,i,i]+=-dt*dt*Vee[p][i]
            A[p,0,1,i,i]+=-dt*dt*Vge[p][i]
            A[p,1,0,i,i]+=-dt*dt*Vge[p][i]
    # boundaries: characteristic propagation, no separate semantic branch in runtime.
    # left (ingoing toward decreasing r*): y = cur + dt*f*(cur[1]-cur[0])/dr
    # right (outgoing): y = cur - dt*f*(cur[-1]-cur[-2])/dr
    for ch in (0,1):
        A[0,ch,ch,0,0]+=1-dt*f0[0]/dr; A[0,ch,ch,0,1]+=dt*f0[0]/dr
        A[2,ch,ch,0,0]+=-dt*f2[0]/dr; A[2,ch,ch,0,1]+=dt*f2[0]/dr
        A[0,ch,ch,-1,-1]+=1-dt*f0[-1]/dr; A[0,ch,ch,-1,-2]+=dt*f0[-1]/dr
        A[2,ch,ch,-1,-1]+=-dt*f2[-1]/dr; A[2,ch,ch,-1,-2]+=dt*f2[-1]/dr
    return A


def add_sample_rule(rules,vo,ho,v1,h1,v2,h2,coef):
    if coef!=0.0:
        rules.append((vo,ho,v1,h1,v2,h2,complex(coef)))


def compile_sparse_quadratic(layout:Layout,r,dt,ell=2):
    """Compile sample-space bilinear rules through exact DFT similarity into
    one fixed sparse product-sum interaction array M on monomials z_i*z_j.
    Returns M, pair_i, pair_j. Runtime is exactly: out = M @ (z[pair_i]*z[pair_j]).
    """
    K=layout.K; nx=layout.nx; last=K-1; prev=K-2; q0=layout.q_var(0)
    rules=[]
    # history shift for EVERY state variable: s'[h]=s[h+1]*q0_current
    for v in range(layout.nvar):
        for h in range(K-1):
            add_sample_rule(rules,v,h,v,h+1,q0,last,1.0)
    # q^p state channels preserved by same interaction
    for p in range(5):
        qp=layout.q_var(p)
        add_sample_rule(rules,qp,last,qp,last,q0,last,1.0)
    # wave newest-sample recurrence. q powers are explicit state channels.
    A=radial_operator_bases(r,dt,ell)
    for co in (0,1):
        for i in range(nx):
            vo=layout.field_var(co,i)
            # -previous interior only
            if 0<i<nx-1:
                add_sample_rule(rules,vo,last,vo,prev,q0,last,-1.0)
            for p in range(5):
                qp=layout.q_var(p)
                block=A[p,co]
                for ci in (0,1):
                    row=block[ci,i]
                    js=np.nonzero(row)[0]
                    for j in js:
                        add_sample_rule(rules,vo,last,layout.field_var(ci,j),last,qp,last,row[j])
    # Compile each sample rule into harmonic coefficients.
    pair_to_col={}; rows=[]; cols=[]; vals=[]
    roots=(np.array([1+0j,0+1j,-1+0j,0-1j]) if K==4 else np.exp(2j*np.pi*np.arange(K)/K))
    for vo,ho,v1,h1,v2,h2,coef in rules:
        for m1 in range(K):
            e1=roots[m1]**h1
            i1=layout.idx(v1,m1)
            for m2 in range(K):
                e12=e1*(roots[m2]**h2)
                i2=layout.idx(v2,m2)
                pair=(i1,i2)
                col=pair_to_col.get(pair)
                if col is None:
                    col=len(pair_to_col); pair_to_col[pair]=col
                for mo in range(K):
                    c=coef*e12*np.exp(-2j*np.pi*mo*ho/K)/K
                    if abs(c)>1e-15:
                        rows.append(layout.idx(vo,mo)); cols.append(col); vals.append(c)
    pairs=[None]*len(pair_to_col)
    for p,c in pair_to_col.items(): pairs[c]=p
    pi=np.array([p[0] for p in pairs],dtype=np.int32)
    pj=np.array([p[1] for p in pairs],dtype=np.int32)
    M=sparse.coo_matrix((np.array(vals,complex),(np.array(rows),np.array(cols))),shape=(layout.nstate,len(pairs))).tocsr()
    M.sum_duplicates()
    # remove Fourier-cancellation roundoff so exact invariant channels do not square numerical noise
    M.data[np.abs(M.data)<1e-13]=0
    M.eliminate_zeros()
    return M,pi,pj,len(rules)


def transition(z,M,pi,pj):
    # THE ONLY state transition operation in the execution loop.
    mon=z[pi]*z[pj]
    return M @ mon


def init_state(layout,r,q,x0=12.0,sigma=4.0,amp=1e-4,channel=0):
    K=layout.K; s=np.zeros((layout.nvar,K),complex)
    g=amp*np.exp(-0.5*((r-x0)/sigma)**2)
    # constant prehistory = displacement with zero initial time derivative
    for i in range(layout.nx): s[layout.field_var(channel,i),:]=g[i]
    qp=[1.0,q,q*q,q*q*q,q*q*q*q]
    for p,val in enumerate(qp): s[layout.q_var(p),:]=val
    return dft_encode(s).reshape(-1)


def decode_state(z,layout):
    return dft_decode(z.reshape(layout.nvar,layout.K))


def readout(z,zprev,layout,r,dt):
    """Observation only. Never used by transition."""
    s=decode_state(z,layout); sp=decode_state(zprev,layout)
    last=layout.K-1
    q=float(np.real(s[layout.q_var(1),last]))
    qp=[float(np.real(s[layout.q_var(p),last])) for p in range(5)]
    fields=np.empty((2,layout.nx),float); prevfields=np.empty_like(fields)
    for ch in (0,1):
        fields[ch]=[s[layout.field_var(ch,i),last].real for i in range(layout.nx)]
        prevfields[ch]=[sp[layout.field_var(ch,i),last].real for i in range(layout.nx)]
    ft=(fields-prevfields)/dt
    fr=np.gradient(fields,r,axis=1)
    # canonical directional flux proxy in r* = f dr convention: S=-psi_t * f psi_r
    f=1-2/r+q*q/r**2
    S=-ft*(f[None,:]*fr)
    # q-power closure/readout
    qerr=max(abs(qp[p]-q**p) for p in range(5))
    return q,qp,qerr,fields,ft,S


def run_case(q,channel=0,nx=61,K=4,rmin=2.4,rmax=62.4,dt=0.12,steps=700,x0=13.0,sigma=4.0,amp=1e-4,sample_stride=20):
    r=np.linspace(rmin,rmax,nx); layout=Layout(nx=nx,K=K)
    M,pi,pj,nrules=compile_sparse_quadratic(layout,r,dt)
    z=init_state(layout,r,q,x0,sigma,amp,channel)
    zprev=z.copy()  # readout-only snapshot; never enters transition
    # audit initial harmonic state; only z crosses update boundary as physical state.
    flux=np.zeros((2,2),float) # [channel, inner/outer]
    samples=[]; max_qerr=0.; max_nonreal=0.
    for n in range(steps):
        old=z
        z=transition(old,M,pi,pj)
        if not np.isfinite(z).all(): raise FloatingPointError('nonfinite state')
        # observation has no feedback
        if n%sample_stride==0 or n==steps-1:
            qread,qp,qerr,fields,ft,S=readout(z,old,layout,r,dt)
            max_qerr=max(max_qerr,qerr)
            max_nonreal=max(max_nonreal,float(np.max(np.abs(decode_state(z,layout).imag))))
            samples.append({'step':n+1,'q':qread,'q_power_error':qerr,
                            'inner_flux':[float(max(0,-S[ch,1])*dt) for ch in (0,1)],
                            'outer_flux':[float(max(0,S[ch,-2])*dt) for ch in (0,1)]})
        # flux integration is observation-only, never used by z transition
        qread,qp,qerr,fields,ft,S=readout(z,old,layout,r,dt)
        for ch in (0,1):
            flux[ch,0]+=max(0.0,-S[ch,1])*dt
            flux[ch,1]+=max(0.0,S[ch,-2])*dt
        zprev=old
    qread,qp,qerr,fields,ft,S=readout(z,zprev,layout,r,dt)
    # state-channel energies as readout only
    dr=r[1]-r[0]
    energy=[float(0.5*np.sum(ft[ch]**2+( (1-2/r+qread*qread/r**2)*np.gradient(fields[ch],r) )**2)*dr) for ch in (0,1)]
    generated=1-channel
    total_flux=float(np.sum(flux))
    gen_flux=float(np.sum(flux[generated]))
    return {
      'model':'strict harmonic-only single quadratic product-sum RN coupled benchmark',
      'persistent_state':'complex harmonic coefficients only',
      'transition':'one fixed sparse quadratic interaction array; same array for all q cases',
      'generation_channel_names':'anonymous channel 0/1; physical grav/EM interpretation is readout only',
      'parameters':{'q_initial_state':q,'channel_initial_state':channel,'nx':nx,'K':K,'rmin':rmin,'rmax':rmax,'dt':dt,'steps':steps,'x0':x0,'sigma':sigma,'amp':amp},
      'interaction':{'shape':list(M.shape),'nnz':int(M.nnz),'monomial_pairs':int(len(pi)),'sample_rules':int(nrules)},
      'audit':{'q_power_max_error':max_qerr,'max_decoded_imaginary':max_nonreal},
      'flux':{'ch0_inner':float(flux[0,0]),'ch0_outer':float(flux[0,1]),'ch1_inner':float(flux[1,0]),'ch1_outer':float(flux[1,1]),'total':total_flux,
              'generated_other_channel_fraction_of_boundary_flux':float(gen_flux/total_flux if total_flux else 0)},
      'final_energy_proxy':energy,
      'samples':samples
    }


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--q',type=float,required=True); ap.add_argument('--channel',type=int,choices=[0,1],default=0); ap.add_argument('--out',required=True)
    ap.add_argument('--steps',type=int,default=700); ap.add_argument('--nx',type=int,default=61); ap.add_argument('--K',type=int,default=4); ap.add_argument('--x0',type=float,default=13.0)
    a=ap.parse_args(); res=run_case(a.q,a.channel,nx=a.nx,K=a.K,steps=a.steps,x0=a.x0)
    Path(a.out).write_text(json.dumps(res,indent=2,ensure_ascii=False))
    print(json.dumps({k:res[k] for k in ['parameters','interaction','audit','flux','final_energy_proxy']},indent=2,ensure_ascii=False))
if __name__=='__main__': main()
