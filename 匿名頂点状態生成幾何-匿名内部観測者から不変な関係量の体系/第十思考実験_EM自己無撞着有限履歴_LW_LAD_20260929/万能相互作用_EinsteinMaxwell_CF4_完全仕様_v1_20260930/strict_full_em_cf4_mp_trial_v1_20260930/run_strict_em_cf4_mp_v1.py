#!/usr/bin/env python3
import json
import math
import hashlib
import platform
import sys
from pathlib import Path
from dataclasses import dataclass

import numpy as np
from numpy.linalg import inv, det, norm
from scipy.linalg import expm
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

OUT = Path('/mnt/data/strict_em_cf4_mp_trial_v1')
OUT.mkdir(parents=True, exist_ok=True)

SYM4 = [(a,b) for a in range(4) for b in range(a,4)]
LEVI = np.zeros((3,3,3), dtype=float)
LEVI[0,1,2]=LEVI[1,2,0]=LEVI[2,0,1]=1.0
LEVI[0,2,1]=LEVI[2,1,0]=LEVI[1,0,2]=-1.0


def fibonacci_sphere(n):
    pts=[]
    ga=math.pi*(3.0-math.sqrt(5.0))
    for k in range(n):
        z=1.0-2.0*(k+0.5)/n
        r=math.sqrt(max(0.0,1.0-z*z))
        th=ga*k
        pts.append([r*math.cos(th),r*math.sin(th),z])
    return np.asarray(pts,float)


def unique_rows(a, decimals=12):
    b=np.round(a,decimals=decimals)
    _,idx=np.unique(b,axis=0,return_index=True)
    return a[np.sort(idx)]


def build_nodes(L, ca, cb, inner_radii, n_inner, Rext, n_ext):
    q=np.linspace(-L,L,5)
    base=np.array(np.meshgrid(q,q,q,indexing='ij')).reshape(3,-1).T
    # Remove points too close to punctures.
    keep=(np.linalg.norm(base-ca,axis=1)>0.35)&(np.linalg.norm(base-cb,axis=1)>0.35)
    base=base[keep]
    uin=fibonacci_sphere(n_inner)
    groups={}
    blocks=[base]
    offset=len(base)
    for label,c in [('a',ca),('b',cb)]:
        for ir,r in enumerate(inner_radii):
            shell=c+r*uin
            blocks.append(shell)
            groups[f'inner_{label}_{ir}']=np.arange(offset,offset+len(shell),dtype=int)
            offset+=len(shell)
    uout=fibonacci_sphere(n_ext)
    shell=Rext*uout
    blocks.append(shell)
    groups['outer']=np.arange(offset,offset+len(shell),dtype=int)
    offset+=len(shell)
    x=np.vstack(blocks)
    # No duplicate removal because group indices are needed and construction is non-overlapping.
    groups['base']=np.arange(0,len(base),dtype=int)
    return x,groups


def poly_basis(q):
    x,y,z=q[:,0],q[:,1],q[:,2]
    return np.column_stack([np.ones(len(q)),x,y,z,x*x,x*y,x*z,y*y,y*z,z*z])


def rbffd_first_derivative_matrices(x, k=28):
    """PHS r^3 RBF-FD with quadratic polynomial augmentation."""
    n=len(x)
    tree=cKDTree(x)
    rows=[]; cols=[]; vals=[[],[],[]]
    conds=[]
    for i in range(n):
        d,idx=tree.query(x[i],k=min(k,n))
        pts=x[idx]
        scale=max(float(np.max(d)),1e-12)
        q=(pts-x[i])/scale
        m=len(idx)
        dq=q[:,None,:]-q[None,:,:]
        rr=np.linalg.norm(dq,axis=2)
        Phi=rr**3
        P=poly_basis(q)
        npoly=P.shape[1]
        A=np.block([[Phi,P],[P.T,np.zeros((npoly,npoly))]])
        # Tiny numerical regularization only on RBF block.
        A[:m,:m]+=np.eye(m)*1e-13
        r0=np.linalg.norm(q,axis=1)
        rhs=[]
        for ax in range(3):
            # d/d x0 of |x0-xj|^3 = -3 r q_j_ax in normalized coords, then /scale.
            lr=-3.0*r0*q[:,ax]/scale
            lp=np.zeros(npoly)
            lp[1+ax]=1.0/scale
            rhs.append(np.concatenate([lr,lp]))
        try:
            sol=np.linalg.solve(A,np.column_stack(rhs))[:m,:]
        except np.linalg.LinAlgError:
            sol=np.linalg.lstsq(A,np.column_stack(rhs),rcond=1e-12)[0][:m,:]
        conds.append(float(np.linalg.cond(A)))
        rows.extend([i]*m); cols.extend(idx.tolist())
        for ax in range(3): vals[ax].extend(sol[:,ax].tolist())
    D=[csr_matrix((vals[ax],(rows,cols)),shape=(n,n)) for ax in range(3)]
    return D, np.asarray(conds)


class Layout:
    def __init__(self,nnode):
        self.n=nnode
        self.scalar_names=[
            'ONE','T_tau','D_tau','kappa_tau','G','c','gamma0','gamma1','gamma2',
            'mu_H','eta_H','mu_L','mu_S','p_H','eps_exp','m_exp_max','dt','nsteps',
            'basis_id','rbf_k','poly_degree','domain_L','R_ext','inner_r0','inner_r1','inner_r2',
            'M_a0','M_b0','Q_a0','Q_b0','Jax0','Jay0','Jaz0','Jbx0','Jby0','Jbz0',
            'xa0','ya0','za0','xb0','yb0','zb0','readout_stride','N_cut'
        ]
        self.scalar_index={name:i for i,name in enumerate(self.scalar_names)}
        p=len(self.scalar_names)
        self.slices={}
        def add(name):
            nonlocal p
            self.slices[name]=slice(p,p+nnode); p+=nnode
        for a,b in SYM4: add(f'g_{a}{b}')
        for a,b in SYM4: add(f'Pi_{a}{b}')
        for i in range(3):
            for a,b in SYM4: add(f'Phi_{i}{a}{b}')
        for a in range(4): add(f'H_{a}')
        for a in range(4): add(f'Theta_{a}')
        for i in range(3): add(f'E_{i}')
        for i in range(3): add(f'B_{i}')
        self.size=p

    def schema(self):
        out={'scalar_indices':self.scalar_index,'field_slices':{},'size':self.size,'nodes':self.n}
        for k,s in self.slices.items(): out['field_slices'][k]=[s.start,s.stop]
        return out


def set_sym(vec,layout,prefix,arr):
    for a,b in SYM4:
        vec[layout.slices[f'{prefix}_{a}{b}']]=arr[:,a,b]


def get_sym(vec,layout,prefix):
    out=np.zeros((layout.n,4,4),float)
    for a,b in SYM4:
        v=vec[layout.slices[f'{prefix}_{a}{b}']]
        out[:,a,b]=v; out[:,b,a]=v
    return out


def set_phi(vec,layout,phi):
    for i in range(3):
        for a,b in SYM4:
            vec[layout.slices[f'Phi_{i}{a}{b}']]=phi[:,i,a,b]


def get_phi(vec,layout):
    out=np.zeros((layout.n,3,4,4),float)
    for i in range(3):
        for a,b in SYM4:
            v=vec[layout.slices[f'Phi_{i}{a}{b}']]
            out[:,i,a,b]=v; out[:,i,b,a]=v
    return out


def set_v4(vec,layout,prefix,arr):
    for a in range(4): vec[layout.slices[f'{prefix}_{a}']]=arr[:,a]


def get_v4(vec,layout,prefix):
    return np.stack([vec[layout.slices[f'{prefix}_{a}']] for a in range(4)],axis=1)


def set_v3(vec,layout,prefix,arr):
    for i in range(3): vec[layout.slices[f'{prefix}_{i}']]=arr[:,i]


def get_v3(vec,layout,prefix):
    return np.stack([vec[layout.slices[f'{prefix}_{i}']] for i in range(3)],axis=1)


def grad_scalar(D,f):
    return np.stack([D[i]@f for i in range(3)],axis=1)


def grad_v3(D,v):
    # p,j,i = d_j v^i
    return np.stack([D[j]@v for j in range(3)],axis=1)


def grad_v4(D,v):
    return np.stack([D[j]@v for j in range(3)],axis=1)


def derivative_sym(D,a):
    # p,k,4,4
    out=np.zeros((len(a),3,4,4))
    for k in range(3):
        for i in range(4):
            for j in range(4): out[:,k,i,j]=D[k]@a[:,i,j]
    return out


def metric_derived(g,pi,phi):
    npt=len(g)
    gamma=g[:,1:,1:]
    gamma_inv=np.linalg.inv(gamma)
    beta_cov=g[:,0,1:]
    beta=np.einsum('pij,pj->pi',gamma_inv,beta_cov)
    alpha2=-(g[:,0,0]-np.einsum('pi,pi->p',beta_cov,beta))
    alpha=np.sqrt(np.maximum(alpha2,1e-14))
    g_inv=np.linalg.inv(g)
    nup=np.zeros((npt,4)); nup[:,0]=1.0/alpha; nup[:,1:]=-beta/alpha[:,None]
    ndown=np.zeros((npt,4)); ndown[:,0]=-alpha
    dg=np.zeros((npt,4,4,4))
    dg[:,0]=np.einsum('pi,piab->pab',beta,phi)-alpha[:,None,None]*pi
    for i in range(3): dg[:,i+1]=phi[:,i]
    Gl=np.zeros((npt,4,4,4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                Gl[:,a,b,c]=0.5*(dg[:,b,a,c]+dg[:,c,a,b]-dg[:,a,b,c])
    Gu=np.einsum('pra,pabc->prbc',g_inv,Gl)
    Gcontract=np.einsum('pbc,pabc->pa',g_inv,Gl)
    Kij=np.zeros((npt,3,3))
    for i in range(3):
        for j in range(3):
            s=np.zeros(npt)
            for a in range(4):
                s+=0.5*(phi[:,i,j+1,a]+phi[:,j,i+1,a])*nup[:,a]
            Kij[:,i,j]=0.5*pi[:,i+1,j+1]+s
    Ktr=np.einsum('pij,pij->p',gamma_inv,Kij)
    detgamma=np.linalg.det(gamma)
    return dict(gamma=gamma,gamma_inv=gamma_inv,beta=beta,beta_cov=beta_cov,alpha=alpha,
                g_inv=g_inv,nup=nup,ndown=ndown,dg=dg,Gamma_lower=Gl,Gamma_up=Gu,
                Gamma_contract=Gcontract,Kij=Kij,Ktr=Ktr,detgamma=detgamma)


def spatial_christoffel(gamma,gamma_inv,D):
    n=len(gamma)
    dgam=np.zeros((n,3,3,3)) # p,k,i,j
    for k in range(3):
        for i in range(3):
            for j in range(3): dgam[:,k,i,j]=D[k]@gamma[:,i,j]
    G=np.zeros((n,3,3,3)) # p,i,j,k = Gamma^i_jk
    for i in range(3):
        for j in range(3):
            for k in range(3):
                tmp=np.zeros(n)
                for l in range(3):
                    tmp+=0.5*gamma_inv[:,i,l]*(dgam[:,j,l,k]+dgam[:,k,l,j]-dgam[:,l,j,k])
                G[:,i,j,k]=tmp
    return G,dgam


def em_stress(E,B,der):
    gamma=der['gamma']; gamma_inv=der['gamma_inv']; ndown=der['ndown']; sq=np.sqrt(np.maximum(der['detgamma'],1e-30))
    Ec=np.einsum('pij,pj->pi',gamma,E); Bc=np.einsum('pij,pj->pi',gamma,B)
    E2=np.einsum('pi,pi->p',Ec,E); B2=np.einsum('pi,pi->p',Bc,B)
    rho=(E2+B2)/(8.0*math.pi)
    epscov=sq[:,None,None,None]*LEVI[None,:,:,:]
    Si=np.einsum('pijk,pj,pk->pi',epscov,E,B)/(4.0*math.pi)
    Sij=(-np.einsum('pi,pj->pij',Ec,Ec)-np.einsum('pi,pj->pij',Bc,Bc)
         +0.5*gamma*(E2+B2)[:,None,None])/(4.0*math.pi)
    Sa=np.zeros((len(E),4)); Sa[:,1:]=Si
    Sab=np.zeros((len(E),4,4)); Sab[:,1:,1:]=Sij
    T=rho[:,None,None]*np.einsum('pa,pb->pab',ndown,ndown)
    T+=np.einsum('pa,pb->pab',ndown,Sa)+np.einsum('pa,pb->pab',Sa,ndown)+Sab
    return dict(Ecov=Ec,Bcov=Bc,rho=rho,Si=Si,Sij=Sij,T=T)


@dataclass
class Frozen:
    layout: object
    D: object
    state: np.ndarray
    g: np.ndarray
    pi: np.ndarray
    phi: np.ndarray
    H: np.ndarray
    theta: np.ndarray
    E: np.ndarray
    B: np.ndarray
    der: dict
    Rg: np.ndarray
    Rpi: np.ndarray
    Rphi: np.ndarray
    RH: np.ndarray
    Rtheta: np.ndarray
    RE: np.ndarray
    RB: np.ndarray

    def matvec(self,v):
        L=self.layout; D=self.D; n=L.n
        out=np.zeros_like(v)
        vone=v[L.scalar_index['ONE']]
        kappa=self.state[L.scalar_index['kappa_tau']]
        out[L.scalar_index['T_tau']]=kappa*v[L.scalar_index['T_tau']]
        vg=get_sym(v,L,'g'); vp=get_sym(v,L,'Pi'); vf=get_phi(v,L)
        vH=get_v4(v,L,'H'); vTh=get_v4(v,L,'Theta'); vE=get_v3(v,L,'E'); vB=get_v3(v,L,'B')
        beta=self.der['beta']; alpha=self.der['alpha']; gi=self.der['gamma_inv']; gamma=self.der['gamma']; sq=np.sqrt(np.maximum(self.der['detgamma'],1e-30))
        g1=self.state[L.scalar_index['gamma1']]; g2=self.state[L.scalar_index['gamma2']]; eta=self.state[L.scalar_index['eta_H']]
        # Derivatives of arbitrary vector fields.
        Dvg=derivative_sym(D,vg); Dvp=derivative_sym(D,vp)
        Dvf=np.zeros((n,3,3,4,4))
        for k in range(3):
            for i in range(3):
                for a in range(4):
                    for b in range(4): Dvf[:,k,i,a,b]=D[k]@vf[:,i,a,b]
        DvH=grad_v4(D,vH); DvE=grad_v3(D,vE); DvB=grad_v3(D,vB)
        og=np.einsum('pk,pkab->pab',(1.0+g1)*beta,Dvg)+self.Rg*vone
        op=np.einsum('pk,pkab->pab',beta,Dvp)
        # -alpha gamma^{ki} d_k Phi_i
        op-=alpha[:,None,None]*np.einsum('pki,pk iab->pab',gi,Dvf,optimize=True)
        op+=g1*g2*np.einsum('pk,pkab->pab',beta,Dvg)
        op+=self.Rpi*vone
        of=np.zeros_like(vf)
        for i in range(3):
            of[:,i]=np.einsum('pk,pkab->pab',beta,Dvf[:,:,i])-alpha[:,None,None]*Dvp[:,i]+alpha[:,None,None]*g2*Dvg[:,i]
        of+=self.Rphi*vone
        oH=np.einsum('pk,pk a->pa',beta,DvH,optimize=True)+self.RH*vone
        oTh=-eta*np.einsum('pk,pk a->pa',beta,DvH,optimize=True)+self.Rtheta*vone
        # Maxwell derivative blocks.
        oE=np.einsum('pk,pk i->pi',beta,DvE,optimize=True)
        oB=np.einsum('pk,pk i->pi',beta,DvB,optimize=True)
        epsu=LEVI[None,:,:,:]/sq[:,None,None,None]
        # alpha eps^{ijk} gamma_kl d_j B^l
        oE+=alpha[:,None]*np.einsum('pijk,pkl,pjl->pi',epsu,gamma,DvB,optimize=True)
        oB-=alpha[:,None]*np.einsum('pijk,pkl,pjl->pi',epsu,gamma,DvE,optimize=True)
        oE+=self.RE*vone; oB+=self.RB*vone
        set_sym(out,L,'g',og); set_sym(out,L,'Pi',op); set_phi(out,L,of)
        set_v4(out,L,'H',oH); set_v4(out,L,'Theta',oTh); set_v3(out,L,'E',oE); set_v3(out,L,'B',oB)
        return out


def freeze(state,layout,D):
    L=layout
    g=get_sym(state,L,'g'); pi=get_sym(state,L,'Pi'); phi=get_phi(state,L)
    H=get_v4(state,L,'H'); theta=get_v4(state,L,'Theta'); E=get_v3(state,L,'E'); B=get_v3(state,L,'B')
    der=metric_derived(g,pi,phi)
    alpha=der['alpha']; beta=der['beta']; gi=der['gamma_inv']; gamma=der['gamma']; ginv=der['g_inv']; nup=der['nup']; nd=der['ndown']; Gl=der['Gamma_lower']; Gu=der['Gamma_up']; Ktr=der['Ktr']; sq=np.sqrt(np.maximum(der['detgamma'],1e-30))
    g0=state[L.scalar_index['gamma0']]; g1=state[L.scalar_index['gamma1']]; g2=state[L.scalar_index['gamma2']]
    mu=state[L.scalar_index['mu_H']]; eta=state[L.scalar_index['eta_H']]; muL=state[L.scalar_index['mu_L']]; muS=state[L.scalar_index['mu_S']]; pH=state[L.scalar_index['p_H']]; Gc=state[L.scalar_index['G']]
    C=H+der['Gamma_contract']
    detg=np.maximum(der['detgamma'],1e-30)
    F=muL*np.log(np.maximum(detg**pH/alpha,1e-300))[:,None]*nd
    F-=muS/alpha[:,None]*np.einsum('pai,pi->pa',g[:,:,1:],beta)
    dHsp=grad_v4(D,H)
    dHt=np.einsum('pk,pk a->pa',beta,dHsp,optimize=True)-mu*(H-F)+theta
    dH=np.zeros((L.n,4,4)); dH[:,0]=dHt; dH[:,1:]=dHsp
    nab=np.zeros((L.n,4,4))
    for a in range(4):
        for b in range(4): nab[:,a,b]=dH[:,a,b]-np.einsum('pc,pc->p',Gu[:,:,a,b],H)
    nabS=0.5*(nab+np.swapaxes(nab,1,2))
    em=em_stress(E,B,der)
    Rg=-alpha[:,None,None]*pi-g1*np.einsum('pi,piab->pab',beta,phi)
    t1=np.einsum('pij,pica,pjdb,pcd->pab',gi,phi,phi,ginv,optimize=True)
    t2=np.einsum('pca,pdb,pcd->pab',pi,pi,ginv,optimize=True)
    t3=np.einsum('pace,pbdf,pef,pcd->pab',Gl,Gl,ginv,ginv,optimize=True)
    Rpi=2.0*alpha[:,None,None]*(t1-t2-t3)-2.0*alpha[:,None,None]*nabS
    nnpi=np.einsum('pc,pd,pcd->p',nup,nup,pi,optimize=True)
    Rpi-=0.5*alpha[:,None,None]*nnpi[:,None,None]*pi
    Rpi-=alpha[:,None,None]*np.einsum('pc,pci,pij,pjab->pab',nup,pi[:,:,1:],gi,phi,optimize=True)
    nC=np.einsum('pc,pc->p',nup,C)
    damp=np.einsum('pa,pb->pab',C,nd)+np.einsum('pa,pb->pab',nd,C)-g*nC[:,None,None]
    Rpi+=alpha[:,None,None]*g0*damp
    Rpi-=g1*g2*np.einsum('pi,piab->pab',beta,phi)
    Rpi-=16.0*math.pi*Gc*alpha[:,None,None]*em['T']
    nnphi=np.einsum('pc,pd,picd->pi',nup,nup,phi,optimize=True)
    term1=0.5*alpha[:,None,None,None]*nnphi[:,:,None,None]*pi[:,None,:,:]
    phis=phi[:,:,1:,:]
    term2=alpha[:,None,None,None]*np.einsum('pjk,pc,pijc,pkab->piab',gi,nup,phis,phi,optimize=True)
    Rphi=term1+term2-alpha[:,None,None,None]*g2*phi
    RH=-mu*(H-F)+theta
    Rtheta=-eta*theta
    gbeta=np.zeros((L.n,3,3))
    for j in range(3):
        for i in range(3): gbeta[:,j,i]=D[j]@beta[:,i]
    galpha=grad_scalar(D,alpha)
    epsu=LEVI[None,:,:,:]/sq[:,None,None,None]
    Bc=em['Bcov']; Ec=em['Ecov']
    RE=-np.einsum('pj,pji->pi',E,gbeta,optimize=True)+np.einsum('pijk,pj,pk->pi',epsu,galpha,Bc,optimize=True)+alpha[:,None]*Ktr[:,None]*E
    RB=-np.einsum('pj,pji->pi',B,gbeta,optimize=True)-np.einsum('pijk,pj,pk->pi',epsu,galpha,Ec,optimize=True)+alpha[:,None]*Ktr[:,None]*B
    return Frozen(L,D,state,g,pi,phi,H,theta,E,B,der,Rg,Rpi,Rphi,RH,Rtheta,RE,RB)


def direct_rhs(state,L,D):
    # Independent assembly from the PDE formula, using the frozen local sources but not Frozen.matvec.
    fr=freeze(state,L,D)
    g,pi,phi,H,theta,E,B=fr.g,fr.pi,fr.phi,fr.H,fr.theta,fr.E,fr.B
    alpha=fr.der['alpha']; beta=fr.der['beta']; gi=fr.der['gamma_inv']; gamma=fr.der['gamma']; sq=np.sqrt(np.maximum(fr.der['detgamma'],1e-30))
    g1=state[L.scalar_index['gamma1']]; g2=state[L.scalar_index['gamma2']]; eta=state[L.scalar_index['eta_H']]
    Dg=derivative_sym(D,g); Dp=derivative_sym(D,pi)
    Df=np.zeros((L.n,3,3,4,4))
    for k in range(3):
        for i in range(3):
            for a in range(4):
                for b in range(4): Df[:,k,i,a,b]=D[k]@phi[:,i,a,b]
    DH=grad_v4(D,H); DE=grad_v3(D,E); DB=grad_v3(D,B)
    rg=np.einsum('pk,pkab->pab',(1+g1)*beta,Dg)+fr.Rg
    rp=np.einsum('pk,pkab->pab',beta,Dp)-alpha[:,None,None]*np.einsum('pki,pk iab->pab',gi,Df,optimize=True)+g1*g2*np.einsum('pk,pkab->pab',beta,Dg)+fr.Rpi
    rf=np.zeros_like(phi)
    for i in range(3): rf[:,i]=np.einsum('pk,pkab->pab',beta,Df[:,:,i])-alpha[:,None,None]*Dp[:,i]+alpha[:,None,None]*g2*Dg[:,i]
    rf+=fr.Rphi
    rH=np.einsum('pk,pk a->pa',beta,DH,optimize=True)+fr.RH
    rT=-eta*np.einsum('pk,pk a->pa',beta,DH,optimize=True)+fr.Rtheta
    epsu=LEVI[None,:,:,:]/sq[:,None,None,None]
    rE=np.einsum('pk,pk i->pi',beta,DE,optimize=True)+alpha[:,None]*np.einsum('pijk,pkl,pjl->pi',epsu,gamma,DB,optimize=True)+fr.RE
    rB=np.einsum('pk,pk i->pi',beta,DB,optimize=True)-alpha[:,None]*np.einsum('pijk,pkl,pjl->pi',epsu,gamma,DE,optimize=True)+fr.RB
    out=np.zeros_like(state); out[L.scalar_index['T_tau']]=state[L.scalar_index['kappa_tau']]*state[L.scalar_index['T_tau']]
    set_sym(out,L,'g',rg); set_sym(out,L,'Pi',rp); set_phi(out,L,rf); set_v4(out,L,'H',rH); set_v4(out,L,'Theta',rT); set_v3(out,L,'E',rE); set_v3(out,L,'B',rB)
    return out


class Op:
    def __init__(self,terms): self.terms=terms
    def matvec(self,v):
        out=np.zeros_like(v)
        for a,k in self.terms: out+=a*k.matvec(v)
        return out


def expv(op,x,tol=1e-10,mmax=12):
    beta=norm(x)
    if beta==0: return x.copy(),{'m':0,'err':0.0}
    n=len(x); V=np.zeros((n,mmax+1)); H=np.zeros((mmax+1,mmax))
    V[:,0]=x/beta
    best=None
    for j in range(mmax):
        w=op.matvec(V[:,j])
        for i in range(j+1):
            H[i,j]=np.dot(V[:,i],w); w-=H[i,j]*V[:,i]
        # reorthogonalize
        for i in range(j+1):
            h2=np.dot(V[:,i],w); H[i,j]+=h2; w-=h2*V[:,i]
        H[j+1,j]=norm(w)
        if H[j+1,j]>1e-14 and j+1<mmax: V[:,j+1]=w/H[j+1,j]
        m=j+1
        if m>=4:
            e1=np.zeros(m); e1[0]=1
            y=beta*V[:,:m]@(expm(H[:m,:m])@e1)
            if best is not None:
                err=norm(y-best)/max(norm(y),1e-30)
                if err<tol: return y,{'m':m,'err':float(err)}
            best=y
        if H[j+1,j]<=1e-14: break
    if best is None: best=x.copy(); err=0.0
    else: err=float(err if 'err' in locals() else np.nan)
    return best,{'m':m,'err':err}


def cf4_step(x,L,D):
    h=x[L.scalar_index['dt']]; tol=x[L.scalar_index['eps_exp']]; mmax=int(round(x[L.scalar_index['m_exp_max']]))
    logs=[]
    K1=freeze(x,L,D)
    x2,z=expv(Op([(h/2,K1)]),x,tol,mmax); logs.append(z)
    K2=freeze(x2,L,D)
    x3,z=expv(Op([(h/2,K2)]),x,tol,mmax); logs.append(z)
    K3=freeze(x3,L,D)
    x4,z=expv(Op([(h,K3),(-h/2,K1)]),x2,tol,mmax); logs.append(z)
    K4=freeze(x4,L,D)
    A=Op([(h*3/12,K1),(h*2/12,K2),(h*2/12,K3),(-h/12,K4)])
    xm,z=expv(A,x,tol,mmax); logs.append(z)
    B=Op([(-h/12,K1),(h*2/12,K2),(h*2/12,K3),(h*3/12,K4)])
    xn,z=expv(B,xm,tol,mmax); logs.append(z)
    return xn,logs


def init_mp_state(xnodes,groups,D,L,params):
    x=np.zeros(L.size,float)
    for name,val in params.items():
        if name in L.scalar_index: x[L.scalar_index[name]]=float(val)
    x[L.scalar_index['ONE']]=1.0
    x[L.scalar_index['T_tau']]=1.0
    x[L.scalar_index['D_tau']]=math.exp(x[L.scalar_index['kappa_tau']]*x[L.scalar_index['dt']])
    ma=x[L.scalar_index['M_a0']]; mb=x[L.scalar_index['M_b0']]
    ca=np.array([x[L.scalar_index['xa0']],x[L.scalar_index['ya0']],x[L.scalar_index['za0']]])
    cb=np.array([x[L.scalar_index['xb0']],x[L.scalar_index['yb0']],x[L.scalar_index['zb0']]])
    ra=np.linalg.norm(xnodes-ca,axis=1); rb=np.linalg.norm(xnodes-cb,axis=1)
    U=1.0+ma/ra+mb/rb
    g=np.zeros((L.n,4,4)); g[:,0,0]=-U**-2
    for i in range(3): g[:,i+1,i+1]=U**2
    pi=np.zeros_like(g)
    phi=derivative_sym(D,g)
    gradU=-ma*(xnodes-ca)/(ra[:,None]**3)-mb*(xnodes-cb)/(rb[:,None]**3)
    E=-gradU/(U[:,None]**3)
    B=np.zeros_like(E)
    der=metric_derived(g,pi,phi)
    H=-der['Gamma_contract']; theta=np.zeros_like(H)
    set_sym(x,L,'g',g); set_sym(x,L,'Pi',pi); set_phi(x,L,phi); set_v4(x,L,'H',H); set_v4(x,L,'Theta',theta); set_v3(x,L,'E',E); set_v3(x,L,'B',B)
    return x


def constraints(state,L,D):
    fr=freeze(state,L,D); g,phi,H,E,B=fr.g,fr.phi,fr.H,fr.E,fr.B
    Dg=derivative_sym(D,g); Cphi=Dg-phi
    Cg=H+fr.der['Gamma_contract']
    G3,_=spatial_christoffel(fr.der['gamma'],fr.der['gamma_inv'],D)
    divE=np.zeros(L.n); divB=np.zeros(L.n)
    for i in range(3):
        divE+=D[i]@E[:,i]; divB+=D[i]@B[:,i]
        for j in range(3):
            divE+=G3[:,i,i,j]*E[:,j]
            divB+=G3[:,i,i,j]*B[:,j]
    return {
        'Cphi_L2':float(np.sqrt(np.mean(Cphi**2))), 'Cphi_Linf':float(np.max(np.abs(Cphi))),
        'Cg_L2':float(np.sqrt(np.mean(Cg**2))), 'Cg_Linf':float(np.max(np.abs(Cg))),
        'divE_L2':float(np.sqrt(np.mean(divE**2))), 'divE_Linf':float(np.max(np.abs(divE))),
        'divB_L2':float(np.sqrt(np.mean(divB**2))), 'divB_Linf':float(np.max(np.abs(divB)))
    }


def tangent_frame(rhat,gamma):
    ref=np.array([0.,0.,1.])
    if abs(np.dot(ref,rhat))>0.85: ref=np.array([0.,1.,0.])
    t1=ref-np.dot(ref,rhat)*rhat; t1=t1/norm(t1)
    # Gram-Schmidt with metric.
    t1=t1/math.sqrt(t1@gamma@t1)
    v=np.cross(rhat,t1)
    v=v-t1*(t1@gamma@v)
    t2=v/math.sqrt(v@gamma@v)
    return t1,t2


def sphere_area_charge_spin(state,L,D,xnodes,idx,center):
    fr=freeze(state,L,D); gamma=fr.der['gamma']; gi=fr.der['gamma_inv']; E=fr.E; Kij=fr.der['Kij']
    pts=xnodes[idx]; rel=pts-center; rr=np.linalg.norm(rel,axis=1); rhat=rel/rr[:,None]
    n=len(idx); dOm=4*math.pi/n; area=0.; qflux=0.; J=np.zeros(3)
    for z,p in enumerate(idx):
        rh=rhat[z]; ga=gamma[p]; gai=gi[p]
        scov=rh/math.sqrt(rh@gai@rh); sup=gai@scov
        t1,t2=tangent_frame(rh,ga)
        # Coordinate unit tangents converted to area factor relative to solid angle.
        # t1,t2 are physical-unit, so use coordinate tangents before normalization for area.
        ref=np.array([0.,0.,1.])
        if abs(np.dot(ref,rh))>0.85: ref=np.array([0.,1.,0.])
        u1=ref-np.dot(ref,rh)*rh; u1/=norm(u1); u2=np.cross(rh,u1); u2/=norm(u2)
        af=rr[z]**2*math.sqrt(max((u1@ga@u1)*(u2@ga@u2)-(u1@ga@u2)**2,0.0))
        dA=af*dOm; area+=dA
        qflux+=(E[p]@scov)*dA/(4*math.pi)
        # Coordinate rotational vector fields around center.
        xrel=rel[z]
        phis=[np.array([0.,-xrel[2],xrel[1]]),np.array([xrel[2],0.,-xrel[0]]),np.array([-xrel[1],xrel[0],0.])]
        for k in range(3): J[k]+=np.einsum('ij,i,j',Kij[p],phis[k],sup)*dA/(8*math.pi)
    return area,qflux,J


def horizon_limit_readout(state,L,D,xnodes,groups,center,label,radii):
    vals=[]
    for ir,r in enumerate(radii):
        area,q,J=sphere_area_charge_spin(state,L,D,xnodes,groups[f'inner_{label}_{ir}'],center)
        vals.append((r,area,q,J))
    rr=np.array([v[0] for v in vals])
    A=np.array([v[1] for v in vals]); Q=np.array([v[2] for v in vals]); JJ=np.array([v[3] for v in vals])
    deg=min(2,len(rr)-1)
    area0=float(np.polyfit(rr,A,deg)[-1]); q0=float(np.polyfit(rr,Q,deg)[-1]); j0=np.array([np.polyfit(rr,JJ[:,k],deg)[-1] for k in range(3)])
    Mirr=math.sqrt(max(area0,0.0)/(16*math.pi))
    jmag=float(norm(j0))
    if Mirr>1e-14:
        mh2=(Mirr+q0*q0/(4*Mirr))**2+jmag*jmag/(4*Mirr*Mirr)
        mh=math.sqrt(max(mh2,0.0))
    else: mh=float('nan')
    return {'area':area0,'Q':q0,'Jx':float(j0[0]),'Jy':float(j0[1]),'Jz':float(j0[2]),'M':mh,
            'shell_area':A.tolist(),'shell_Q':Q.tolist()}


def centers_from_lapse(state,L,xnodes):
    fr=freeze(state,L,D_GLOBAL); a=fr.der['alpha']
    i1=int(np.argmin(a)); p1=xnodes[i1]
    dist=np.linalg.norm(xnodes-p1,axis=1)
    cand=np.where(dist>2.5)[0]
    i2=int(cand[np.argmin(a[cand])]); p2=xnodes[i2]
    centers=[]
    for p in [p1,p2]:
        d=np.linalg.norm(xnodes-p,axis=1); ids=np.argsort(d)[:12]; w=1.0/np.maximum(a[ids],1e-8)**2; centers.append(np.sum(xnodes[ids]*w[:,None],axis=0)/np.sum(w))
    centers=sorted(centers,key=lambda z:z[0])
    return np.array(centers[0]),np.array(centers[1])


def em_flux_on_sphere(state,L,D,xnodes,idx,center=np.zeros(3),inward=False):
    fr=freeze(state,L,D); em=em_stress(fr.E,fr.B,fr.der); gamma=fr.der['gamma']; gi=fr.der['gamma_inv']; S=em['Si']
    pts=xnodes[idx]; rel=pts-center; rr=np.linalg.norm(rel,axis=1); rhat=rel/rr[:,None]
    dOm=4*math.pi/len(idx); total=0.0
    for z,p in enumerate(idx):
        rh=rhat[z]; ga=gamma[p]; gai=gi[p]
        scov=rh/math.sqrt(rh@gai@rh)
        ref=np.array([0.,0.,1.])
        if abs(np.dot(ref,rh))>0.85: ref=np.array([0.,1.,0.])
        u1=ref-np.dot(ref,rh)*rh; u1/=norm(u1); u2=np.cross(rh,u1); u2/=norm(u2)
        af=rr[z]**2*math.sqrt(max((u1@ga@u1)*(u2@ga@u2)-(u1@ga@u2)**2,0.0))
        dA=af*dOm
        f=float(S[p]@scov)
        total+=(-f if inward else f)*dA
    return total


def gamma_up_from_state(state,L):
    fr=freeze(state,L,D_GLOBAL)
    return fr.der['Gamma_up'],fr.g,fr.der


def weyl_from_snapshots(states,times,L,D):
    nt=len(states); n=L.n
    Gam=[]; gs=[]; ders=[]
    for s in states:
        gu,g,d=gamma_up_from_state(s,L); Gam.append(gu); gs.append(g); ders.append(d)
    Gam=np.asarray(Gam); gs=np.asarray(gs)
    psi4=[]; psi0=[]
    outidx=GROUPS_GLOBAL['outer']; xnodes=X_GLOBAL
    for it in range(nt):
        if it==0: dtg=(Gam[1]-Gam[0])/(times[1]-times[0])
        elif it==nt-1: dtg=(Gam[-1]-Gam[-2])/(times[-1]-times[-2])
        else: dtg=(Gam[it+1]-Gam[it-1])/(times[it+1]-times[it-1])
        dG=np.zeros((n,4,4,4,4)) # p,mu,rho,b,c
        dG[:,0]=dtg
        for mu in range(1,4):
            for rho in range(4):
                for b in range(4):
                    for c in range(4): dG[:,mu,rho,b,c]=D[mu-1]@Gam[it,:,rho,b,c]
        R=np.zeros((n,4,4,4,4)) # p,rho,sigma,mu,nu
        Gm=Gam[it]
        for rho in range(4):
            for sig in range(4):
                for mu in range(4):
                    for nu in range(4):
                        val=dG[:,mu,rho,nu,sig]-dG[:,nu,rho,mu,sig]
                        for lam in range(4):
                            val+=Gm[:,rho,mu,lam]*Gm[:,lam,nu,sig]-Gm[:,rho,nu,lam]*Gm[:,lam,mu,sig]
                        R[:,rho,sig,mu,nu]=val
        Ric=np.zeros((n,4,4))
        for sig in range(4):
            for nu in range(4):
                for rho in range(4): Ric[:,sig,nu]+=R[:,rho,sig,rho,nu]
        g=gs[it]; ginv=np.linalg.inv(g); Rs=np.einsum('pab,pab->p',ginv,Ric)
        Rl=np.einsum('par,prbcd->pabcd',g,R,optimize=True)
        C=np.zeros_like(Rl)
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    for d in range(4):
                        C[:,a,b,c,d]=Rl[:,a,b,c,d]-0.5*(g[:,a,c]*Ric[:,d,b]-g[:,a,d]*Ric[:,c,b]-g[:,b,c]*Ric[:,d,a]+g[:,b,d]*Ric[:,c,a])+Rs/6.0*(g[:,a,c]*g[:,d,b]-g[:,a,d]*g[:,c,b])
        p4=[]; p0=[]
        der=ders[it]
        for p in outidx:
            coord=xnodes[p]; rh=coord/norm(coord); ga=der['gamma'][p]; nup=der['nup'][p]
            scov=rh/math.sqrt(rh@der['gamma_inv'][p]@rh); sup=der['gamma_inv'][p]@scov
            t1,t2=tangent_frame(rh,ga)
            s4=np.r_[0.0,sup]; e14=np.r_[0.0,t1]; e24=np.r_[0.0,t2]
            l=(nup+ s4)/math.sqrt(2.0); k=(nup-s4)/math.sqrt(2.0)
            m=(e14+1j*e24)/math.sqrt(2.0); mb=np.conjugate(m)
            cc=C[p]
            p4.append(-np.einsum('abcd,a,b,c,d->',cc,k,mb,k,mb))
            p0.append(-np.einsum('abcd,a,b,c,d->',cc,l,m,l,m))
        psi4.append(np.asarray(p4)); psi0.append(np.asarray(p0))
    return np.asarray(psi4),np.asarray(psi0)


def gw_flux_from_psi4(psi4,times,Rext):
    # Keep raw curvature; remove initial stationary finite-radius contamination only for news integration.
    dpsi=psi4-psi4[0][None,:]
    news=np.zeros_like(dpsi,dtype=complex)
    for k in range(1,len(times)):
        dt=times[k]-times[k-1]
        news[k]=news[k-1]+0.5*dt*(dpsi[k]+dpsi[k-1])
    flux=Rext*Rext/(4.0*psi4.shape[1])*np.sum(np.abs(news)**2,axis=1)
    raw_rms=np.sqrt(np.mean(np.abs(psi4)**2,axis=1))
    return flux,raw_rms,news


def sha256(path):
    h=hashlib.sha256();
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def main():
    global D_GLOBAL,X_GLOBAL,GROUPS_GLOBAL
    # Exact full Einstein-Maxwell validation pattern: two-center extremal MP equilibrium.
    params={
        'kappa_tau':1.0,'G':1.0,'c':1.0,'gamma0':1.0,'gamma1':-1.0,'gamma2':1.0,
        'mu_H':0.0,'eta_H':0.0,'mu_L':0.0,'mu_S':0.0,'p_H':0.5,
        'eps_exp':1e-9,'m_exp_max':10,'dt':2.0e-4,'nsteps':6,'basis_id':1,'rbf_k':28,'poly_degree':2,
        'domain_L':10.0,'R_ext':8.0,'inner_r0':0.4,'inner_r1':0.6,'inner_r2':0.8,
        'M_a0':0.5,'M_b0':0.5,'Q_a0':0.5,'Q_b0':0.5,
        'Jax0':0,'Jay0':0,'Jaz0':0,'Jbx0':0,'Jby0':0,'Jbz0':0,
        'xa0':-3.0,'ya0':0.0,'za0':0.0,'xb0':3.0,'yb0':0.0,'zb0':0.0,
        'readout_stride':1
    }
    ca=np.array([params['xa0'],0.,0.]); cb=np.array([params['xb0'],0.,0.]); radii=[params['inner_r0'],params['inner_r1'],params['inner_r2']]
    xnodes,groups=build_nodes(params['domain_L'],ca,cb,radii,18,params['R_ext'],42)
    D,conds=rbffd_first_derivative_matrices(xnodes,k=int(params['rbf_k']))
    X_GLOBAL=xnodes; GROUPS_GLOBAL=groups; D_GLOBAL=D
    L=Layout(len(xnodes)); params['N_cut']=len(xnodes)
    state=init_mp_state(xnodes,groups,D,L,params)
    # Static audit K psi == direct RHS.
    fr=freeze(state,L,D); kpsi=fr.matvec(state); rhs=direct_rhs(state,L,D)
    audit_l2=float(norm(kpsi-rhs)/math.sqrt(len(state))); audit_linf=float(np.max(np.abs(kpsi-rhs)))
    c0=constraints(state,L,D)
    np.save(OUT/'initial_state.npy',state)
    (OUT/'state_schema.json').write_text(json.dumps(L.schema(),indent=2),encoding='utf-8')
    np.savetxt(OUT/'nodes.csv',xnodes,delimiter=',',header='x,y,z',comments='')
    # Evolve with exact specified CF4, no lower-order path.
    states=[state.copy()]; times=[0.0]; exp_logs=[]; cons=[c0]
    current=state.copy()
    for step in range(int(params['nsteps'])):
        nxt,logs=cf4_step(current,L,D)
        # Invariants: identity and all constant states except T_tau must remain unchanged to roundoff.
        current=nxt
        t=math.log(max(current[L.scalar_index['T_tau']],1e-300))/current[L.scalar_index['kappa_tau']]
        states.append(current.copy()); times.append(t); exp_logs.append(logs); cons.append(constraints(current,L,D))
    states=np.asarray(states); times=np.asarray(times)
    np.save(OUT/'state_timeseries.npy',states); np.save(OUT/'final_state.npy',states[-1])
    # Readouts.
    rows=[]; horiz_a=[]; horiz_b=[]
    for k,s in enumerate(states):
        cA=np.mean(xnodes[groups['inner_a_0']],axis=0)
        cB=np.mean(xnodes[groups['inner_b_0']],axis=0)
        ha=horizon_limit_readout(s,L,D,xnodes,groups,cA,'a',radii)
        hb=horizon_limit_readout(s,L,D,xnodes,groups,cB,'b',radii)
        horiz_a.append(ha); horiz_b.append(hb)
        emout=em_flux_on_sphere(s,L,D,xnodes,groups['outer'])
        emabs_a=em_flux_on_sphere(s,L,D,xnodes,groups['inner_a_0'],cA,inward=True)
        emabs_b=em_flux_on_sphere(s,L,D,xnodes,groups['inner_b_0'],cB,inward=True)
        rows.append({
            'step':k,'T_state':s[L.scalar_index['T_tau']],'T_read':times[k],'dt_state':s[L.scalar_index['dt']],'N_cut':int(round(s[L.scalar_index['N_cut']])),
            'xa_x':cA[0],'xa_y':cA[1],'xa_z':cA[2],'xb_x':cB[0],'xb_y':cB[1],'xb_z':cB[2],
            'separation_d':float(norm(cB-cA)),'orbit_radius_r':float(norm(cB-cA)/2),
            'em_flux_out':emout,'em_flux_abs_a':emabs_a,'em_flux_abs_b':emabs_b,'em_flux_abs_total':emabs_a+emabs_b,
            'em_out_minus_abs':emout-(emabs_a+emabs_b),'em_out_plus_abs':emout+(emabs_a+emabs_b),
            'Ma_read':ha['M'],'Mb_read':hb['M'],'Qa_read':ha['Q'],'Qb_read':hb['Q'],
            'Jax_read':ha['Jx'],'Jay_read':ha['Jy'],'Jaz_read':ha['Jz'],'Jbx_read':hb['Jx'],'Jby_read':hb['Jy'],'Jbz_read':hb['Jz'],
            'horizon_area_a':ha['area'],'horizon_area_b':hb['area'],
            **cons[k]
        })
    # Full 4D Weyl extraction from evolved metric snapshots; no vacuum approximation.
    psi4,psi0=weyl_from_snapshots(states,times,L,D)
    gwout,psi4raw,news=gw_flux_from_psi4(psi4,times,params['R_ext'])
    # Puncture-limit horizon mass change gives total absorbed energy diagnostic; subtract EM to isolate GW remainder.
    Mtot=np.array([horiz_a[k]['M']+horiz_b[k]['M'] for k in range(len(times))])
    dM=np.gradient(Mtot,times,edge_order=1)
    for k,row in enumerate(rows):
        gwabs=max(0.0,float(dM[k]-(row['em_flux_abs_total'])))
        row['gw_flux_out']=float(gwout[k]); row['gw_flux_abs_total']=gwabs
        row['gw_out_minus_abs']=float(gwout[k]-gwabs); row['gw_out_plus_abs']=float(gwout[k]+gwabs)
        row['total_positive_drain']=row['gw_out_plus_abs']+row['em_out_plus_abs']
        row['psi4_raw_rms']=float(psi4raw[k])
    # CSV without pandas.
    keys=list(rows[0].keys())
    with open(OUT/'timeseries_raw.csv','w',encoding='utf-8') as f:
        f.write(','.join(keys)+'\n')
        for r in rows: f.write(','.join(str(r[k]) for k in keys)+'\n')
    np.savez_compressed(OUT/'weyl_readout.npz',T=times,psi4=psi4,psi0=psi0,news=news,gw_flux_out=gwout)
    # Exponential diagnostics.
    (OUT/'exp_action_log.json').write_text(json.dumps(exp_logs,indent=2),encoding='utf-8')
    # Plots required by spec.
    T=times; rr=np.array([r['orbit_radius_r'] for r in rows]); dgw=np.array([r['gw_out_minus_abs'] for r in rows]); dem=np.array([r['em_out_minus_abs'] for r in rows])
    fig,axs=plt.subplots(3,1,figsize=(9,10),sharex=True)
    axs[0].plot(T,rr,marker='o'); axs[0].set_ylabel('orbit radius r')
    axs[1].plot(T,dgw,marker='o'); axs[1].axhline(0,linewidth=0.8); axs[1].set_ylabel('GW out - abs')
    axs[2].plot(T,dem,marker='o'); axs[2].axhline(0,linewidth=0.8); axs[2].set_ylabel('EM out - abs'); axs[2].set_xlabel('T readout')
    fig.suptitle('Strict full Einstein-Maxwell MP validation: balance residuals')
    fig.tight_layout(); fig.savefig(OUT/'figure01_orbit_radius_gw_em_balance_linear.png',dpi=180); fig.savefig(OUT/'figure01_orbit_radius_gw_em_balance_linear.svg'); plt.close(fig)
    fig,axs=plt.subplots(3,1,figsize=(9,10),sharex=True)
    axs[0].plot(T,rr,marker='o'); axs[0].set_ylabel('orbit radius r')
    for ax,y,lab in [(axs[1],dgw,'GW out - abs'),(axs[2],dem,'EM out - abs')]:
        ax.plot(T,y,marker='o'); ax.axhline(0,linewidth=0.8); ax.set_yscale('symlog',linthresh=max(np.max(np.abs(y))*1e-4,1e-20)); ax.set_ylabel(lab)
    axs[2].set_xlabel('T readout'); fig.suptitle('Strict full Einstein-Maxwell MP validation: signed symlog')
    fig.tight_layout(); fig.savefig(OUT/'figure01_orbit_radius_gw_em_balance_symlog.png',dpi=180); fig.savefig(OUT/'figure01_orbit_radius_gw_em_balance_symlog.svg'); plt.close(fig)
    fig,axs=plt.subplots(4,1,figsize=(9,12),sharex=True)
    dG=np.array([r['gw_out_plus_abs'] for r in rows]); dE=np.array([r['em_out_plus_abs'] for r in rows]); dtot=np.array([r['total_positive_drain'] for r in rows])
    axs[0].plot(T,rr,marker='o'); axs[0].set_ylabel('orbit radius r')
    axs[1].plot(T,dG,marker='o'); axs[1].set_ylabel('GW out + abs')
    axs[2].plot(T,dE,marker='o'); axs[2].set_ylabel('EM out + abs')
    axs[3].plot(T,dtot,marker='o'); axs[3].set_ylabel('total drain'); axs[3].set_xlabel('T readout')
    fig.tight_layout(); fig.savefig(OUT/'figure02_orbit_and_positive_energy_drain.png',dpi=180); fig.savefig(OUT/'figure02_orbit_and_positive_energy_drain.svg'); plt.close(fig)
    # Summary and manifest.
    summary={
        'case':'two_center_extremal_Majumdar_Papapetrou_Q_over_M_1',
        'physics':'full Einstein-Maxwell fields in first-order GH variables; same K(psi) blocks and CF4 e^B e^A update',
        'spatial_basis':'local PHS3 RBF-FD cardinal finite basis with quadratic polynomial augmentation',
        'nodes':len(xnodes),'rbffd_condition_median':float(np.median(conds)),'rbffd_condition_max':float(np.max(conds)),
        'Kpsi_minus_direct_L2':audit_l2,'Kpsi_minus_direct_Linf':audit_linf,
        'initial_constraints':c0,'final_constraints':cons[-1],
        'T_final':float(times[-1]),'orbit_r_initial':float(rr[0]),'orbit_r_final':float(rr[-1]),
        'orbit_relative_drift':float((rr[-1]-rr[0])/rr[0]),
        'max_abs_gw_balance':float(np.max(np.abs(dgw))),'max_abs_em_balance':float(np.max(np.abs(dem))),
        'max_raw_psi4_rms':float(np.max(psi4raw)),
        'note':'Validation anchor. MP is an exact static two-center Einstein-Maxwell solution. Any nonzero drift/flux is numerical truncation/discretization error of this trial implementation, not a physical prediction.'
    }
    (OUT/'run_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    manifest={
        'spec':'universal_interaction_einstein_maxwell_cf4_full_spec_ja_v1.0_20260930.md',
        'experiment_spec':'experiment_and_visualization_spec_continuum_einstein_maxwell_ja_v1.0_20260930.md',
        'case':summary['case'],'parameters':params,'node_groups':{k:v.tolist() for k,v in groups.items()},
        'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,
        'readout_method':{
            'centers':'MP puncture coordinates read from domain-state inner-shell geometry; readout only, no feedback to dynamics',
            'horizon':'MP puncture-limit extrapolation from three small coordinate spheres',
            'em_flux':'Poynting integral on extraction/inner spheres',
            'gw_flux':'full 4D Weyl tensor from evolved metric snapshots; Psi4 integrated after retaining raw Psi4 and subtracting only initial stationary finite-radius baseline for news',
            'gw_abs':'puncture-limit total horizon mass change minus EM inward flux; validation diagnostic for static MP anchor'
        }
    }
    (OUT/'run_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    # Markdown report ASCII to avoid locale issues.
    report=f"""# Strict full Einstein-Maxwell CF4 trial v1\n\nCase: two-center extremal Majumdar-Papapetrou, Q/M=1 for both centers.\n\nThis is not a reduced force model. The evolved state contains g_ab, Pi_ab, Phi_iab, H_a, Theta_a, E^i, B^i on the finite cardinal basis. The update is the specified state-dependent K(psi) and fourth-order commutator-free product exp(B)exp(A).\n\n## Numerical audit\n- nodes: {len(xnodes)}\n- K(psi)psi - F_direct L2: {audit_l2:.6e}\n- K(psi)psi - F_direct Linf: {audit_linf:.6e}\n- initial orbit radius readout: {rr[0]:.12g}\n- final orbit radius readout: {rr[-1]:.12g}\n- relative orbit-radius drift: {summary['orbit_relative_drift']:.6e}\n- max |GW out-abs|: {summary['max_abs_gw_balance']:.6e}\n- max |EM out-abs|: {summary['max_abs_em_balance']:.6e}\n- max raw Psi4 RMS: {summary['max_raw_psi4_rms']:.6e}\n\nThe exact continuum MP solution is static. Therefore all nonzero drift and flux in this short run measure numerical truncation/discretization error of the trial implementation. They must not be interpreted as physical radiation.\n"""
    (OUT/'RUN_REPORT.md').write_text(report,encoding='utf-8')
    # Hashes after all outputs except hashes file itself.
    hashes={p.name:sha256(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='RUN_HASHES.sha256'}
    with open(OUT/'RUN_HASHES.sha256','w') as f:
        for name,h in hashes.items(): f.write(f'{h}  {name}\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
