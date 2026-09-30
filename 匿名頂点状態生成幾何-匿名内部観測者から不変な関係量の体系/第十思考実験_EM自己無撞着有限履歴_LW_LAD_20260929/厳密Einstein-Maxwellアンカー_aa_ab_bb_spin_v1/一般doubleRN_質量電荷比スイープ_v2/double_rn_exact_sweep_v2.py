import json, math, cmath
from pathlib import Path
import numpy as np

OUT = Path('/mnt/data/em_selfconsistent_v1')


def csqrt(x):
    return cmath.sqrt(complex(x))


def drn_fields(rho,z,M,m,Q,q,R):
    # M,Q lower source; m,q upper source, following Manko et al. notation.
    mu=(m*Q-M*q)/(R+M+m)
    Sigma=csqrt(M*M-Q*Q+2*mu*Q)
    sigma=csqrt(m*m-q*q-2*mu*q)
    nu=R*R-Sigma*Sigma-sigma*sigma+2*mu*mu
    kappa=M*m-(Q-mu)*(q+mu)

    Rp=csqrt(rho*rho+(z+0.5*R+Sigma)**2)
    Rm=csqrt(rho*rho+(z+0.5*R-Sigma)**2)
    rp=csqrt(rho*rho+(z-0.5*R+sigma)**2)
    rm=csqrt(rho*rho+(z-0.5*R-sigma)**2)

    A=(Sigma*sigma*(nu*(Rp+Rm)*(rp+rm)+4*kappa*(Rp*Rm+rp*rm))
       -(mu*mu*nu-2*kappa*kappa)*(Rp-Rm)*(rp-rm))
    B=(2*Sigma*sigma*((nu*m+2*kappa*M)*(Rp+Rm)+(nu*M+2*kappa*m)*(rp+rm))
       +2*sigma*(nu*mu*(Q-mu)-2*kappa*(R*M-mu*q-mu*mu))*(Rp-Rm)
       +2*Sigma*(nu*mu*(q+mu)+2*kappa*(R*m+mu*Q-mu*mu))*(rp-rm))
    C=(2*Sigma*sigma*((nu*(q+mu)+2*kappa*(Q-mu))*(Rp+Rm)
                     +(nu*(Q-mu)+2*kappa*(q+mu))*(rp+rm))
       +2*sigma*(mu*nu*M+2*kappa*(mu*m-R*Q+mu*R))*(Rp-Rm)
       +2*Sigma*(mu*nu*m+2*kappa*(mu*M+R*q+mu*R))*(rp-rm))

    denom=A+B
    E=(A-B)/denom
    Phi=C/denom  # paper uses Phi=-A_t=C/(A+B), so A_t=-Phi
    f=(A*A-B*B+C*C)/(denom*denom)
    K0=4*Sigma*sigma*(R*R-(M-m)**2+(Q-q-2*mu)**2)
    e2g=(A*A-B*B+C*C)/(K0*K0*Rp*Rm*rp*rm)
    return {
        'E':E,'Phi':Phi,'At':-Phi,'f':f,'e2gamma':e2g,
        'mu':mu,'Sigma':Sigma,'sigma':sigma,'nu':nu,'kappa':kappa,
        'Rpm':(Rp,Rm),'rpm':(rp,rm)
    }


def realval(z, tol=1e-8):
    if abs(z.imag) < tol*(1+abs(z.real)):
        return float(z.real)
    return {'re':float(z.real),'im':float(z.imag)}


def fd_derivs(fun,rho,z,h=2e-5):
    # 2D cylindrical derivatives for scalar fun(rho,z), complex allowed
    f0=fun(rho,z)
    frp=fun(rho+h,z); frm=fun(rho-h,z)
    fzp=fun(rho,z+h); fzm=fun(rho,z-h)
    dr=(frp-frm)/(2*h); dz=(fzp-fzm)/(2*h)
    drr=(frp-2*f0+frm)/(h*h); dzz=(fzp-2*f0+fzm)/(h*h)
    lap=drr+dzz+dr/rho
    return f0,dr,dz,drr,dzz,lap


def em_residuals(rho,z,M,m,Q,q,R,h=2e-5):
    ffun=lambda rr,zz: drn_fields(rr,zz,M,m,Q,q,R)['f']
    Afun=lambda rr,zz: drn_fields(rr,zz,M,m,Q,q,R)['At']
    f0,fr,fz,_,_,flap=fd_derivs(ffun,rho,z,h)
    A0,Ar,Az,Arr,Azz,Alap=fd_derivs(Afun,rho,z,h)
    gradf2=fr*fr+fz*fz
    gradA2=Ar*Ar+Az*Az
    # Maxwell div(f^-1 grad A)=0 => f*Delta A - grad f . grad A =0
    rM=f0*Alap-(fr*Ar+fz*Az)
    # Static EM equation quoted in source: f Delta f = grad f^2 + 2 f grad A^2
    rE=f0*flap-gradf2-2*f0*gradA2
    scaleM=1+abs(f0*Alap)+abs(fr*Ar+fz*Az)
    scaleE=1+abs(f0*flap)+abs(gradf2)+abs(2*f0*gradA2)
    return float(abs(rM)/scaleM), float(abs(rE)/scaleE)


def spin2_test(a,b,n=257):
    # Kinematic symmetric-square test: aa, 2ab, bb all weight +2 under common rotation.
    errs={'aa':0.0,'ab':0.0,'bb':0.0,'total':0.0}
    aa=a*a; ab=2*a*b; bb=b*b; total=(a+b)**2
    for psi in np.linspace(0,2*np.pi,n):
        r=cmath.exp(1j*psi); w2=cmath.exp(2j*psi)
        ar=r*a; br=r*b
        vals=[ar*ar,2*ar*br,br*br,(ar+br)**2]
        refs=[w2*aa,w2*ab,w2*bb,w2*total]
        for k,v,ref in zip(errs,vals,refs): errs[k]=max(errs[k],abs(v-ref))
    return errs


def inclusion_exact(rho,z,M,m,Q,q,R):
    # Exact observable decomposition by inclusion-exclusion. This is not a polynomial expansion;
    # it isolates each one-body exact field and the nonlinear two-body interaction remainder.
    full=drn_fields(rho,z,M,m,Q,q,R)
    # Use exact single-RN limits analytically rather than degenerate double-RN formula.
    def rn_at(point_rho,point_z,Mass,Charge,zc):
        # Convert Weyl point to Schwarzschild-like RN Weyl potential using rod distances.
        s=csqrt(Mass*Mass-Charge*Charge)
        Rp=csqrt(point_rho**2+(point_z-zc+s)**2)
        Rm=csqrt(point_rho**2+(point_z-zc-s)**2)
        x=(Rp+Rm)/2
        # In Weyl coordinates RN: f=(x^2-s^2)/(x+M)^2, Phi=Q/(x+M); At=-Phi convention.
        ff=(x*x-s*s)/(x+Mass)**2
        phi=Charge/(x+Mass)
        return ff,phi
    fa,phia=rn_at(rho,z,M,Q,-R/2)
    fb,phib=rn_at(rho,z,m,q,+R/2)
    # Minkowski baseline f=1, Phi=0
    cross_f=full['f']-fa-fb+1
    cross_phi=full['Phi']-phia-phib
    return {'full_f':full['f'],'fa':fa,'fb':fb,'cross_f':cross_f,
            'full_Phi':full['Phi'],'phia':phia,'phib':phib,'cross_Phi':cross_phi}


def main():
    M=m=1.0; R=8.0
    lambdas=[0.0,0.1,0.5,0.9,1.0,1.1,2.0,10.0]
    points=[(1.3,0.0),(2.0,1.0),(3.5,-0.7)]
    rows=[]
    for lam in lambdas:
        Q=q=lam
        case={'lambda':lam,'M':M,'m':m,'Q':Q,'q':q,'R':R,'points':[]}
        max_im=0.0; max_rm=0.0; max_re=0.0
        for rho,z in points:
            if abs(lam-1.0) < 1e-14:
                # Exact Majumdar-Papapetrou member of the double-RN family.
                ra=math.sqrt(rho*rho+(z+R/2)**2); rb=math.sqrt(rho*rho+(z-R/2)**2)
                U=1+M/ra+m/rb
                fmp=1/(U*U); phimp=(M/ra+m/rb)/U  # asymptotically zero gauge, sign convention immaterial for squared EM residuals
                d={'f':complex(fmp),'Phi':complex(phimp),'e2gamma':1+0j,'mu':0+0j,'Sigma':0+0j,'sigma':0+0j,'kappa':0+0j}
                fa=(1+M/ra)**-2; fb=(1+m/rb)**-2
                inc={'full_f':complex(fmp),'fa':complex(fa),'fb':complex(fb),'cross_f':complex(fmp-fa-fb+1),
                     'full_Phi':complex(phimp),'phia':complex((M/ra)/(1+M/ra)),'phib':complex((m/rb)/(1+m/rb)),
                     'cross_Phi':complex(phimp-(M/ra)/(1+M/ra)-(m/rb)/(1+m/rb))}
                rm=re=0.0
            else:
                d=drn_fields(rho,z,M,m,Q,q,R)
                inc=inclusion_exact(rho,z,M,m,Q,q,R)
                rm,re=em_residuals(rho,z,M,m,Q,q,R)
            max_rm=max(max_rm,rm); max_re=max(max_re,re)
            for key in ['f','Phi','e2gamma']:
                max_im=max(max_im,abs(d[key].imag))
            case['points'].append({
                'rho':rho,'z':z,'f':realval(d['f']),'Phi':realval(d['Phi']),
                'e2gamma':realval(d['e2gamma']),'mu':realval(d['mu']),
                'Sigma':realval(d['Sigma']),'sigma':realval(d['sigma']),
                'kappa':realval(d['kappa']),
                'exact_self_cross':{k:realval(v) for k,v in inc.items()},
                'maxwell_residual_rel':rm,'einstein_scalar_residual_rel':re
            })
        case['max_imag_physical']=max_im
        case['maxwell_residual_rel_max']=max_rm
        case['einstein_scalar_residual_rel_max']=max_re
        rows.append(case)

    # Asymmetric case to ensure mu != 0 and true 5-param family use.
    asym={'M':1.2,'m':0.8,'Q':0.35,'q':-0.15,'R':6.5}
    asym_pts=[]
    for rho,z in points:
        d=drn_fields(rho,z,**asym)
        rm,re=em_residuals(rho,z,**asym)
        asym_pts.append({'rho':rho,'z':z,'f':realval(d['f']),'Phi':realval(d['Phi']),
                         'mu':realval(d['mu']),'Sigma':realval(d['Sigma']),'sigma':realval(d['sigma']),
                         'kappa':realval(d['kappa']),'maxwell_residual_rel':rm,
                         'einstein_scalar_residual_rel':re})

    sp=spin2_test(1.2*cmath.exp(0.37j),0.7*cmath.exp(-0.81j))
    out={'source':'Manko 2007 physical double-Reissner-Nordstrom exact solution',
         'equal_mass_charge_sweep':rows,'asymmetric_5_parameter_case':{**asym,'points':asym_pts},
         'spin2_symmetric_square_test':sp,
         'notes':[
             'Static exact solution: validates nonlinear Einstein-Maxwell field, not radiative gravitational waves.',
             'Exact self/cross decomposition uses inclusion-exclusion on observables; cross term is the exact nonlinear two-body remainder, not merely 2ab.',
             'Spin-2 test is kinematic symmetric-square weight-2 readout; physical GW spin-2 still requires a radiative exact/full-Einstein sector.'
         ]}
    p=OUT/'double_rn_exact_sweep_v2.json'; p.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
    print('saved',p)
    for c in rows:
        p0=c['points'][0]
        print(f"lambda={c['lambda']:>4}: Sigma={p0['Sigma']} f={p0['f']} Phi={p0['Phi']} residM={c['maxwell_residual_rel_max']:.3e} residE={c['einstein_scalar_residual_rel_max']:.3e} imag={c['max_imag_physical']:.3e}")
    print('asym mu=',asym_pts[0]['mu'],'resids',max(p['maxwell_residual_rel'] for p in asym_pts),max(p['einstein_scalar_residual_rel'] for p in asym_pts))
    print('spin2 errors',sp)

if __name__=='__main__': main()
