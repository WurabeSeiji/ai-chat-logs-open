# coefficients_relational_full_list_ja_v3_20261002.md の検算スクリプト（2026-10-02）
# 第1部: 各項の値と円軌道・α の確認   第2部: 単項式展開（22 行）
# 実行: python3 coefficients_relational_full_list_v3_check_20261002.py   （mpmath, sympy）

# ===== 第1部 =====
# Numerical check of the relational Hamiltonian coefficient list (v3).
# Units G = q0 = c = 1, Gauss, Q_a=-3, Q_b=+3, J_a=J_b=9/(2 alpha).
from mpmath import mp, mpf, findroot, diff, log, pi
mp.dps = 30
alpha = mpf('7.2973525643e-3')
ma = mpf('7.0e-4'); ratio = mpf('1836.15267343'); mb = ma*ratio
M = ma+mb; mu = ma*mb/M; nu = mu/M
Ja = 9/(2*alpha); aa = Ja/ma            # ring radius of a
kC = 9/(mu*aa); kM = M/aa; kG = ma*mb/(mu*aa)
C0m1 = 9/(ma*mb)
print("kC=",kC," kM=",kM," kG=",kG," kG==kM:",abs(kG-kM)<1e-25," C0-1=",C0m1," nu=",nu)
print("check kC/(C0-1)=",kC/C0m1)
def vm(m,xi): return 2*m*(ma/mu)/xi
# Hamiltonian per mu, circular (pi=0): P^2 = vm^2, (n.P)=0
def eps_terms(m,xi,sa=1,sb=1):
    v2 = vm(m,xi)**2; P2=v2; pi_=mpf(0)
    T={}
    T['rot kinetic 1/2 v_m^2']        = P2/2
    T['Coulomb -kC/xi']               = -kC/xi
    T['Newton -kG/xi']                = -kG/xi
    T['SR kinetic (3nu-1)P^4/8']      = (3*nu-1)*P2**2/8
    T['Darwin -nu kC (P^2+pi^2)/(2xi)']= -nu*kC*(P2+pi_**2)/(2*xi)
    T['EIH -kM[(3+nu)P^2+nu pi^2]/(2xi)'] = -kM*((3+nu)*P2+nu*pi_**2)/(2*xi)
    T['EIH +kM^2/(2xi^2)']            = kM**2/(2*xi**2)
    T['mixed static 3 kC kM/(2 xi^2)']= 3*kC*kM/(2*xi**2)
    T['EM2PC -kC/(8xi)[3nu^2 pi^4+nu(3nu-2)P^4+2nu(nu-1)P^2pi^2]'] = -kC/(8*xi)*(3*nu**2*pi_**4+nu*(3*nu-2)*P2**2+2*nu*(nu-1)*P2*pi_**2)
    T['EM2PC (1-5nu+5nu^2)P^6/16']    = (1-5*nu+5*nu**2)*P2**3/16
    T['EM2PC nu kC^2 P^2/(4xi^2)']    = nu*kC**2*P2/(4*xi**2)
    T['EM2PC nu kC^3/(4xi^3)']        = nu*kC**3/(4*xi**3)
    T['SO_a kC(1+2ma/mb) m sa/xi^3']  = kC*(1+2*ma/mb)*m*sa/xi**3
    T['SO_b kC(ma/mb)^2(1+2mb/ma) m sb/xi^3'] = kC*(ma/mb)**2*(1+2*mb/ma)*m*sb/xi**3
    T['SO^G_a kG(3+4ma/mb) m sa/xi^3']= kG*(3+4*ma/mb)*m*sa/xi**3
    T['SO^G_b kG(4ma/mb+3(ma/mb)^2) m sb/xi^3'] = kG*(4*ma/mb+3*(ma/mb)**2)*m*sb/xi**3
    T['SS -kC(ma/mb) sa sb/xi^3']     = -kC*(ma/mb)*sa*sb/xi**3
    T['SS^G -kG(ma/mb) sa sb/xi^3']   = -kG*(ma/mb)*sa*sb/xi**3
    T['KN Q2_a -kC/(2xi^3)']          = -kC/(2*xi**3)
    T['KN Q2_b -kC(ma/mb)^2/(2xi^3)'] = -kC*(ma/mb)**2/(2*xi**3)
    T['KN M2_a -kG/(2xi^3)']          = -kG/(2*xi**3)
    T['KN M2_b -kG(ma/mb)^2/(2xi^3)'] = -kG*(ma/mb)**2/(2*xi**3)
    return T
def eps(m,xi,**kw): return sum(eps_terms(m,xi,**kw).values())
# Newtonian circular radius
kappa=kC+kG
xi1N = 4*(ma/mu)**2/kappa
print("xi1 Newtonian =",xi1N,"  2/alpha=",2/alpha)
# full circular orbit: d eps/d xi = 0 at m=1 (spins 0 to isolate orbital terms)
f=lambda x: diff(lambda y: eps(1,y,sa=0,sb=0), x)
xi1 = findroot(f, xi1N)
print("xi1 with all conservative terms (spins 0) =",xi1,"  rel shift=",(xi1-xi1N)/xi1N)
v1 = vm(1,xi1); print("v1 =",v1," alpha=",alpha," v1/alpha-1=",v1/alpha-1)
print("L/J_a = mu/ma*xi1^2*rho_1 =", (mu/ma)*xi1*vm(1,xi1))
print("cycles per orbit 2pi/rho_1 =", 2*pi*xi1/v1, "  2/alpha^2=",2/alpha**2)
T=eps_terms(1,xi1)
coul=abs(T['Coulomb -kC/xi'])
print("\nterm                                              value(per mu)        /|Coulomb|")
for k,v in T.items():
    print(f"{k:50s} {mp.nstr(v,8):>18s}  {mp.nstr(v/coul,5):>12s}")
E1=eps(1,xi1,sa=0,sb=0); print("\neps_1(xi1) =",E1,"  -alpha^2/2 =",-alpha**2/2, "  ratio=",E1/(-alpha**2/2))
# single-mode m: xi_m, eps_m
for m in (1,2,3):
    xm = findroot(lambda x: diff(lambda y: eps(m,y,sa=0,sb=0), x), m*m*xi1N)
    print(f"m={m}: xi_m={mp.nstr(xm,12)}  xi_m/(m^2 xi1)={mp.nstr(xm/(m*m*xi1),10)}  eps_m={mp.nstr(eps(m,xm,sa=0,sb=0),12)}  eps_m*m^2/eps_1={mp.nstr(eps(m,xm,sa=0,sb=0)*m*m/E1,10)}")
# phase advance per internal rad
print("\nE_m = eps_m*mu/(2ma): m=1 ->", E1*mu/(2*ma))
# node step: 62 nodes per orbit
rho1=v1/xi1; dth=2*pi/62/rho1
print("rho_1=",rho1," Delta theta_a,0 (62 nodes/orbit)=",dth," internal cycles=",dth/(2*pi))
# K_n kernel closed form at xi1 (appendix)
for n in range(0,13):
    s=sum(mpf(1)/(2*k-1) for k in range(1,n+1))
    print(f"xi*K_{n} ~ {mp.nstr((log(8*xi1)-2*s)/pi,6)}", end='; ')
print()

# ===== 第2部 =====
# Expand the relational Hamiltonian eps(xi, pi, m, s_a, s_b) into monomials  c * pi^a * m^b * s_a^c1 s_b^c2 * xi^(-p)
import sympy as sp  # noqa
xi, pi_, m, sa, sb = sp.symbols('xi pi m s_a s_b')
kC, kG, nu, g, al = sp.symbols('kappa_C kappa_G nu g alpha')   # g = m_a/mu = 1 + m_a/m_b ; ma/mb = g-1
r = g-1   # m_a/m_b
vm = 2*m*g/xi
P2 = pi_**2 + vm**2
terms = {
 'T1 rot+radial kinetic':      P2/2,
 'T2 Coulomb':                 -kC/xi,
 'T3 Newton':                  -kG/xi,
 'T4 SR kinetic':              (3*nu-1)*P2**2/8,
 'T5 Darwin':                  -nu*kC*(P2+pi_**2)/(2*xi),
 'T6 EIH velocity':            -kG*((3+nu)*P2+nu*pi_**2)/(2*xi),
 'T7 EIH static':              kG**2/(2*xi**2),
 'T8 mixed static (RN+Coulomb-gravity)': sp.Rational(3,2)*kC*kG/xi**2,
 'T9 EM 2PC velocity':         -kC/(8*xi)*(3*nu**2*pi_**4+nu*(3*nu-2)*P2**2+2*nu*(nu-1)*P2*pi_**2),
 'T10 EM 2PC kinetic':         (1-5*nu+5*nu**2)*P2**3/16,
 'T11 EM 2PC kC^2':            nu*kC**2*P2/(4*xi**2),
 'T12 EM 2PC kC^3':            nu*kC**3/(4*xi**3),
 'T13 SO_a (EM, Thomas incl.)':kC*(1+2*r)*m*sa/xi**3,
 'T14 SO_b (EM)':              kC*r**2*(1+2/r)*m*sb/xi**3,
 'T15 SO_a (grav 1PN)':        kG*(3+4*r)*m*sa/xi**3,
 'T16 SO_b (grav 1PN)':        kG*(4*r+3*r**2)*m*sb/xi**3,
 'T17 SS (EM)':                -kC*r*sa*sb/xi**3,
 'T18 SS (grav)':              -kG*r*sa*sb/xi**3,
 'T19 KN Q2 of a':             -kC/(2*xi**3),
 'T20 KN Q2 of b':             -kC*r**2/(2*xi**3),
 'T21 KN M2 of a':             -kG/(2*xi**3),
 'T22 KN M2 of b':             -kG*r**2/(2*xi**3),
}
num = {kC: sp.Float('0.0146026536534366284',20), kG: sp.Float('1.45980263966841231e-6',20),
       nu: sp.Float('5.44024290348585712e-4',20), g: 1+sp.Float('1',20)/sp.Float('1836.15267343',20), al: sp.Float('7.2973525643e-3',20)}
total = sum(terms.values())
poly = sp.Poly(sp.expand(total*xi**6), pi_, m, sa, sb, xi)   # multiply by xi^6 to make polynomial
rows=[]
for mon, coef in poly.terms():
    a,b,c1,c2,q = mon
    p = 6 - q
    rows.append((p,a,b,c1,c2,sp.simplify(coef), sp.N(coef.subs(num),10)))
rows.sort(key=lambda t:(t[0],t[1],t[2],t[3],t[4]))
print("monomial  c * pi^a * m^b * s_a^c1 * s_b^c2 * xi^(-p)")
print(f"{'p':>2} {'a':>2} {'b':>2} {'sa':>2} {'sb':>2}  {'coefficient (symbolic)':60s} numeric")
for p,a,b,c1,c2,coef,val in rows:
    print(f"{p:>2} {a:>2} {b:>2} {c1:>2} {c2:>2}  {str(coef):60s} {val}")
print("\nnumber of monomials:",len(rows))
# derivatives needed by the flow
print("\nd eps/d xi and d eps/d pi are monomials of the same table with exponents shifted (p->p+1, a->a-1).")
