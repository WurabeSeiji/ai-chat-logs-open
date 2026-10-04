"""Check values for the dictionary doc: derivatives of the v7 epsilon (22 monomials, v7 §1.2)
at the circular orbit m=1, pi=0, xi1, spins up/down. Pure arithmetic, no dynamics."""
import math
alpha = 7.2973525643e-3
m_a = 7.0e-4
ratio = 1836.15267343
m_b = m_a*ratio
M = m_a+m_b
mu = m_a*m_b/M
nu = mu/M
g = m_a/mu
kC = 2*alpha*g
C0m1 = 9/(m_a*m_b)
kG = kC/C0m1
k = kC+kG
# monomials: (p, a, b, sa, sb, coeff)
mono = [
 (0,2,0,0,0, 0.5),
 (0,4,0,0,0, (3*nu-1)/8),
 (0,6,0,0,0, (1-5*nu+5*nu**2)/16),
 (1,0,0,0,0, -kC-kG),
 (1,2,0,0,0, -nu*(kC+kG)-1.5*kG),
 (1,4,0,0,0, kC*nu*(1-2*nu)/2),
 (2,0,0,0,0, kG*(3*kC+kG)/2),
 (2,0,2,0,0, 2*g**2),
 (2,2,0,0,0, kC**2*nu/4),
 (2,2,2,0,0, g**2*(3*nu-1)),
 (2,4,2,0,0, 3*g**2*(1-5*nu+5*nu**2)/4),
 (3,0,0,0,0, -(kC+kG)*(g**2-2*g+2)/2 + kC**3*nu/4),
 (3,0,0,1,1, -(g-1)*(kC+kG)),
 (3,0,1,0,1, kC*(g**2-1)+kG*(3*g**2-2*g-1)),
 (3,0,1,1,0, kC*(2*g-1)+kG*(4*g-1)),
 (3,0,2,0,0, -2*g**2*(nu*(kC+kG)+3*kG)),
 (3,2,2,0,0, g**2*kC*nu*(3-4*nu)),
 (4,0,2,0,0, g**2*kC**2*nu),
 (4,0,4,0,0, g**4*(6*nu-2)),
 (4,2,4,0,0, g**4*(3-15*nu+15*nu**2)),
 (5,0,4,0,0, 2*g**4*kC*nu*(2-3*nu)),
 (6,0,6,0,0, g**6*(4-20*nu+20*nu**2)),
]
def eps(xi,pi,m,sa,sb):
    return sum(c*pi**a*m**b*sa**ca*sb**cb*xi**(-p) for p,a,b,ca,cb,c in mono)
def d_m(xi,pi,m,sa,sb):
    return sum(c*pi**a*b*m**(b-1)*sa**ca*sb**cb*xi**(-p) for p,a,b,ca,cb,c in mono if b>0)
def d_xi(xi,pi,m,sa,sb):
    return sum(-p*c*pi**a*m**b*sa**ca*sb**cb*xi**(-p-1) for p,a,b,ca,cb,c in mono if p>0)
def d_xixi(xi,pi,m,sa,sb):
    return sum(p*(p+1)*c*pi**a*m**b*sa**ca*sb**cb*xi**(-p-2) for p,a,b,ca,cb,c in mono if p>0)
def d_pipi(xi,pi,m,sa,sb):
    return sum(c*a*(a-1)*pi**(a-2)*m**b*sa**ca*sb**cb*xi**(-p) for p,a,b,ca,cb,c in mono if a>=2)

xi1 = 274.192027; m=1.0; pi=0.0; sa=1; sb=-1
print("kC,kG,nu,g", kC,kG,nu,g)
print("eps(xi1)", eps(xi1,pi,m,sa,sb))
print("d_xi (should be ~0 at circular)", d_xi(xi1,pi,m,sa,sb))
dm = d_m(xi1,pi,m,sa,sb); print("d_m", dm, " leading 4g^2 m/xi^2 =", 4*g**2*m/xi1**2)
exx = d_xixi(xi1,pi,m,sa,sb); epp = d_pipi(xi1,pi,m,sa,sb)
omX = math.sqrt(exx*epp); print("eps_xixi, eps_pipi, omega_X", exx, epp, omX, " v7: 2.661597e-5")
rho = dm/(2*g); print("rho = d_m/(2g)", rho, " gamma_cl = rho/omX", rho/omX, " -1 =", rho/omX-1)
# v7 node width
dth = 10153.0
th = dth/(4*g)*dm; thX = dth/2*omX
print("theta (v7 node)", th, " 2theta", 2*th, " deg", math.degrees(2*th), " nodes/orbit", 2*math.pi/(2*th))
print("theta_X", thX, " theta/theta_X", th/thX, " sin^2 theta", math.sin(th)**2)
# radiation: dipole Larmor circular |D''|^2 = xi^2 omega^4, omega = orbital angular freq = rho
P_D = (2/3)*kC*xi1**2*rho**4
dE_node = P_D*dth
e1 = eps(xi1,pi,m,sa,sb)
print("P_D", P_D, " dE/node", dE_node, " sin^2 theta_rad (/|eps|)", dE_node/abs(e1), " theta_rad", math.sqrt(dE_node/abs(e1)))
print("per orbit fraction", dE_node/abs(e1)*2*math.pi/(2*th))
# per-term d/dm contributions at circular orbit, grouped by v7 T-terms (analytic)
v2 = (2*m*g/xi1)**2; P2 = pi**2+v2
terms = {
 "T1": 4*g**2*m/xi1**2,
 "T4": 2*(3*nu-1)*g**2*m*P2/xi1**2,
 "T5": -4*nu*kC*g**2*m/xi1**3,
 "T6": -4*(3+nu)*kG*g**2*m/xi1**3,
 "T9": -(2*nu*kC*g**2*m/xi1**3)*((3*nu-2)*P2+(nu-1)*pi**2),
 "T10": 1.5*(1-5*nu+5*nu**2)*g**2*m*P2**2/xi1**2,
 "T11": 2*nu*kC**2*g**2*m/xi1**4,
 "T13": kC*(2*g-1)*sa/xi1**3,
 "T14": kC*((g-1)**2+2*(g-1))*sb/xi1**3,
 "T15": kG*(4*g-1)*sa/xi1**3,
 "T16": kG*(4*(g-1)+3*(g-1)**2)*sb/xi1**3,
}
tot = sum(terms.values())
print("sum of per-term d/dm", tot, " vs monomial d_m", dm, " diff", tot-dm)
for kk,v in terms.items():
    print(f"  {kk}: {v:.6e}  ratio to T1 {v/terms['T1']:.3e}")
