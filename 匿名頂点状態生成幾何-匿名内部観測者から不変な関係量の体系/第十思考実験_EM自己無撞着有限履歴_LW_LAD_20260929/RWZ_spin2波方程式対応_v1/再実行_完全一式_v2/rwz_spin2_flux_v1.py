import argparse, json, math
import numpy as np
from scipy.special import lambertw


def r_of_rstar(rs, M):
    z = np.exp(np.clip(rs/(2*M)-1.0, -700, 700))
    return 2*M*(1.0 + np.real(lambertw(z)))


def potentials(r, M, ell):
    f = 1.0 - 2.0*M/r
    L = ell*(ell+1)
    Vrw = f*(L/r**2 - 6.0*M/r**3)
    n = (ell-1)*(ell+2)/2.0
    Vz = 2.0*f*(n*n*(n+1)*r**3 + 3*n*n*M*r**2 + 9*n*M*M*r + 9*M**3) / (r**3*(n*r + 3*M)**2)
    return Vrw, Vz


def energy(psi_t, psi_x, psi, V, dx, pref):
    return pref*np.sum(0.5*(psi_t**2 + psi_x**2 + V*psi**2))*dx


def run(parity='odd', ell=2, M=1.0, xmin=-120.0, xmax=220.0, nx=6801,
        cfl=0.45, tmax=420.0, x0=5.0, sigma=7.0, amp=1e-3,
        initial='displacement', sample_stride=50):
    x = np.linspace(xmin, xmax, nx)
    dx = x[1]-x[0]
    dt = cfl*dx
    nsteps = int(tmax/dt)
    r = r_of_rstar(x, M)
    Vrw, Vz = potentials(r, M, ell)
    V = Vrw if parity == 'odd' else Vz

    # Moncrief master field normalization factor for fluxes at horizon/infinity.
    pref = math.factorial(ell+2)/math.factorial(ell-2)/(64.0*math.pi)

    g = amp*np.exp(-0.5*((x-x0)/sigma)**2)
    if initial == 'displacement':
        psi0 = g.copy()
        pi0 = np.zeros_like(x)
    elif initial == 'outgoing':
        psi0 = g.copy()
        gx = -(x-x0)/(sigma*sigma)*g
        pi0 = -gx  # right-moving: psi_t = -psi_x
    elif initial == 'ingoing':
        psi0 = g.copy()
        gx = -(x-x0)/(sigma*sigma)*g
        pi0 = gx   # left-moving: psi_t = +psi_x
    else:
        raise ValueError(initial)

    # second-order initialization
    lap0 = np.zeros_like(x)
    lap0[1:-1] = (psi0[2:] - 2*psi0[1:-1] + psi0[:-2])/dx**2
    acc0 = lap0 - V*psi0
    psi_prev = psi0 - dt*pi0 + 0.5*dt*dt*acc0
    psi = psi0.copy()

    # Sommerfeld-compatible ghost initialization
    psi_prev[0] = psi0[0] - dt*pi0[0]
    psi_prev[-1] = psi0[-1] - dt*pi0[-1]

    flux_h = 0.0
    flux_inf = 0.0
    samples=[]

    # Initial energy computed from prescribed pi0.
    psix = np.gradient(psi0, dx)
    E0 = energy(pi0, psix, psi0, V, dx, pref)

    for nstep in range(nsteps):
        psi_next = np.empty_like(psi)
        psi_next[1:-1] = (2*psi[1:-1] - psi_prev[1:-1]
                          + (dt/dx)**2*(psi[2:] - 2*psi[1:-1] + psi[:-2])
                          - dt*dt*V[1:-1]*psi[1:-1])
        # First-order Sommerfeld at boundaries: left ingoing to horizon, right outgoing to infinity.
        c = dt/dx
        psi_next[0] = psi[0] + c*(psi[1]-psi[0])
        psi_next[-1] = psi[-1] - c*(psi[-1]-psi[-2])

        psi_t = (psi_next - psi_prev)/(2*dt)
        psi_x = np.gradient(psi, dx)
        S = -pref*psi_t*psi_x
        # outward fluxes from computational domain
        Fh = max(0.0, -S[1])
        Finf = max(0.0, S[-2])
        flux_h += Fh*dt
        flux_inf += Finf*dt

        if nstep % sample_stride == 0 or nstep == nsteps-1:
            Ein = energy(psi_t, psi_x, psi, V, dx, pref)
            samples.append({
                't': float(nstep*dt),
                'E_domain': float(Ein),
                'E_horizon_cum': float(flux_h),
                'E_infinity_cum': float(flux_inf),
                'balance_error': float((Ein+flux_h+flux_inf-E0)/E0 if E0 else 0.0),
                'psi_h': float(psi[1]),
                'psi_inf': float(psi[-2]),
            })
        psi_prev, psi = psi, psi_next

    psi_t = (psi - psi_prev)/dt
    psi_x = np.gradient(psi, dx)
    Efinal = energy(psi_t, psi_x, psi, V, dx, pref)
    result={
        'equation': 'Regge-Wheeler' if parity=='odd' else 'Zerilli-Moncrief',
        'parity': parity, 'ell': ell, 'M': M,
        'grid': {'xmin':xmin,'xmax':xmax,'nx':nx,'dx':dx,'dt':dt,'tmax':tmax,'steps':nsteps},
        'initial': {'type':initial,'x0':x0,'sigma':sigma,'amp':amp},
        'flux_prefactor': pref,
        'E0': float(E0), 'Efinal_domain': float(Efinal),
        'E_horizon': float(flux_h), 'E_infinity': float(flux_inf),
        'fraction_horizon': float(flux_h/E0), 'fraction_infinity': float(flux_inf/E0),
        'fraction_remaining': float(Efinal/E0),
        'closure': float((Efinal+flux_h+flux_inf)/E0),
        'samples': samples,
    }
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--parity', choices=['odd','even'], default='odd')
    ap.add_argument('--x0', type=float, default=5.0)
    ap.add_argument('--initial', choices=['displacement','outgoing','ingoing'], default='displacement')
    ap.add_argument('--out', required=True)
    args=ap.parse_args()
    res=run(parity=args.parity,x0=args.x0,initial=args.initial)
    with open(args.out,'w') as f: json.dump(res,f,indent=2)
    print(json.dumps({k:res[k] for k in ['equation','E0','E_horizon','E_infinity','fraction_horizon','fraction_infinity','fraction_remaining','closure']},indent=2))

if __name__=='__main__': main()
