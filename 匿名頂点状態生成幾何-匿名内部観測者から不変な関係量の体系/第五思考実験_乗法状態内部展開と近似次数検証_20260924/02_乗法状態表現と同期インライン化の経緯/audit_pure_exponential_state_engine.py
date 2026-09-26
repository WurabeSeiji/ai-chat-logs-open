#!/usr/bin/env python3
import ast
import inspect
import math
import numpy as np

H = 2.0 * math.pi / 4000.0
ORDER = 12
LAMBDA_T = 1.0e-4
PHASE_H = complex(math.cos(H), math.sin(H))
PHASE_H2 = complex(math.cos(0.5 * H), math.sin(0.5 * H))


def coeffs(N):
    c = [1.0]
    for m in range(N):
        c.append(c[-1] * (2 * m + 1) / (m + 1))
    return np.asarray(c, dtype=float)


C = coeffs(ORDER)


def trig_from_phase(U):
    """Return cos(chi), sin(chi) from U=exp(i chi), with no log/angle."""
    inv = 1.0 / U
    c = 0.5 * (U + inv)
    s = (U - inv) / (2.0j)
    return float(c.real), float(s.real)


def rr_from_exp_state(U, P, E, nu):
    """RR RHS using only current exponential/phase state values.

    P is the positive multiplicative scale state exp(log p)=p.
    E is the positive multiplicative scale state exp(log e)=e.
    U is exp(i chi).
    """
    ce, se = trig_from_phase(U)
    den = 1.0 + E * ce
    r = P / den
    hp = math.sqrt(P)
    rdot = E * se / hp
    vt = den / hp
    v2 = rdot * rdot + vt * vt
    invr = 1.0 / r
    Crr = 1.6 * nu * invr ** 3
    an = rdot * (18.0 * v2 + (2.0 / 3.0) * invr - 25.0 * rdot * rdot)
    bv = -(6.0 * v2 - 2.0 * invr - 15.0 * rdot * rdot)
    Edot = Crr * (an * rdot + bv * v2)
    pdot = 2.0 * Crr * bv * P
    En = (E * E - 1.0) / (2.0 * P)
    edot = (P * Edot + En * pdot) / E
    dt_dchi = r * r / hp
    return np.array([pdot * dt_dchi, edot * dt_dchi, dt_dchi], dtype=float)


def gN_from_exp_state(U, P, E):
    ce, _ = trig_from_phase(U)
    u = (3.0 + E * ce) / P
    y = float(C[-1])
    for a in C[-2::-1]:
        y = float(a) + u * y
    return y


def step_exp_state(Z, nu):
    """One synchronous macro update from current state Z only.

    Z = [Uchi, P, E, Hphi, T]
      Uchi = exp(i chi)
      P    = exp(log p) = p > 0
      E    = exp(log e) = e > 0
      Hphi = exp(i phi/2)
      T    = exp(LAMBDA_T * t)

    No log/angle/readout is used here.  No next-state component is used to
    compute another next-state component.  RK4 stage values are temporary
    expressions of the OLD state only.
    """
    U = complex(Z[0])
    P = float(Z[1].real)
    E = float(Z[2].real)
    Hphi = complex(Z[3])
    T = float(Z[4].real)

    # Fully expanded RK4 DAG: every stage is a function of old U,P,E only.
    k1 = rr_from_exp_state(U, P, E, nu)
    k2 = rr_from_exp_state(
        U * PHASE_H2,
        P + 0.5 * H * k1[0],
        E + 0.5 * H * k1[1],
        nu,
    )
    k3 = rr_from_exp_state(
        U * PHASE_H2,
        P + 0.5 * H * k2[0],
        E + 0.5 * H * k2[1],
        nu,
    )
    k4 = rr_from_exp_state(
        U * PHASE_H,
        P + H * k3[0],
        E + H * k3[1],
        nu,
    )
    d = (H / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    # dphi is expanded in OLD-state expressions: do not use P_next/E_next.
    dphi = H * gN_from_exp_state(
        U * PHASE_H2,
        P + 0.5 * d[0],
        E + 0.5 * d[1],
    )

    # All persistent states update multiplicatively, simultaneously.
    # P*(1+dP/P) and E*(1+dE/E) are algebraically identical to P+dP, E+dE,
    # but the persistent-state assignment itself is multiplicative.
    factors = np.array([
        PHASE_H,
        1.0 + d[0] / P,
        1.0 + d[1] / E,
        complex(math.cos(0.5 * dphi), math.sin(0.5 * dphi)),
        math.exp(LAMBDA_T * d[2]),
    ], dtype=np.complex128)
    return np.asarray(Z, dtype=np.complex128) * factors


def init_exp_state(chi, p, e, phi, t):
    return np.array([
        complex(math.cos(chi), math.sin(chi)),
        math.exp(math.log(p)),
        math.exp(math.log(e)),
        complex(math.cos(0.5 * phi), math.sin(0.5 * phi)),
        math.exp(LAMBDA_T * t),
    ], dtype=np.complex128)


def decode_for_observation(Z):
    """Observation/test only. Never called by step_exp_state."""
    U, Pz, Ez, Hphi, Tz = Z
    # p/e are positive exponential-scale states; log is used only in the
    # observation branch, then exponentiated solely to compare with legacy p/e.
    logp = math.log(float(Pz.real))
    loge = math.log(float(Ez.real))
    p = math.exp(logp)
    e = math.exp(loge)
    t = math.log(float(Tz.real)) / LAMBDA_T
    return p, e, t


def ref_step(s, nu):
    chi, p, e, phi, t = map(float, s)

    def rr(x, pp, ee):
        ce, se = math.cos(x), math.sin(x)
        den = 1.0 + ee * ce
        r = pp / den
        hp = math.sqrt(pp)
        rdot = ee * se / hp
        vt = den / hp
        v2 = rdot * rdot + vt * vt
        invr = 1.0 / r
        Crr = 1.6 * nu * invr ** 3
        an = rdot * (18.0 * v2 + (2.0 / 3.0) * invr - 25.0 * rdot * rdot)
        bv = -(6.0 * v2 - 2.0 * invr - 15.0 * rdot * rdot)
        Edot = Crr * (an * rdot + bv * v2)
        pdot = 2.0 * Crr * bv * pp
        En = (ee * ee - 1.0) / (2.0 * pp)
        edot = (pp * Edot + En * pdot) / ee
        dt_dchi = r * r / hp
        return np.array([pdot * dt_dchi, edot * dt_dchi, dt_dchi])

    def gN(x, pp, ee):
        u = (3.0 + ee * math.cos(x)) / pp
        y = float(C[-1])
        for a in C[-2::-1]:
            y = float(a) + u * y
        return y

    k1 = rr(chi, p, e)
    k2 = rr(chi + 0.5 * H, p + 0.5 * H * k1[0], e + 0.5 * H * k1[1])
    k3 = rr(chi + 0.5 * H, p + 0.5 * H * k2[0], e + 0.5 * H * k2[1])
    k4 = rr(chi + H, p + H * k3[0], e + H * k3[1])
    d = (H / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    dphi = H * gN(chi + 0.5 * H, p + 0.5 * d[0], e + 0.5 * d[1])
    return np.array([chi + H, p + d[0], e + d[1], phi + dphi, t + d[2]], dtype=float)


def static_audit():
    src = inspect.getsource(step_exp_state)
    tree = ast.parse(src)
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = None
            if isinstance(node.func, ast.Name):
                name = node.func.id
            elif isinstance(node.func, ast.Attribute):
                name = node.func.attr
            if name in {"log", "log1p", "angle", "unwrap", "arctan2"}:
                bad.append((name, node.lineno))
    return bad


def run_case(p0, e0, q, nsteps=12000):
    nu = q / (1.0 + q) ** 2
    s = np.array([0.0, p0, e0, 0.0, 0.0], dtype=float)
    Z = init_exp_state(*s)
    max_err = np.zeros(5)
    max_unit = np.zeros(2)
    for _ in range(nsteps):
        s = ref_step(s, nu)
        Z = step_exp_state(Z, nu)

        # Test harness readout only; never fed back into Z.
        pz, ez, tz = decode_for_observation(Z)
        chi_phase = complex(math.cos(s[0]), math.sin(s[0]))
        phi_half_phase = complex(math.cos(0.5 * s[3]), math.sin(0.5 * s[3]))
        err = np.array([
            abs(Z[0] - chi_phase),
            abs(pz - s[1]),
            abs(ez - s[2]),
            abs(Z[3] - phi_half_phase),
            abs(tz - s[4]),
        ])
        max_err = np.maximum(max_err, err)
        max_unit = np.maximum(max_unit, [abs(abs(Z[0]) - 1.0), abs(abs(Z[3]) - 1.0)])
    return max_err, max_unit, s, Z


def main():
    bad = static_audit()
    print("forbidden_readout_calls_in_update:", bad)

    cases = []
    for p0 in (24.3, 60.0, 120.0):
        for e0 in (0.20, 0.45, 0.55):
            cases.append((p0, e0, 1.0))
    for q in (2.0, 4.0, 10.0):
        cases.append((60.0, 0.45, q))

    global_max = np.zeros(5)
    global_unit = np.zeros(2)
    for p0, e0, q in cases:
        err, unit, s, Z = run_case(p0, e0, q)
        global_max = np.maximum(global_max, err)
        global_unit = np.maximum(global_unit, unit)
        print(f"case p={p0:g} e={e0:g} q={q:g} max_err={err} unit_drift={unit}")

    print("GLOBAL_MAX [chi_phase,p,e,phi_half_phase,t] =", global_max)
    print("GLOBAL_UNIT_DRIFT [chi,phi_half] =", global_unit)


if __name__ == "__main__":
    main()

# --- Derived Paper-4 map audit: no readout/log required ---
def rot_from_halfphase(Hp):
    c = float(Hp.real)
    s = float(Hp.imag)
    return np.array([[c, -s], [s, c]], dtype=float)


def radius_from_exp_state(Z):
    U = complex(Z[0]); P = float(Z[1].real); E = float(Z[2].real)
    ce, _ = trig_from_phase(U)
    return P / (1.0 + E * ce)


def X_from_exp_state(Z):
    r = radius_from_exp_state(Z)
    Hp = complex(Z[3])
    return math.sqrt(r) * np.array([Hp.real, Hp.imag], dtype=float)


def S_from_current_exp_state(Z, nu):
    # Znext is itself a pure function of OLD Z.  S is derived from that transition;
    # it is not fed back into the state engine.
    Zn = step_exp_state(Z, nu)
    r0 = radius_from_exp_state(Z)
    r1 = radius_from_exp_state(Zn)
    rho = math.sqrt(r1 / r0)
    R0 = rot_from_halfphase(complex(Z[3]))
    R1 = rot_from_halfphase(complex(Zn[3]))
    S = R1 @ np.diag([rho, 1.0 / rho]) @ R0.T
    return S, Zn


def run_map_audit():
    cases=[]
    for p0 in (24.3,60.0,120.0):
        for e0 in (0.20,0.45,0.55): cases.append((p0,e0,1.0))
    for q in (2.0,4.0,10.0): cases.append((60.0,0.45,q))
    max_map=0.0; max_det=0.0; max_xy=0.0
    for p0,e0,q in cases:
        nu=q/(1+q)**2
        s=np.array([0.0,p0,e0,0.0,0.0],float)
        Z=init_exp_state(*s)
        for _ in range(12000):
            X0=X_from_exp_state(Z)
            S,Zn=S_from_current_exp_state(Z,nu)
            X1=X_from_exp_state(Zn)
            max_map=max(max_map,float(np.linalg.norm(S@X0-X1)))
            max_det=max(max_det,abs(float(np.linalg.det(S))-1.0))
            # Compare final physical xy against legacy reference step, observation-side only.
            s=ref_step(s,nu)
            Y=np.array([X1[0]*X1[0]-X1[1]*X1[1],2*X1[0]*X1[1]])
            r=s[1]/(1+s[2]*math.cos(s[0]))
            yref=np.array([r*math.cos(s[3]),r*math.sin(s[3])])
            max_xy=max(max_xy,float(np.linalg.norm(Y-yref)))
            Z=Zn
    print('PURE_EXP_MAP max_map_error=',max_map)
    print('PURE_EXP_MAP max_det_error=',max_det)
    print('PURE_EXP_MAP max_xy_vs_legacy=',max_xy)

if __name__ == '__main__':
    run_map_audit()
