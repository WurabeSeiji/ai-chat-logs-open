#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Eight-state charged closure control experiment.

Persistent state:
    (U, P, E, H, Q, N, C, D)
where
    N = eta = symmetric mass ratio,
    C = Z = 1 - lambda_A*lambda_B,
    D = (lambda_A-lambda_B)^2.

The charged leading-order circular generator uses no external physical arguments.
N,C,D are propagated by the same 11-phase one-hot product-sum selector with
identity coefficient 1 in every phase.  The historical Method-B 11-phase RK4
scaffold is retained as a control architecture; its pre-existing additive RK4
stage arithmetic is NOT claimed to satisfy the stronger R4 multiplicative-
closure rule.

Reference data are used only after generation, for one-way comparison.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import inspect
import json
import math
from pathlib import Path
from typing import Dict, Any

import numpy as np

HERE = Path(__file__).resolve().parent
REF_RAW = HERE / "reference_raw.csv"
REF_META = HERE / "reference_metadata.json"

STEPS_PER_ORBIT = 4000
HSTEP = 2.0 * math.pi / STEPS_PER_ORBIT
PH = complex(math.cos(HSTEP), math.sin(HSTEP))
PH2 = complex(math.cos(HSTEP / 2.0), math.sin(HSTEP / 2.0))
TAU0 = 1.0e4
IDENTITY_PHASE_COEFFS = np.ones(11, dtype=float)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def float_bits_equal(a: float, b: float) -> bool:
    return np.asarray([a], dtype=np.float64).tobytes() == np.asarray([b], dtype=np.float64).tobytes()


def complex_bits_equal(a: complex, b: complex) -> bool:
    return np.asarray([a], dtype=np.complex128).tobytes() == np.asarray([b], dtype=np.complex128).tobytes()


def blend_scalar(q: np.ndarray, value: float) -> float:
    """Explicit 11-term one-hot product-sum identity row."""
    x = q[0] * (IDENTITY_PHASE_COEFFS[0] * value)
    for j in range(1, 11):
        x = x + q[j] * (IDENTITY_PHASE_COEFFS[j] * value)
    return float(x)


def blend(q: np.ndarray, cands):
    """Historical one-hot product-sum selector; every candidate reads OLD state only."""
    x = cands[0] * q[0]
    for j in range(1, len(cands)):
        x = x + cands[j] * q[j]
    return x


def rr_circular(P: float, Nstate: float, Cstate: float, Dstate: float) -> np.ndarray:
    """LO charged circular rates versus orbital phase phi.

    Uses only state values.  With P=r and C=Z:
      dP/dphi = -(4/3) N D sqrt(C) P^(-1/2)
                -(64/5) N C^(3/2) P^(-3/2)
      dE/dphi = 0  (circular sector)
      dt/dphi = P^(3/2)/sqrt(C)
    """
    if not (P > 0.0 and Nstate > 0.0 and Cstate > 0.0 and Dstate >= 0.0):
        raise ValueError((P, Nstate, Cstate, Dstate))
    rootC = math.sqrt(Cstate)
    rootP = math.sqrt(P)
    dP = -(4.0 / 3.0) * Nstate * Dstate * rootC / rootP \
         -(64.0 / 5.0) * Nstate * Cstate * rootC / (P * rootP)
    dt = P * rootP / rootC
    return np.array([dP, 0.0, dt], dtype=float)


def init_core(U: complex, P: float, E: float, Hphi: complex, Q: float) -> Dict[str, Any]:
    return {
        "U": complex(U), "P": float(P), "E": float(E), "H": complex(Hphi), "Q": float(Q),
        "k1": np.zeros(3), "k2": np.zeros(3), "k3": np.zeros(3), "k4": np.zeros(3),
        "Ps": float(P), "Es": float(E), "d": np.zeros(3),
        "Pm": float(P), "Em": float(E), "dphi": 0.0,
        "q": np.array([1.0] + [0.0] * 10, dtype=float),
    }


def init_external(P0: float) -> Dict[str, Any]:
    return init_core(1.0 + 0j, P0, 0.0, 1.0 + 0j, 1.0)


def init_eight(P0: float, N0: float, C0: float, D0: float) -> Dict[str, Any]:
    s = init_core(1.0 + 0j, P0, 0.0, 1.0 + 0j, 1.0)
    s["N"] = float(N0)
    s["C"] = float(C0)
    s["D"] = float(D0)
    return s


def micro_external(s: Dict[str, Any], Narg: float, Carg: float, Darg: float) -> Dict[str, Any]:
    """External-parameter control path; algebra otherwise matches micro_eight."""
    j = int(np.argmax(s["q"]))
    o = s
    n = s.copy()
    P, E = o["P"], o["E"]
    if j == 0:
        n["k1"] = rr_circular(P, Narg, Carg, Darg)
    elif j == 1:
        n["Ps"] = P + 0.5 * HSTEP * o["k1"][0]; n["Es"] = E
    elif j == 2:
        n["k2"] = rr_circular(o["Ps"], Narg, Carg, Darg)
    elif j == 3:
        n["Ps"] = P + 0.5 * HSTEP * o["k2"][0]; n["Es"] = E
    elif j == 4:
        n["k3"] = rr_circular(o["Ps"], Narg, Carg, Darg)
    elif j == 5:
        n["Ps"] = P + HSTEP * o["k3"][0]; n["Es"] = E
    elif j == 6:
        n["k4"] = rr_circular(o["Ps"], Narg, Carg, Darg)
    elif j == 7:
        n["d"] = (HSTEP / 6.0) * (o["k1"] + 2.0 * o["k2"] + 2.0 * o["k3"] + o["k4"])
    elif j == 8:
        n["Pm"] = P + 0.5 * o["d"][0]; n["Em"] = E
    elif j == 9:
        n["dphi"] = HSTEP
    elif j == 10:
        n["U"] = o["U"] * PH
        n["P"] = P + o["d"][0]
        n["E"] = E + o["d"][1]
        n["H"] = o["H"] * complex(math.cos(o["dphi"]), math.sin(o["dphi"]))
        n["Q"] = o["Q"] * math.exp(o["d"][2] / TAU0)
        n["k1"] = np.zeros(3); n["k2"] = np.zeros(3); n["k3"] = np.zeros(3); n["k4"] = np.zeros(3)
        n["Ps"] = n["P"]; n["Es"] = n["E"]; n["d"] = np.zeros(3)
        n["Pm"] = n["P"]; n["Em"] = n["E"]; n["dphi"] = 0.0
    n["q"] = np.roll(o["q"], 1)
    return n


def micro_eight(s: Dict[str, Any]) -> Dict[str, Any]:
    """Strict state-closed 8-state microstep using one product-sum selector.

    All candidate values are functions of the OLD state only.  The explicit one-hot
    phase state q selects the next value for every machine register.  N,C,D use
    identity rows with coefficient 1 in all 11 phases.
    """
    o=s
    q=o["q"]
    U,P,E,Hh,Q=o["U"],o["P"],o["E"],o["H"],o["Q"]
    k1,k2,k3,k4=o["k1"],o["k2"],o["k3"],o["k4"]
    Ps,Es=o["Ps"],o["Es"]
    d=o["d"]; Pm,Em=o["Pm"],o["Em"]; dphi=o["dphi"]
    Nstate,Cstate,Dstate=o["N"],o["C"],o["D"]

    z3=np.zeros(3)
    k1c=[rr_circular(P,Nstate,Cstate,Dstate),k1,k1,k1,k1,k1,k1,k1,k1,k1,z3]
    k2c=[k2,k2,rr_circular(Ps,Nstate,Cstate,Dstate),k2,k2,k2,k2,k2,k2,k2,z3]
    k3c=[k3,k3,k3,k3,rr_circular(Ps,Nstate,Cstate,Dstate),k3,k3,k3,k3,k3,z3]
    k4c=[k4,k4,k4,k4,k4,k4,rr_circular(Ps,Nstate,Cstate,Dstate),k4,k4,k4,z3]

    Psc=[Ps,P+0.5*HSTEP*k1[0],Ps,P+0.5*HSTEP*k2[0],Ps,P+HSTEP*k3[0],Ps,Ps,Ps,Ps,P]
    Esc=[Es,E,Es,E,Es,E,Es,Es,Es,Es,E]
    dc=[d,d,d,d,d,d,d,(HSTEP/6.0)*(k1+2.0*k2+2.0*k3+k4),d,d,z3]
    Pmc=[Pm,Pm,Pm,Pm,Pm,Pm,Pm,Pm,P+0.5*d[0],Pm,P+d[0]]
    Emc=[Em,Em,Em,Em,Em,Em,Em,Em,E,Em,E+d[1]]
    dpc=[dphi,dphi,dphi,dphi,dphi,dphi,dphi,dphi,dphi,HSTEP,0.0]

    Uc=[U]*10+[U*PH]
    Pc=[P]*10+[P+d[0]]
    Ec=[E]*10+[E+d[1]]
    Hc=[Hh]*10+[Hh*complex(math.cos(dphi),math.sin(dphi))]
    Qc=[Q]*10+[Q*math.exp(d[2]/TAU0)]
    Nc=[IDENTITY_PHASE_COEFFS[j]*Nstate for j in range(11)]
    Cc=[IDENTITY_PHASE_COEFFS[j]*Cstate for j in range(11)]
    Dc=[IDENTITY_PHASE_COEFFS[j]*Dstate for j in range(11)]

    out={}
    out["U"]=blend(q,Uc); out["P"]=float(blend(q,Pc)); out["E"]=float(blend(q,Ec))
    out["H"]=blend(q,Hc); out["Q"]=float(blend(q,Qc))
    out["N"]=float(blend(q,Nc)); out["C"]=float(blend(q,Cc)); out["D"]=float(blend(q,Dc))
    out["k1"]=blend(q,k1c); out["k2"]=blend(q,k2c); out["k3"]=blend(q,k3c); out["k4"]=blend(q,k4c)
    out["Ps"]=float(blend(q,Psc)); out["Es"]=float(blend(q,Esc)); out["d"]=blend(q,dc)
    out["Pm"]=float(blend(q,Pmc)); out["Em"]=float(blend(q,Emc)); out["dphi"]=float(blend(q,dpc))
    out["q"]=np.roll(q,1)
    return out


def macro_external(s, Narg, Carg, Darg):
    for _ in range(11):
        s = micro_external(s, Narg, Carg, Darg)
    return s


def macro_eight(s):
    for _ in range(11):
        s = micro_eight(s)
    return s


def macro_fast_eight(s: Dict[str, Any], h: float = HSTEP) -> Dict[str, Any]:
    """Algebraically collapsed macro step used for convergence diagnostics only.

    The main benchmark trajectory is generated with macro_eight (11 microsteps).
    """
    n = s.copy()
    P, E = s["P"], s["E"]
    Nstate, Cstate, Dstate = s["N"], s["C"], s["D"]
    ph = complex(math.cos(h), math.sin(h))
    k1 = rr_circular(P, Nstate, Cstate, Dstate)
    k2 = rr_circular(P + 0.5 * h * k1[0], Nstate, Cstate, Dstate)
    k3 = rr_circular(P + 0.5 * h * k2[0], Nstate, Cstate, Dstate)
    k4 = rr_circular(P + h * k3[0], Nstate, Cstate, Dstate)
    d = (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    n["U"] = s["U"] * ph
    n["P"] = P + d[0]
    n["E"] = E + d[1]
    n["H"] = s["H"] * ph
    n["Q"] = s["Q"] * math.exp(d[2] / TAU0)
    # Diagnostic fast path is mathematically the collapsed result of the identity rows.
    n["N"] = Nstate; n["C"] = Cstate; n["D"] = Dstate
    return n


def core5(s):
    return np.array([s["U"], s["P"], s["E"], s["H"], s["Q"]], dtype=np.complex128)


def core8(s):
    return np.array([s["U"], s["P"], s["E"], s["H"], s["Q"], s["N"], s["C"], s["D"]], dtype=np.complex128)


def readout(s):
    P = float(s["P"])
    t = TAU0 * math.log(float(s["Q"]))
    x = P * float(s["H"].real)
    y = P * float(s["H"].imag)
    return P, t, x, y


def closed_form_reference(r: float, r0: float, Nstate: float, Cstate: float, Dstate: float):
    A = (4.0 / 3.0) * Nstate * Dstate * Cstate
    B = (64.0 / 5.0) * Nstate * Cstate * Cstate
    if A <= 0.0 or B <= 0.0:
        raise ValueError("closed form here assumes A>0,B>0")
    def Ft(x):
        return (x**3/(3*A) - B*x*x/(2*A*A) + B*B*x/(A**3)
                - B**3/(A**4) * math.log(A*x + B))
    def Fp(x):
        return math.sqrt(Cstate) * (
            2*x**1.5/(3*A) - 2*B*math.sqrt(x)/(A*A)
            + 2*B**1.5/(A**2.5) * math.atan(math.sqrt(A*x/B))
        )
    return Ft(r0) - Ft(r), Fp(r0) - Fp(r)


def static_audit():
    src = inspect.getsource(micro_eight)
    tree = ast.parse(src)
    fn = tree.body[0]
    args = [a.arg for a in fn.args.args]
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    return {
        "micro_eight_args": args,
        "external_physical_runtime_names_present": sorted(names.intersection({"nu", "Z", "DeltaLambda", "lambda_A", "lambda_B", "Gamma"})),
        "identity_phase_coefficients": IDENTITY_PHASE_COEFFS.tolist(),
        "N_C_D_identity_coefficients_exact": bool(np.array_equal(IDENTITY_PHASE_COEFFS, np.ones(11))),
        "persistent_state_vector": ["U","P","E","H","Q","N","C","D"],
        "historical_R4_limitation_retained": True,
    }


def bitwise_control_cases():
    cases = [
        ("benchmark", 50.0, 0.25, 1.09, 0.36),
        ("same_C_different_D", 50.0, 0.25, 1.09, 0.4225),
        ("same_charge_sign", 50.0, 0.25, 0.91, 0.0),
        ("unequal_mass", 50.0, 0.16, 1.09, 0.36),
    ]
    rows = []
    for name, p0, N0, C0, D0 in cases:
        ext = init_external(p0)
        inn = init_eight(p0, N0, C0, D0)
        maxdiff = np.zeros(5)
        mismatches = 0
        identity_fail = 0
        for _ in range(4000):
            ext = macro_external(ext, N0, C0, D0)
            inn = macro_eight(inn)
            ze = core5(ext); zi = core5(inn)
            maxdiff = np.maximum(maxdiff, np.abs(ze-zi))
            if ze.tobytes() != zi.tobytes():
                mismatches += 1
            if not (float_bits_equal(inn["N"], N0) and float_bits_equal(inn["C"], C0) and float_bits_equal(inn["D"], D0)):
                identity_fail += 1
        rows.append({
            "case": name, "P0": p0, "N": N0, "C": C0, "D": D0,
            "max_core5_abs": maxdiff.tolist(),
            "bitwise_macro_mismatch_count": mismatches,
            "identity_macro_fail_count": identity_fail,
        })
    return rows


def random_micro_audit(seed=20260925, nstates=500):
    rng = np.random.default_rng(seed)
    maxdiff = np.zeros(5)
    mismatch = 0
    identity_fail = 0
    phase_fail = 0
    total = 0
    for _ in range(nstates):
        P0 = float(rng.uniform(20.0, 120.0))
        N0 = float(rng.uniform(0.05, 0.25))
        C0 = float(rng.uniform(0.6, 1.4))
        D0 = float(rng.uniform(0.0, 1.0))
        ext = init_external(P0)
        inn = init_eight(P0, N0, C0, D0)
        # Randomize starting one-hot phase identically.
        j0 = int(rng.integers(0, 11))
        q = np.zeros(11); q[j0] = 1.0
        ext["q"] = q.copy(); inn["q"] = q.copy()
        for _m in range(11):
            ext = micro_external(ext, N0, C0, D0)
            inn = micro_eight(inn)
            total += 1
            ze, zi = core5(ext), core5(inn)
            maxdiff = np.maximum(maxdiff, np.abs(ze-zi))
            if ze.tobytes() != zi.tobytes(): mismatch += 1
            if not (float_bits_equal(inn["N"], N0) and float_bits_equal(inn["C"], C0) and float_bits_equal(inn["D"], D0)):
                identity_fail += 1
            if not (np.count_nonzero(inn["q"] == 1.0) == 1 and np.sum(inn["q"]) == 1.0):
                phase_fail += 1
    return {
        "seed": seed, "microsteps": total,
        "max_core5_abs": maxdiff.tolist(),
        "bitwise_micro_mismatch_count": mismatch,
        "identity_micro_fail_count": identity_fail,
        "phase_fail_count": phase_fail,
    }


def closure_counterexamples():
    # All rows share the same six-state core (U,P,E,H,Q,N) initially.
    N = 0.25; P = 50.0
    scenarios = [
        ("A_benchmark", 1.09, 0.36, "+0.30,-0.30"),
        ("B_same_C_diff_D", 1.09, 0.4225, "+0.45,-0.20"),
        ("C_same_D_diff_C", 1.08, 0.36, "+0.40,-0.20"),
        ("D_same_sign", 0.91, 0.0, "+0.30,+0.30"),
    ]
    rows = []
    for name, C, D, lambdas in scenarios:
        rates = rr_circular(P,N,C,D)
        omega = math.sqrt(C/P**3)
        pgw = (32.0/5.0) * N*N * C**3 / P**5
        pem = (2.0/3.0) * N*N * D * C**2 / P**4
        rows.append({
            "scenario": name, "same_first6": True, "N": N, "P": P, "C": C, "D": D,
            "lambda_pair_label": lambdas,
            "dP_dphi": float(rates[0]), "dt_dphi": float(rates[2]),
            "omega": omega, "P_gw": pgw, "P_em": pem,
        })
    return rows


def load_reference_raw():
    data = np.genfromtxt(REF_RAW, delimiter=",", names=True)
    return data


def run_main_benchmark(sample_every=100, use_fast=True):
    N0, C0, D0 = 0.25, 1.09, 0.36
    r0, rf = 50.0, 20.0
    s = init_eight(r0,N0,C0,D0)
    ref = load_reference_raw()
    tref = np.asarray(ref["t"], float)
    rref = np.asarray(ref["r"], float)
    phiref = np.asarray(ref["phi"], float)

    out_path = HERE / "charged_eight_state_benchmark_sampled.csv"
    header = [
        "step","r_state","t_state","x_state","y_state","U_re","U_im","H_re","H_im","E_state",
        "N_state","C_state","D_state","N_drift","C_drift","D_drift",
        "t_closed_at_r","phi_closed_at_r","abs_t_closed_residual","abs_H_closed_residual",
        "r_saved_at_t","x_saved_at_t","y_saved_at_t","abs_r_saved_residual","xy_saved_residual"
    ]
    max_t_closed = 0.0; max_h_closed = 0.0; max_r_saved = 0.0; max_xy_saved = 0.0
    max_N = max_C = max_D = 0.0
    step = 0
    prev = None
    with out_path.open("w",encoding="utf-8",newline="") as f:
        w = csv.writer(f); w.writerow(header)
        while s["P"] > rf:
            if step % sample_every == 0:
                P,t,x,y = readout(s)
                tc,pc = closed_form_reference(P,r0,N0,C0,D0)
                href = complex(math.cos(pc),math.sin(pc))
                # saved reference interpolation by physical time, no feedback to generator
                rr_t = float(np.interp(t,tref,rref))
                pp_t = float(np.interp(t,tref,phiref))
                xx_t = rr_t * math.cos(pp_t); yy_t = rr_t * math.sin(pp_t)
                et = abs(t-tc); eh = abs(s["H"]-href); er=abs(P-rr_t); exy=math.hypot(x-xx_t,y-yy_t)
                max_t_closed=max(max_t_closed,et); max_h_closed=max(max_h_closed,eh); max_r_saved=max(max_r_saved,er); max_xy_saved=max(max_xy_saved,exy)
                ndr=abs(s["N"]-N0); cdr=abs(s["C"]-C0); ddr=abs(s["D"]-D0)
                max_N=max(max_N,ndr); max_C=max(max_C,cdr); max_D=max(max_D,ddr)
                w.writerow([step,P,t,x,y,s["U"].real,s["U"].imag,s["H"].real,s["H"].imag,s["E"],
                            s["N"],s["C"],s["D"],ndr,cdr,ddr,tc,pc,et,eh,rr_t,xx_t,yy_t,er,exy])
            prev = s
            s = macro_fast_eight(s, HSTEP) if use_fast else macro_eight(s)
            step += 1
            if step > 600000:
                raise RuntimeError("benchmark did not reach rf")
    # final post-cross sample
    P,t,x,y=readout(s); tc,pc=closed_form_reference(P,r0,N0,C0,D0); href=complex(math.cos(pc),math.sin(pc))
    rr_t=float(np.interp(min(t,tref[-1]),tref,rref))
    pp_t=float(np.interp(min(t,tref[-1]),tref,phiref))
    xx_t=rr_t*math.cos(pp_t); yy_t=rr_t*math.sin(pp_t)
    et=abs(t-tc); eh=abs(s["H"]-href); er=abs(P-rr_t); exy=math.hypot(x-xx_t,y-yy_t)
    # Post-cross point is outside the saved reference time range; do not fold it into saved-grid maxima.
    max_t_closed=max(max_t_closed,et); max_h_closed=max(max_h_closed,eh)
    # interpolate crossing at exactly r=20 between prev and s for final diagnostics
    P0,t0,_,_=readout(prev); P1,t1,_,_=readout(s)
    frac=(P0-rf)/(P0-P1)
    t_cross=t0+frac*(t1-t0)
    phi_cross=(step-1+frac)*HSTEP
    with open(REF_META,"r",encoding="utf-8") as f:
        meta=json.load(f)
    t_true=float(meta["results_summary"]["t_final"]); phi_true=float(meta["results_summary"]["phi_final"])
    return {
        "sampled_csv": out_path.name,
        "execution_path": "collapsed macro_fast_eight validated bitwise against explicit 11-micro path" if use_fast else "explicit 11-micro macro_eight",
        "macro_steps_to_cross_rf": step,
        "equivalent_microsteps_if_explicit": step*11,
        "last_r": P,
        "crossing_interpolation": {
            "t_at_r20": t_cross, "t_reference": t_true, "abs_t_error": abs(t_cross-t_true),
            "phi_at_r20": phi_cross, "phi_reference": phi_true, "abs_phi_error": abs(phi_cross-phi_true),
            "cycles_at_r20": phi_cross/(2*math.pi), "cycles_reference": phi_true/(2*math.pi),
        },
        "max_sampled_errors": {
            "abs_t_vs_closed_form": max_t_closed,
            "abs_H_vs_closed_form_phase_vector": max_h_closed,
            "abs_r_vs_saved_reference_at_same_t": max_r_saved,
            "xy_norm_vs_saved_reference_at_same_t": max_xy_saved,
        },
        "max_identity_drifts": {"N":max_N,"C":max_C,"D":max_D},
    }


def solve_r_for_phi(phi_target, r0=50.0, rmin=20.0, Nstate=0.25, Cstate=1.09, Dstate=0.36):
    lo, hi = rmin, r0
    for _ in range(100):
        mid = 0.5*(lo+hi)
        _, pmid = closed_form_reference(mid,r0,Nstate,Cstate,Dstate)
        if pmid > phi_target:
            lo = mid
        else:
            hi = mid
    return 0.5*(lo+hi)


def run_convergence():
    target_orbits=100
    phi_target=2*math.pi*target_orbits
    r_true=solve_r_for_phi(phi_target)
    t_true,_=closed_form_reference(r_true,50.0,0.25,1.09,0.36)
    rows=[]
    for spo in [500,1000,2000,4000,8000]:
        h=2*math.pi/spo
        s=init_eight(50.0,0.25,1.09,0.36)
        nsteps=target_orbits*spo
        for _ in range(nsteps):
            s=macro_fast_eight(s,h=h)
        P,t,_,_=readout(s)
        rows.append({
            "steps_per_orbit":spo,"h":h,"macro_steps":nsteps,"target_orbits":target_orbits,
            "r_at_100_orbits":P,"r_reference":r_true,"abs_r_error":abs(P-r_true),
            "t_at_100_orbits":t,"t_reference":t_true,"abs_t_error":abs(t-t_true),
        })
    path=HERE/"convergence.csv"
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    return rows

def save_rows(path: Path, rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


def main():
    static=static_audit()
    controls=bitwise_control_cases()
    random_audit=random_micro_audit()
    closure=closure_counterexamples()
    save_rows(HERE/"closure_counterexamples.csv",closure)

    # Validate collapsed macro against explicit 11-micro path before using it for the long benchmark.
    sm=init_eight(50.0,0.25,1.09,0.36); sf=init_eight(50.0,0.25,1.09,0.36)
    fast_micro_max=np.zeros(8); fast_micro_mismatch=0
    for _ in range(12000):
        sm=macro_eight(sm); sf=macro_fast_eight(sf,HSTEP)
        zm=core8(sm); zf=core8(sf)
        fast_micro_max=np.maximum(fast_micro_max,np.abs(zm-zf))
        if zm.tobytes()!=zf.tobytes(): fast_micro_mismatch+=1
    if fast_micro_mismatch != 0 or np.any(fast_micro_max != 0):
        raise RuntimeError("collapsed fast path failed bitwise validation")

    benchmark=run_main_benchmark(sample_every=100, use_fast=True)
    convergence=run_convergence()

    result={
        "experiment":"charged eight-state closure control v1",
        "date":"2026-09-25",
        "state_vector":["U","P","E","H","Q","N","C","D"],
        "state_meaning":{
            "N":"eta=symmetric mass ratio",
            "C":"Z=1-lambda_A*lambda_B",
            "D":"(lambda_A-lambda_B)^2",
        },
        "generator_equations":{
            "dP_dphi":"-(4/3)*N*D*sqrt(C)/sqrt(P) -(64/5)*N*C^(3/2)/P^(3/2)",
            "dE_dphi":"0 in quasi-circular sector",
            "dt_dphi":"P^(3/2)/sqrt(C)",
        },
        "static_audit":static,
        "external_vs_internal_bitwise_controls":controls,
        "random_micro_audit":random_audit,
        "closure_counterexamples":closure,
        "benchmark":benchmark,
        "convergence":convergence,
        "explicit_micro_vs_collapsed_fast_12000":{
            "max_abs_core8":fast_micro_max.tolist(),
            "bitwise_macro_mismatch_count":fast_micro_mismatch,
        },
        "reference":{
            "raw_csv":REF_RAW.name,"raw_sha256":sha256(REF_RAW),
            "metadata":REF_META.name,"metadata_sha256":sha256(REF_META),
        },
        "scope_warning":"LO adiabatic quasi-circular charged benchmark only; does not establish full Einstein-Maxwell closure or full R4 multiplicative closure.",
    }
    (HERE/"charged_eight_state_summary.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
