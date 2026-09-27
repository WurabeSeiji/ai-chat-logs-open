#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 8 strict re-run / 8A2

Goal
----
Execute the minimal aa-ab-bb coupled experiment strictly as a single synchronous
state-dependent array action

    Psi_{k+1} = S(Psi_k) @ Psi_k

with no separate radiation-amplitude layer, no post-hoc normalization, and no
sequential local-update-then-coupling shortcut.

State
-----
Psi = [Psi_aa, Psi_ab, Psi_bb], each channel using the exact 41-real full
microstate layout from Paper 7. Total dimension = 123.

Local blocks
------------
For each 41-state channel z, construct an explicit 41x41 state-dependent matrix
S_local(z) that exactly reproduces the accepted Paper-7 transition(z):

    S_local(z) @ z == transition(z)

The first 30 rows are the existing 11-phase one-hot product-sum written as matrix
rows on the q columns; the q rows are the same cyclic shift.

Cross-channel hypothesis (minimal, not derived)
-----------------------------------------------
- aa := a^2, ab := a*b, bb := b^2 only for the initial phase assignment.
- common orbit normal (coplanar, same orientation).
- initial relation phases: (phi_aa, phi_ab, phi_bb) = (delta, 0, -delta).
- scalar cross coefficient k_ij = cos(phi_i - phi_j), evaluated from the CURRENT
  H phase in Psi at every microstep.
- the same scalar coefficient multiplies the entire source local-action block.

Thus the full 123x123 action is built directly as

    S_full[block i, block j] = k_ij(Psi) * S_local(Psi_j)

and one synchronous step is exactly

    Psi' = S_full(Psi) @ Psi.

No row normalization, no clipping, no renormalization, no reference trajectory,
no separate GW/EM amplitude cancellation calculation.

Readout-only diagnostics
------------------------
P, C, D, H phase, q sums, finite/domain checks are recorded but never fed back.
The Paper-6/7 interpretation that radiative dissipation appears as orbital P
shrink is retained; no new radiation physics is added.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(parents=True, exist_ok=True)

STEPS_PER_ORBIT = 4000
HSTEP = 2.0 * math.pi / STEPS_PER_ORBIT
PH_RE = math.cos(HSTEP)
PH_IM = math.sin(HSTEP)
TAU0 = 1.0e4
P0 = 50.0
N0 = 0.25

UR,UI,P,E,HR,HI,Q,N,C,D = range(10)
K1=10; K2=13; K3=16; K4=19
PS=22; ES=23; DV=24; PM=27; EM=28; DPHI=29; QB=30
NST=41
NCH=3
NTOT=NST*NCH
CHANNELS=("aa","ab","bb")

# n=1 opposite-sign ab source, with self terms from Paper 7.
CHANNEL_CD=((0.91,0.0),(1.09,0.36),(0.91,0.0))


def rates(p: float, n: float, c: float, d: float):
    if not (math.isfinite(p) and math.isfinite(c) and p > 0.0 and c > 0.0):
        raise ValueError(f"rates domain failure: p={p}, c={c}")
    rp=math.sqrt(p); rc=math.sqrt(c)
    dp=-(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
    dt=p*rp/rc
    return dp,0.0,dt


def transition_reference(z: np.ndarray) -> np.ndarray:
    """Exact direct Paper-7 transition, used only for audit/control comparison."""
    q=z[QB:QB+11]
    p=z[P]; e=z[E]; n=z[N]; c=z[C]; dd=z[D]
    k1=z[K1:K1+3]; k2=z[K2:K2+3]; k3=z[K3:K3+3]; k4=z[K4:K4+3]
    ps=z[PS]; es=z[ES]; dv=z[DV:DV+3]; pm=z[PM]; em=z[EM]; dphi=z[DPHI]
    rP,rE,rT=rates(p,n,c,dd)
    rPsP,rPsE,rPsT=rates(ps,n,c,dd)
    cand=np.empty((30,11),dtype=np.float64)
    for i in range(30):
        for j in range(11): cand[i,j]=z[i]
    cand[K1,0]=rP; cand[K1+1,0]=rE; cand[K1+2,0]=rT
    for a in range(3): cand[K1+a,10]=0.0
    cand[PS,1]=p+0.5*HSTEP*k1[0]; cand[ES,1]=e+0.5*HSTEP*k1[1]
    cand[K2,2]=rPsP; cand[K2+1,2]=rPsE; cand[K2+2,2]=rPsT
    for a in range(3): cand[K2+a,10]=0.0
    cand[PS,3]=p+0.5*HSTEP*k2[0]; cand[ES,3]=e+0.5*HSTEP*k2[1]
    cand[K3,4]=rPsP; cand[K3+1,4]=rPsE; cand[K3+2,4]=rPsT
    for a in range(3): cand[K3+a,10]=0.0
    cand[PS,5]=p+HSTEP*k3[0]; cand[ES,5]=e+HSTEP*k3[1]
    cand[K4,6]=rPsP; cand[K4+1,6]=rPsE; cand[K4+2,6]=rPsT
    for a in range(3): cand[K4+a,10]=0.0
    for a in range(3):
        cand[DV+a,7]=(HSTEP/6.0)*(k1[a]+2.0*k2[a]+2.0*k3[a]+k4[a])
        cand[DV+a,10]=0.0
    cand[PM,8]=p+0.5*dv[0]; cand[EM,8]=e+0.5*dv[1]
    cand[PM,10]=p+dv[0]; cand[EM,10]=e+dv[1]
    cand[DPHI,9]=HSTEP; cand[DPHI,10]=0.0
    ur=z[UR]; ui=z[UI]; hr=z[HR]; hi=z[HI]
    cand[UR,10]=ur*PH_RE-ui*PH_IM; cand[UI,10]=ur*PH_IM+ui*PH_RE
    cand[P,10]=p+dv[0]; cand[E,10]=e+dv[1]
    cd=math.cos(dphi); sd=math.sin(dphi)
    cand[HR,10]=hr*cd-hi*sd; cand[HI,10]=hr*sd+hi*cd
    cand[Q,10]=z[Q]*math.exp(dv[2]/TAU0)
    cand[PS,10]=p+dv[0]; cand[ES,10]=e+dv[1]
    out=np.empty(NST,dtype=np.float64)
    for i in range(30): out[i]=float(np.dot(q,cand[i]))
    out[QB]=q[10]
    for j in range(1,11): out[QB+j]=q[j-1]
    return out


def local_action_matrix(z: np.ndarray) -> np.ndarray:
    """Return explicit 41x41 S_local(z) with S_local(z)@z == transition_reference(z)."""
    q=z[QB:QB+11]
    p=z[P]; e=z[E]; n=z[N]; c=z[C]; dd=z[D]
    k1=z[K1:K1+3]; k2=z[K2:K2+3]; k3=z[K3:K3+3]; k4=z[K4:K4+3]
    ps=z[PS]; es=z[ES]; dv=z[DV:DV+3]; pm=z[PM]; em=z[EM]; dphi=z[DPHI]
    rP,rE,rT=rates(p,n,c,dd)
    rPsP,rPsE,rPsT=rates(ps,n,c,dd)
    cand=np.empty((30,11),dtype=np.float64)
    for i in range(30):
        for j in range(11): cand[i,j]=z[i]
    cand[K1,0]=rP; cand[K1+1,0]=rE; cand[K1+2,0]=rT
    for a in range(3): cand[K1+a,10]=0.0
    cand[PS,1]=p+0.5*HSTEP*k1[0]; cand[ES,1]=e+0.5*HSTEP*k1[1]
    cand[K2,2]=rPsP; cand[K2+1,2]=rPsE; cand[K2+2,2]=rPsT
    for a in range(3): cand[K2+a,10]=0.0
    cand[PS,3]=p+0.5*HSTEP*k2[0]; cand[ES,3]=e+0.5*HSTEP*k2[1]
    cand[K3,4]=rPsP; cand[K3+1,4]=rPsE; cand[K3+2,4]=rPsT
    for a in range(3): cand[K3+a,10]=0.0
    cand[PS,5]=p+HSTEP*k3[0]; cand[ES,5]=e+HSTEP*k3[1]
    cand[K4,6]=rPsP; cand[K4+1,6]=rPsE; cand[K4+2,6]=rPsT
    for a in range(3): cand[K4+a,10]=0.0
    for a in range(3):
        cand[DV+a,7]=(HSTEP/6.0)*(k1[a]+2.0*k2[a]+2.0*k3[a]+k4[a])
        cand[DV+a,10]=0.0
    cand[PM,8]=p+0.5*dv[0]; cand[EM,8]=e+0.5*dv[1]
    cand[PM,10]=p+dv[0]; cand[EM,10]=e+dv[1]
    cand[DPHI,9]=HSTEP; cand[DPHI,10]=0.0
    ur=z[UR]; ui=z[UI]; hr=z[HR]; hi=z[HI]
    cand[UR,10]=ur*PH_RE-ui*PH_IM; cand[UI,10]=ur*PH_IM+ui*PH_RE
    cand[P,10]=p+dv[0]; cand[E,10]=e+dv[1]
    cd=math.cos(dphi); sd=math.sin(dphi)
    cand[HR,10]=hr*cd-hi*sd; cand[HI,10]=hr*sd+hi*cd
    cand[Q,10]=z[Q]*math.exp(dv[2]/TAU0)
    cand[PS,10]=p+dv[0]; cand[ES,10]=e+dv[1]

    S=np.zeros((NST,NST),dtype=np.float64)
    # First 30 outputs: exact one-hot product-sum, written as matrix entries on q columns.
    for i in range(30):
        for j in range(11):
            S[i,QB+j]=cand[i,j]
    # q cyclic shift, exactly as Paper 7.
    S[QB, QB+10]=1.0
    for j in range(1,11):
        S[QB+j, QB+j-1]=1.0
    return S


def phase_from_H(z: np.ndarray) -> float:
    hr=float(z[HR]); hi=float(z[HI])
    if not (math.isfinite(hr) and math.isfinite(hi)) or (hr==0.0 and hi==0.0):
        raise ValueError(f"invalid H phase vector: ({hr},{hi})")
    return math.atan2(hi,hr)


def K_from_current_state(psi: np.ndarray) -> np.ndarray:
    blocks=psi.reshape(NCH,NST)
    ph=np.array([phase_from_H(blocks[i]) for i in range(NCH)],dtype=np.float64)
    K=np.empty((NCH,NCH),dtype=np.float64)
    for i in range(NCH):
        for j in range(NCH):
            K[i,j]=1.0 if i==j else math.cos(ph[i]-ph[j])
    return K


def full_action_matrix(psi: np.ndarray, cross: bool=True) -> tuple[np.ndarray,np.ndarray]:
    blocks=psi.reshape(NCH,NST)
    Sloc=[local_action_matrix(blocks[j]) for j in range(NCH)]
    K=K_from_current_state(psi) if cross else np.eye(NCH,dtype=np.float64)
    S=np.zeros((NTOT,NTOT),dtype=np.float64)
    for i in range(NCH):
        r=slice(i*NST,(i+1)*NST)
        for j in range(NCH):
            c=slice(j*NST,(j+1)*NST)
            S[r,c]=K[i,j]*Sloc[j]
    return S,K


def strict_step(psi: np.ndarray, cross: bool=True) -> tuple[np.ndarray,np.ndarray]:
    S,K=full_action_matrix(psi,cross=cross)
    return S @ psi, K


def init_channel(c0: float,d0: float,phi0: float) -> np.ndarray:
    z=np.zeros(NST,dtype=np.float64)
    z[UR]=1.0; z[P]=P0
    z[HR]=math.cos(phi0); z[HI]=math.sin(phi0)
    z[Q]=1.0; z[N]=N0; z[C]=c0; z[D]=d0
    z[PS]=P0; z[PM]=P0; z[QB]=1.0
    return z


def init_psi(delta: float) -> np.ndarray:
    phases=(delta,0.0,-delta)
    return np.concatenate([
        init_channel(CHANNEL_CD[i][0],CHANNEL_CD[i][1],phases[i]) for i in range(NCH)
    ])


def validity(psi: np.ndarray) -> tuple[bool,str]:
    if not np.all(np.isfinite(psi)):
        return False,"nonfinite_state"
    b=psi.reshape(NCH,NST)
    for i,name in enumerate(CHANNELS):
        if b[i,P] <= 0.0: return False,f"{name}:P<=0"
        if b[i,C] <= 0.0: return False,f"{name}:C<=0"
        if b[i,PS] <= 0.0: return False,f"{name}:PS<=0"
        if b[i,HR]==0.0 and b[i,HI]==0.0: return False,f"{name}:H_zero"
    return True,"ok"


def run_microsteps(psi0: np.ndarray, n_micro: int, cross: bool=True):
    psi=psi0.copy()
    rows=[]
    for m in range(n_micro+1):
        b=psi.reshape(NCH,NST)
        try:
            K=K_from_current_state(psi) if cross else np.eye(NCH)
        except Exception:
            K=np.full((NCH,NCH),np.nan)
        rows.append({
            "micro":m,
            "P_aa":float(b[0,P]),"P_ab":float(b[1,P]),"P_bb":float(b[2,P]),
            "C_aa":float(b[0,C]),"C_ab":float(b[1,C]),"C_bb":float(b[2,C]),
            "D_aa":float(b[0,D]),"D_ab":float(b[1,D]),"D_bb":float(b[2,D]),
            "phase_aa":float(math.atan2(b[0,HI],b[0,HR])) if (b[0,HR]!=0 or b[0,HI]!=0) else math.nan,
            "phase_ab":float(math.atan2(b[1,HI],b[1,HR])) if (b[1,HR]!=0 or b[1,HI]!=0) else math.nan,
            "phase_bb":float(math.atan2(b[2,HI],b[2,HR])) if (b[2,HR]!=0 or b[2,HI]!=0) else math.nan,
            "qsum_aa":float(np.sum(b[0,QB:QB+11])),"qsum_ab":float(np.sum(b[1,QB:QB+11])),"qsum_bb":float(np.sum(b[2,QB:QB+11])),
            "K00":float(K[0,0]),"K01":float(K[0,1]),"K02":float(K[0,2]),
            "K10":float(K[1,0]),"K11":float(K[1,1]),"K12":float(K[1,2]),
            "K20":float(K[2,0]),"K21":float(K[2,1]),"K22":float(K[2,2]),
        })
        if m==n_micro: break
        ok,why=validity(psi)
        if not ok: return psi,rows,False,why,m
        try:
            psi,_=strict_step(psi,cross=cross)
        except Exception as e:
            return psi,rows,False,f"exception:{type(e).__name__}:{e}",m
    ok,why=validity(psi)
    return psi,rows,ok,why,n_micro


def audit_local_matrix() -> dict:
    # Audit on all three initial channels and on all 11 local microsteps without cross coupling.
    delta=0.731
    psi=init_psi(delta)
    blocks=psi.reshape(NCH,NST).copy()
    maxerr=0.0
    cases=[]
    for ch in range(NCH):
        z=blocks[ch].copy()
        for micro in range(11):
            S=local_action_matrix(z)
            a=S@z
            b=transition_reference(z)
            err=float(np.max(np.abs(a-b)))
            maxerr=max(maxerr,err)
            cases.append({"channel":CHANNELS[ch],"micro":micro,"max_abs_error":err})
            z=b
    return {"pass":bool(maxerr<=1e-12),"max_abs_error":maxerr,"cases":cases}


def audit_blockdiag_control() -> dict:
    delta=0.731
    psi=init_psi(delta)
    ref=psi.reshape(NCH,NST).copy()
    maxerr=0.0
    cases=[]
    for micro in range(11):
        # One exact big-matrix step with cross disabled.
        psi,_=strict_step(psi,cross=False)
        # Independent Paper-7 transitions.
        for ch in range(NCH): ref[ch]=transition_reference(ref[ch])
        err=float(np.max(np.abs(psi-ref.reshape(-1))))
        maxerr=max(maxerr,err)
        cases.append({"micro":micro+1,"max_abs_error":err})
    return {"pass":bool(maxerr<=1e-12),"max_abs_error":maxerr,"cases":cases}


def write_csv(path: Path, rows: list[dict]):
    if not rows: return
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest()


def main():
    audit1=audit_local_matrix()
    audit2=audit_blockdiag_control()
    if not (audit1["pass"] and audit2["pass"]):
        raise RuntimeError("matrix audit failed")

    # Strict one-macrostep scan: 0..2pi, 0.5-degree spacing (721 endpoints).
    deltas=np.linspace(0.0,2.0*math.pi,721,endpoint=True)
    scan=[]
    valid_candidates=[]
    for idx,delta in enumerate(deltas):
        psi0=init_psi(float(delta))
        psif,rows,ok,why,completed=run_microsteps(psi0,11,cross=True)
        b0=psi0.reshape(NCH,NST); bf=psif.reshape(NCH,NST)
        dP=bf[:,P]-b0[:,P]
        entry={
            "grid_index":idx,"delta_rad":float(delta),"delta_over_pi":float(delta/math.pi),
            "completed_microsteps":int(completed),"valid_after_macro":bool(ok),"status":why,
            "P_aa":float(bf[0,P]),"P_ab":float(bf[1,P]),"P_bb":float(bf[2,P]),
            "dP_aa":float(dP[0]),"dP_ab":float(dP[1]),"dP_bb":float(dP[2]),
            "dP_L2":float(np.linalg.norm(dP)),"dP_sum":float(np.sum(dP)),
            "C_aa":float(bf[0,C]),"C_ab":float(bf[1,C]),"C_bb":float(bf[2,C]),
            "qsum_aa":float(np.sum(bf[0,QB:QB+11])),"qsum_ab":float(np.sum(bf[1,QB:QB+11])),"qsum_bb":float(np.sum(bf[2,QB:QB+11])),
        }
        scan.append(entry)
        if ok:
            valid_candidates.append(entry)

    write_csv(OUT/"delta_scan_one_macrostep.csv",scan)
    # Candidates sorted by smallest radial change norm, no physics redefinition.
    valid_candidates=sorted(valid_candidates,key=lambda x:x["dP_L2"])
    top=valid_candidates[:10]
    write_csv(OUT/"valid_candidates_top10.csv",top)

    # Detailed trajectories: exact previously discussed 2pi/3 and best valid candidate, if any.
    detail_specs=[("delta_2pi_over_3",2.0*math.pi/3.0)]
    if top:
        detail_specs.append(("best_valid_grid",float(top[0]["delta_rad"])))
    details=[]
    for label,delta in detail_specs:
        psi0=init_psi(delta)
        psif,rows,ok,why,completed=run_microsteps(psi0,11,cross=True)
        write_csv(OUT/f"microtrajectory_{label}.csv",rows)
        np.save(OUT/f"initial_psi123_{label}.npy",psi0)
        np.save(OUT/f"final_psi123_{label}.npy",psif)
        # Also save the actual full S at micro=0.
        S0,K0=full_action_matrix(psi0,cross=True)
        np.save(OUT/f"S123_micro0_{label}.npy",S0)
        np.save(OUT/f"K3_micro0_{label}.npy",K0)
        details.append({"label":label,"delta_rad":delta,"delta_over_pi":delta/math.pi,"valid":ok,"status":why,"completed_microsteps":completed})

    # If a valid candidate exists, continue it for up to 100 macrosteps, always strict full S.
    longrun=None
    if top:
        delta=float(top[0]["delta_rad"])
        psi=init_psi(delta)
        macro_rows=[]
        status="ok"
        completed_macro=0
        for macro in range(101):
            b=psi.reshape(NCH,NST)
            macro_rows.append({
                "macro":macro,"delta_init":delta,
                "P_aa":float(b[0,P]),"P_ab":float(b[1,P]),"P_bb":float(b[2,P]),
                "C_aa":float(b[0,C]),"C_ab":float(b[1,C]),"C_bb":float(b[2,C]),
                "phase_aa":float(math.atan2(b[0,HI],b[0,HR])) if (b[0,HR]!=0 or b[0,HI]!=0) else math.nan,
                "phase_ab":float(math.atan2(b[1,HI],b[1,HR])) if (b[1,HR]!=0 or b[1,HI]!=0) else math.nan,
                "phase_bb":float(math.atan2(b[2,HI],b[2,HR])) if (b[2,HR]!=0 or b[2,HI]!=0) else math.nan,
            })
            if macro==100: break
            ok,why=validity(psi)
            if not ok:
                status=why; break
            psin,rows,ok,why,completed=run_microsteps(psi,11,cross=True)
            psi=psin
            completed_macro=macro+1
            if not ok:
                status=why; break
        write_csv(OUT/"best_candidate_longrun_100macro.csv",macro_rows)
        np.save(OUT/"best_candidate_longrun_final_psi123.npy",psi)
        longrun={"delta_rad":delta,"delta_over_pi":delta/math.pi,"completed_macrosteps":completed_macro,"status":status}

    summary={
        "experiment":"Paper 8 strict full-state S(Psi)->Psi re-run (8A2)",
        "status":"completed",
        "state_dimension":NTOT,"channel_state_dimension":NST,
        "update_definition":"Psi_next = S_full(Psi) @ Psi; dense 123x123 matrix; synchronous",
        "cross_hypothesis":"k_ij=cos(phi_i-phi_j) from current H phases; same scalar multiplies the entire source local-action block; no normalization",
        "initial_phase_definition":"(phi_aa,phi_ab,phi_bb)=(delta,0,-delta)",
        "initial_CD":{"aa":[0.91,0.0],"ab":[1.09,0.36],"bb":[0.91,0.0]},
        "audits":{"local_matrix_equals_Paper7_transition":audit1,"blockdiag_123_equals_three_independent_Paper7_transitions":audit2},
        "scan":{"grid_points":len(scan),"valid_after_one_macrostep":len(valid_candidates),"invalid_after_one_macrostep":len(scan)-len(valid_candidates),"top_valid":top},
        "details":details,"longrun":longrun,
        "guardrail":"No radiation-amplitude layer, no post-hoc normalization, no clipping, no separate K-after-transition step. All coupled evolution is the single dense S_full(Psi) @ Psi action."
    }
    (OUT/"summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")

    # Text result note.
    lines=[]
    lines.append("# Paper 8 strict re-run 8A2 result\n")
    lines.append("The experiment uses only the single synchronous update `Psi_next = S_full(Psi) @ Psi` with a dense 123x123 state-dependent action.\n")
    lines.append(f"Local-matrix audit max error: {audit1['max_abs_error']:.17g}\n")
    lines.append(f"Block-diagonal 123-state control max error: {audit2['max_abs_error']:.17g}\n")
    lines.append(f"One-macrostep delta scan: {len(valid_candidates)} valid / {len(scan)} total.\n")
    if top:
        t=top[0]
        lines.append(f"Best valid grid point by ||dP||: delta/pi={t['delta_over_pi']:.17g}, ||dP||={t['dP_L2']:.17g}, dP=({t['dP_aa']:.17g},{t['dP_ab']:.17g},{t['dP_bb']:.17g}).\n")
    else:
        lines.append("No delta survived one complete macrostep under the minimal full-state coupling hypothesis.\n")
    for d in details:
        lines.append(f"Detail {d['label']}: delta/pi={d['delta_over_pi']:.17g}, completed_microsteps={d['completed_microsteps']}, valid={d['valid']}, status={d['status']}.\n")
    if longrun:
        lines.append(f"Best-candidate long run: completed_macrosteps={longrun['completed_macrosteps']}, status={longrun['status']}.\n")
    (OUT/"RESULTS_ja.md").write_text("".join(lines),encoding="utf-8")

    # Package hashes.
    files=[Path(__file__)]+sorted([p for p in OUT.iterdir() if p.is_file()])
    manifest={p.name:sha256(p) for p in files}
    (HERE/"SHA256SUMS.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")

    print(json.dumps({
        "audit_local_max":audit1["max_abs_error"],
        "audit_blockdiag_max":audit2["max_abs_error"],
        "valid":len(valid_candidates),"total":len(scan),
        "top":top[:3],"details":details,"longrun":longrun
    },ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
