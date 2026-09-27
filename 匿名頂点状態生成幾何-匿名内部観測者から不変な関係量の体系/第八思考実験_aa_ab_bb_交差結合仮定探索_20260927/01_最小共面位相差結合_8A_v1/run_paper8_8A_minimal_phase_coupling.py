#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 8 / Experiment 8A
Minimal coplanar phase-difference coupling scan for relation states (aa, ab, bb).

Hypotheses tested (not derived):
  H1: relation states are represented as aa=a^2, ab=a*b, bb=b^2.
  H2: their orbit normals are identical (coplanar, same orientation).
  H3: phi_a=+delta/2, phi_b=-delta/2, hence
      phi_aa=+delta, phi_ab=0, phi_bb=-delta.
  H4: real symmetric cross-coupling matrix
      K(delta) = [[1, cos d, cos 2d],
                  [cos d, 1, cos d],
                  [cos 2d, cos d, 1]].
  H5: cancellation is tested at the radiation-amplitude layer.
      GW phase follows 2*orbital phase (quadrupole), EM dipole phase follows orbital phase.

No feedback to the Paper-7 local transition is introduced in this first screen.
This is intentionally a hypothesis-screening experiment, not yet the full coupled dynamics.
"""
from __future__ import annotations
import csv, json, hashlib, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

# Paper 6 supplement / Paper 7 normalization and n=1 source pattern.
G = 1.0
Q0 = 1.0
C_LIGHT = 1.0
M_OVER_Q0 = 20.0/3.0
N = 0.25
P0 = 50.0

# n=1, opposite-sign ab source; aa and bb are self terms.
CAA, DAA = 0.91, 0.0
CAB, DAB = 1.09, 0.36
CBB, DBB = 0.91, 0.0

RADII = np.array([50.0, 40.0, 30.0, 20.0], dtype=np.float64)
GRID_N = 360001  # includes 2*pi/3 and 4*pi/3 exactly because GRID_N-1=360000.

# Paper-7 logical full microstate layout, copied only for reproducible 123-state initialization.
NST = 41
UR,UI,P,E,HR,HI,Q,NIDX,C,D = range(10)
K1=10; K2=13; K3=16; K4=19
PS=22; ES=23; DV=24; PM=27; EM=28; DPHI=29; QB=30


def init_channel(c0: float, d0: float, phi0: float) -> np.ndarray:
    z = np.zeros(NST, dtype=np.float64)
    z[UR] = 1.0
    z[P] = P0
    z[HR] = math.cos(phi0)
    z[HI] = math.sin(phi0)
    z[Q] = 1.0
    z[NIDX] = N
    z[C] = c0
    z[D] = d0
    z[PS] = P0
    z[PM] = P0
    z[QB] = 1.0
    return z


def initial_state_123(delta: float) -> np.ndarray:
    # Definition-level phase assignment from H1/H3.
    ph = np.array([delta, 0.0, -delta], dtype=np.float64)
    zaa = init_channel(CAA, DAA, ph[0])
    zab = init_channel(CAB, DAB, ph[1])
    zbb = init_channel(CBB, DBB, ph[2])
    return np.concatenate([zaa, zab, zbb])


def coupling_matrix(delta: float) -> np.ndarray:
    c1 = math.cos(delta)
    c2 = math.cos(2.0*delta)
    return np.array([[1.0,c1,c2],[c1,1.0,c1],[c2,c1,1.0]], dtype=np.float64)


def powers(r: float):
    cs = np.array([CAA, CAB, CBB], dtype=np.float64)
    ds = np.array([DAA, DAB, DBB], dtype=np.float64)
    pgw = (32.0/5.0)*(N*N)*(cs**3)/(r**5)
    pem = (2.0/3.0)*(N*N)*ds*(cs**2)/(r**4)
    return pgw, pem


def radiation_residuals(delta: float, r: float):
    K = coupling_matrix(delta)
    ph = np.array([delta, 0.0, -delta], dtype=np.float64)
    pgw, pem = powers(r)
    xgw = np.sqrt(pgw) * np.exp(2j*ph)   # quadrupole phase convention H5
    xem = np.sqrt(pem) * np.exp(1j*ph)   # dipole phase convention H5
    ygw = K @ xgw
    yem = K @ xem
    Agw = ygw.sum()
    Aem = yem.sum()
    # Two normalizations: source amplitude sum and post-coupling amplitude sum.
    gw_src = float(np.sum(np.abs(xgw)))
    em_src = float(np.sum(np.abs(xem)))
    gw_post = float(np.sum(np.abs(ygw)))
    em_post = float(np.sum(np.abs(yem)))
    rgw_src = abs(Agw)/(gw_src + 1e-300)
    rem_src = abs(Aem)/(em_src + 1e-300)
    rgw_post = abs(Agw)/(gw_post + 1e-300)
    rem_post = abs(Aem)/(em_post + 1e-300)
    joint = math.hypot(rgw_post, rem_post)
    return {
        'A_GW_re':float(Agw.real),'A_GW_im':float(Agw.imag),'A_GW_abs':float(abs(Agw)),
        'A_EM_re':float(Aem.real),'A_EM_im':float(Aem.imag),'A_EM_abs':float(abs(Aem)),
        'GW_resid_source_norm':float(rgw_src),'EM_resid_source_norm':float(rem_src),
        'GW_resid_post_norm':float(rgw_post),'EM_resid_post_norm':float(rem_post),
        'joint_post_norm':float(joint),
        'K_row_sum_0':float(K[0].sum()),'K_row_sum_1':float(K[1].sum()),'K_row_sum_2':float(K[2].sum()),
    }


def apply_K_componentwise_to_123(z123: np.ndarray, delta: float) -> np.ndarray:
    blocks = z123.reshape(3, NST)
    return (coupling_matrix(delta) @ blocks).reshape(3*NST)


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):
            h.update(b)
    return h.hexdigest()


def main():
    # Dense deterministic scan (vectorized; total amplitude uses column sums of K).
    deltas = np.linspace(0.0, 2.0*math.pi, GRID_N, endpoint=True)
    c1=np.cos(deltas); c2=np.cos(2.0*deltas)
    s0=1.0+c1+c2
    s1=1.0+2.0*c1
    s2=s0.copy()
    ph0=deltas; ph1=np.zeros_like(deltas); ph2=-deltas
    scan_path = OUT/'delta_scan.csv'
    best = {float(r): {'joint':(float('inf'),None), 'gw':(float('inf'),None), 'em':(float('inf'),None)} for r in RADII}
    all_rows=[]
    for r in RADII:
        pgw,pem=powers(float(r))
        x0=np.sqrt(pgw[0])*np.exp(2j*ph0)
        x1=np.sqrt(pgw[1])*np.exp(2j*ph1)
        x2=np.sqrt(pgw[2])*np.exp(2j*ph2)
        Agw=s0*x0+s1*x1+s2*x2
        # EM source exists only in ab for this n=1 self/other/self setup.
        e1=np.sqrt(pem[1])*np.exp(1j*ph1)
        Aem=s1*e1
        # Post-coupling normalization needs y=Kx. Vectorized explicitly.
        y0=x0+c1*x1+c2*x2
        y1=c1*x0+x1+c1*x2
        y2=c2*x0+c1*x1+x2
        gwpost=np.abs(y0)+np.abs(y1)+np.abs(y2)
        # EM: K @ [0,e1,0] = [c1 e1,e1,c1 e1]
        empost=np.abs(c1*e1)+np.abs(e1)+np.abs(c1*e1)
        rgw=np.abs(Agw)/(gwpost+1e-300)
        rem=np.abs(Aem)/(empost+1e-300)
        joint=np.hypot(rgw,rem)
        for key,arr in [('joint',joint),('gw',rgw),('em',rem)]:
            idx=int(np.argmin(arr))
            best[float(r)][key]=(float(arr[idx]),float(deltas[idx]))
        all_rows.append(np.column_stack([deltas,deltas/math.pi,np.full_like(deltas,r),rgw,rem,joint]))
    data=np.vstack(all_rows)
    np.savetxt(scan_path,data,delimiter=',',header='delta_rad,delta_over_pi,radius,GW_resid_post_norm,EM_resid_post_norm,joint_post_norm',comments='')

    # Exact candidate values implied by the scan/grid.
    exact_candidates=[]
    for label,delta in [('2pi_over_3',2.0*math.pi/3.0),('4pi_over_3',4.0*math.pi/3.0)]:
        K=coupling_matrix(delta)
        z=initial_state_123(delta)
        zK=apply_K_componentwise_to_123(z,delta)
        # Componentwise channel-sum residual of the 123-state coupling layer.
        sum_before=z.reshape(3,NST).sum(axis=0)
        sum_after=zK.reshape(3,NST).sum(axis=0)
        item={
            'label':label,'delta_rad':delta,'delta_over_pi':delta/math.pi,
            'K':K.tolist(),'K_row_sums':K.sum(axis=1).tolist(),'K_eigenvalues':np.linalg.eigvalsh(K).tolist(),
            'initial_123_nonzero_count':int(np.count_nonzero(z)),
            'channel_sum_before_L2':float(np.linalg.norm(sum_before)),
            'channel_sum_after_L2':float(np.linalg.norm(sum_after)),
            'channel_sum_after_max_abs':float(np.max(np.abs(sum_after))),
            'radii':{}
        }
        for r in RADII:
            item['radii'][str(float(r))]=radiation_residuals(delta,float(r))
        exact_candidates.append(item)
        np.save(OUT/f'initial_state_123_{label}.npy', z)
        np.save(OUT/f'coupled_state_123_{label}.npy', zK)

    # Human-readable initial nonzero entries at the first exact candidate.
    delta=2.0*math.pi/3.0
    z=initial_state_123(delta).reshape(3,NST)
    names=['aa','ab','bb']
    with (OUT/'initial_state_nonzero_2pi_over_3.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['channel','component_index','value'])
        for ci,name in enumerate(names):
            for idx,val in enumerate(z[ci]):
                if val != 0.0:
                    w.writerow([name,idx,repr(float(val))])

    summary={
        'experiment':'Paper 8 / 8A minimal coplanar phase-difference coupling',
        'status':'completed',
        'normalization':{'G':G,'q0':Q0,'c':C_LIGHT,'M_over_q0':M_OVER_Q0,'N':N},
        'n':1,
        'channel_initial_relations':{
            'aa':{'C':CAA,'D':DAA},'ab':{'C':CAB,'D':DAB},'bb':{'C':CBB,'D':DBB}
        },
        'hypotheses':[
            'aa:=a^2, ab:=a*b, bb:=b^2',
            'orbit normals n_aa=n_ab=n_bb',
            'phi_a=+delta/2, phi_b=-delta/2 -> phi_aa=delta, phi_ab=0, phi_bb=-delta',
            'K_ii=1, K_aa_ab=K_ab_bb=cos(delta), K_aa_bb=cos(2delta), symmetric',
            'first screen applies K at radiation-amplitude layer; no feedback into Paper-7 local transition',
            'GW radiation phase=2*relation orbital phase; EM dipole phase=relation orbital phase'
        ],
        'scan':{'delta_min':0.0,'delta_max':2.0*math.pi,'grid_points':GRID_N,'radii':RADII.tolist()},
        'best_by_radius':{
            str(r):{
                k:{'residual':float(v[0]),'delta_rad':float(v[1]),'delta_over_pi':float(v[1]/math.pi)}
                for k,v in best[r].items()
            } for r in best
        },
        'exact_candidates':exact_candidates,
        'interpretation_guardrail':'A zero here establishes only that this explicit trial coupling has a cancellation sector; it does not establish that the coupling law is physical or unique.'
    }
    (OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

    # Reproducibility manifest.
    files=[Path(__file__), scan_path, OUT/'summary.json', OUT/'initial_state_nonzero_2pi_over_3.csv']
    for label in ['2pi_over_3','4pi_over_3']:
        files += [OUT/f'initial_state_123_{label}.npy', OUT/f'coupled_state_123_{label}.npy']
    manifest={p.name:sha256(p) for p in files}
    (OUT/'SHA256SUMS.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

    print(json.dumps({
        'best_by_radius':summary['best_by_radius'],
        'exact_2pi_over_3':exact_candidates[0],
        'exact_4pi_over_3':exact_candidates[1]
    },ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
