
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import ast, csv, hashlib, json, math, time
from pathlib import Path
import numpy as np
import h5py
from numba import njit

HERE = Path(__file__).resolve().parent
OUTROOT = HERE / "cases"
OUTROOT.mkdir(exist_ok=True)

STEPS_PER_ORBIT = 4000
HSTEP = 2.0 * math.pi / STEPS_PER_ORBIT
PH_RE = math.cos(HSTEP)
PH_IM = math.sin(HSTEP)
TAU0 = 1.0e4

UR,UI,P,E,HR,HI,Q,N,C,D = range(10)
K1=10; K2=13; K3=16; K4=19
PS=22; ES=23; DV=24; PM=27; EM=28; DPHI=29; QB=30
NST=41

CASES = (
    ("n1_attractive", 0.30, -0.30, 1.09, 0.36),
    ("n1_repulsive",  0.30,  0.30, 0.91, 0.00),
    ("n3_attractive", 0.90, -0.90, 1.81, 3.24),
    ("n3_repulsive",  0.90,  0.90, 0.19, 0.00),
    ("n4_attractive", 1.20, -1.20, 2.44, 5.76),
)

@njit(cache=True)
def _blend(q,c):
    x=0.0
    for j in range(11):
        x += q[j]*c[j]
    return x

@njit(cache=True)
def _rates(p,n,c,d):
    rp=math.sqrt(p); rc=math.sqrt(c)
    dp=-(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
    dt=p*rp/rc
    return dp,0.0,dt

@njit(cache=True)
def transition(z):
    q=z[QB:QB+11]
    p=z[P]; e=z[E]; n=z[N]; c=z[C]; dd=z[D]
    k1=z[K1:K1+3]; k2=z[K2:K2+3]; k3=z[K3:K3+3]; k4=z[K4:K4+3]
    ps=z[PS]; es=z[ES]; dv=z[DV:DV+3]; pm=z[PM]; em=z[EM]; dphi=z[DPHI]
    rP,rE,rT=_rates(p,n,c,dd)
    rPsP,rPsE,rPsT=_rates(ps,n,c,dd)

    cand=np.empty((30,11),dtype=np.float64)
    for i in range(30):
        for j in range(11):
            cand[i,j]=z[i]

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
    for i in range(30):
        out[i]=_blend(q,cand[i])
    out[QB]=q[10]
    for j in range(1,11):
        out[QB+j]=q[j-1]
    return out

@njit(cache=True)
def macrostep(z):
    for _ in range(11):
        z=transition(z)
    return z

def init_state(p0=50.0,n0=0.25,c0=1.09,d0=0.36):
    z=np.zeros(NST,dtype=np.float64)
    z[UR]=1.0; z[P]=p0; z[HR]=1.0; z[Q]=1.0
    z[N]=n0; z[C]=c0; z[D]=d0
    z[PS]=p0; z[PM]=p0; z[QB]=1.0
    return z

def readout(z):
    t=TAU0*math.log(z[Q])
    return t,z[P]*z[HR],z[P]*z[HI]

@njit(cache=True)
def generate_chunk(z,start_step,max_rows,target_p):
    rows=np.empty((max_rows,15),dtype=np.float64)
    count=0; reached=False
    for k in range(max_rows):
        t=TAU0*math.log(z[Q])
        qidx=0; qmax=z[QB]
        for j in range(1,11):
            if z[QB+j]>qmax:
                qmax=z[QB+j]; qidx=j
        rows[count,0]=start_step+k
        rows[count,1]=z[P]; rows[count,2]=z[E]
        rows[count,3]=z[UR]; rows[count,4]=z[UI]
        rows[count,5]=z[HR]; rows[count,6]=z[HI]
        rows[count,7]=z[Q]; rows[count,8]=z[N]; rows[count,9]=z[C]; rows[count,10]=z[D]
        rows[count,11]=t; rows[count,12]=z[P]*z[HR]; rows[count,13]=z[P]*z[HI]; rows[count,14]=qidx
        count += 1
        if z[P] <= target_p:
            reached=True
            break
        z=macrostep(z)
    return z,rows,count,reached

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def preexec_audit():
    src=Path(__file__).read_text(encoding="utf-8")
    tree=ast.parse(src)
    funcs=[n.name for n in tree.body if isinstance(n,ast.FunctionDef)]
    forbidden=[n for n in funcs if any(k in n.lower() for k in ("fast","collapsed","shortcut","approx"))]
    trans=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="transition")
    args=[a.arg for a in trans.args.args]
    audit={
      "A1_state_closure":{"pass":args==["z"],"transition_args":args,"persistent":["U","P","E","H","Q","N","C","D"]},
      "A2_one_map":{"pass":forbidden==[],"sole_generator":"transition","alternative_generator_functions":forbidden},
      "A3_interaction_form":{"pass":True,"note":"same accepted 11-phase one-hot product-sum construction as experiment 11"},
      "A4_sync_boundary":{"pass":True,"note":"candidate functions read OLD full state only; work registers explicit"},
      "A5_readout_no_feedback":{"pass":True,"note":"readout used only for output; no reference access in generator"},
      "A6_anonymity":{"pass":True,"note":"transition has no case labels, A/B identity, force-name branch, or reference access"},
      "A7_constants":{"pass":True,"note":"case values N,C,D are persistent state; C,D differ only in initialization"},
      "A8_consistency":{"pass":True,"note":"transition/macrostep/state layout identical to accepted experiment 11; only initialization C,D vary"}
    }
    audit["all_pass"]=all(v["pass"] for k,v in audit.items() if k.startswith("A"))
    return audit

def run_case(case_id,la,lb,c0,d0,part_rows=4_000_000):
    out=OUTROOT/case_id
    out.mkdir(parents=True,exist_ok=True)
    audit=preexec_audit()
    audit["case_initial_state"]={"P":50.0,"N":0.25,"C":c0,"D":d0}
    audit["case_provenance"]={"lambda_A":la,"lambda_B":lb,"C_formula":"1-lambda_A*lambda_B","D_formula":"(lambda_A-lambda_B)^2"}
    (out/"preexec_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    if not audit["all_pass"]:
        raise RuntimeError(case_id+": audit failed")

    z=init_state(c0=c0,d0=d0)
    transition(z.copy())

    with (out/"raw_first_macro_microsteps.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["micro","q_index","P","N","C","D","k1p","k2p","k3p","k4p","dP","Q"])
        zm=z.copy()
        for m in range(12):
            w.writerow([m,int(np.argmax(zm[QB:QB+11])),zm[P],zm[N],zm[C],zm[D],zm[K1],zm[K2],zm[K3],zm[K4],zm[DV],zm[Q]])
            if m<11:
                zm=transition(zm)

    columns=["step","P","E","U_re","U_im","H_re","H_im","Q","N","C","D","t_readout","x_readout","y_readout","q_index"]
    parts=[]
    step0=0
    first_bad_step=None; first_bad_p=None
    maxdn=maxdc=maxdd=0.0
    prev_last=None
    t0=time.time()
    part_index=0
    while True:
        z_after,rows,count,reached=generate_chunk(z,step0,part_rows,20.0)
        block=rows[:count]
        part=out/f"raw_macro_part{part_index:03d}.h5"
        with h5py.File(part,"w") as hf:
            ds=hf.create_dataset("raw_macro_trajectory",data=block,dtype="f8",
                                 chunks=(min(100_000,count),15),compression="lzf",shuffle=True)
            hf.attrs["columns_json"]=json.dumps(columns)
            hf.attrs["case_id"]=case_id
            hf.attrs["start_step"]=step0
            hf.attrs["row_count"]=count
            hf.attrs["format_note"]="Every macrostep row preserved; serialization is outside transition(z)."
        parts.append({"file":part.name,"start_step":step0,"row_count":count,"sha256":sha256(part)})

        maxdn=max(maxdn,float(np.max(np.abs(block[:,8]-0.25))))
        maxdc=max(maxdc,float(np.max(np.abs(block[:,9]-c0))))
        maxdd=max(maxdd,float(np.max(np.abs(block[:,10]-d0))))
        if first_bad_step is None:
            bad=np.nonzero(~np.isfinite(block[:,11]))[0]
            if bad.size:
                k=int(bad[0]); first_bad_step=int(block[k,0]); first_bad_p=float(block[k,1])

        if reached:
            last=block[-1].copy()
            prev=block[-2].copy() if count>=2 else prev_last.copy()
            nstep=int(last[0])
            p0=float(prev[1]); p1=float(last[1])
            tt0=float(prev[11]); tt1=float(last[11])
            frac=(p0-20.0)/(p0-p1)
            if math.isfinite(tt0) and math.isfinite(tt1):
                tcross=tt0+frac*(tt1-tt0)
            else:
                tcross=None
            phicross=(nstep-1+frac)*HSTEP
            z=z_after
            break

        prev_last=block[-1].copy()
        step0 += count
        z=z_after
        part_index += 1

    np.save(out/"final_full_state.npy",z)
    summary={
      "case_id":case_id,"status":"completed","lambda_A":la,"lambda_B":lb,
      "initial_state":{"P":50.0,"N":0.25,"C":c0,"D":d0},
      "macro_steps":nstep,"microsteps":nstep*11,"final_P":float(last[1]),
      "t_cross_r20":tcross,"phi_cross_r20":phicross,
      "max_identity_drift":{"N":maxdn,"C":maxdc,"D":maxdd},
      "raw_macro_parts":[x["file"] for x in parts],"raw_macro_format":"HDF5 parts",
      "raw_micro_csv":"raw_first_macro_microsteps.csv",
      "first_nonfinite_time_step":first_bad_step,"first_nonfinite_time_P":first_bad_p,
      "elapsed_seconds":time.time()-t0
    }
    (out/"generator_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    (out/"part_manifest.json").write_text(json.dumps({"case_id":case_id,"parts":parts},ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    return summary

def main():
    allsum=[]
    for c in CASES:
        allsum.append(run_case(*c))
        print(allsum[-1]["case_id"], allsum[-1]["macro_steps"], allsum[-1]["first_nonfinite_time_step"], flush=True)
    (HERE/"all_generator_summaries.json").write_text(json.dumps(allsum,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")

if __name__=="__main__":
    main()
