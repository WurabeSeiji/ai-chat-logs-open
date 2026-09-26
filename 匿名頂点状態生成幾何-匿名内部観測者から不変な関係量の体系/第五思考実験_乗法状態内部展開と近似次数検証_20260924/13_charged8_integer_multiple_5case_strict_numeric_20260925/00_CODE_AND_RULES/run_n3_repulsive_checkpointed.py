
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import argparse, csv, json, math, importlib.util
from pathlib import Path
import numpy as np
import h5py

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("core", HERE/"run_strict_charged8_5cases.py")
core=importlib.util.module_from_spec(spec); spec.loader.exec_module(core)

def save_part(case_id, la, lb, c0, d0, resume=False, part_rows=4_000_000):
    out=HERE/"cases"/case_id
    out.mkdir(parents=True,exist_ok=True)
    cp=out/"checkpoint_state.npz"
    manifest_path=out/"part_manifest.json"

    audit=core.preexec_audit()
    audit["case_initial_state"]={"P":50.0,"N":0.25,"C":c0,"D":d0}
    audit["case_provenance"]={"lambda_A":la,"lambda_B":lb,"C_formula":"1-lambda_A*lambda_B","D_formula":"(lambda_A-lambda_B)^2"}
    audit["checkpoint_note"]="Checkpoint contains exact full 41-component state plus serialization counters only; no new physical state or update rule."
    (out/"preexec_audit.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    if not audit["all_pass"]: raise RuntimeError("audit failed")

    if resume:
        zc=np.load(cp,allow_pickle=False)
        z=zc["z"].astype(np.float64)
        step0=int(zc["next_step"]); part_index=int(zc["part_index"])
        prev_last=zc["last_row"].astype(np.float64)
        first_bad_step=int(zc["first_bad_step"]); first_bad_step=None if first_bad_step<0 else first_bad_step
        first_bad_p=float(zc["first_bad_p"]); first_bad_p=None if not np.isfinite(first_bad_p) else first_bad_p
        maxdn=float(zc["maxdn"]); maxdc=float(zc["maxdc"]); maxdd=float(zc["maxdd"])
    else:
        z=core.init_state(c0=c0,d0=d0)
        step0=0; part_index=0; prev_last=None
        first_bad_step=None; first_bad_p=None
        maxdn=maxdc=maxdd=0.0
        # first macro microstate evidence
        with (out/"raw_first_macro_microsteps.csv").open("w",newline="",encoding="utf-8") as f:
            w=csv.writer(f)
            w.writerow(["micro","q_index","P","N","C","D","k1p","k2p","k3p","k4p","dP","Q"])
            zm=z.copy()
            for m in range(12):
                w.writerow([m,int(np.argmax(zm[core.QB:core.QB+11])),zm[core.P],zm[core.N],zm[core.C],zm[core.D],zm[core.K1],zm[core.K2],zm[core.K3],zm[core.K4],zm[core.DV],zm[core.Q]])
                if m<11: zm=core.transition(zm)

    core.transition(z.copy())  # compile exact transition
    z_after,rows,count,reached=core.generate_chunk(z,step0,part_rows,20.0)
    block=rows[:count]
    cols=["step","P","E","U_re","U_im","H_re","H_im","Q","N","C","D","t_readout","x_readout","y_readout","q_index"]

    part=out/f"raw_macro_part{part_index:03d}.h5"
    with h5py.File(part,"w") as hf:
        ds=hf.create_dataset("raw_macro_trajectory",data=block,dtype="f8",
                             chunks=(min(100_000,count),15),compression="lzf",shuffle=True)
        hf.attrs["columns_json"]=json.dumps(cols)
        hf.attrs["case_id"]=case_id
        hf.attrs["start_step"]=step0
        hf.attrs["row_count"]=count
        hf.attrs["format_note"]="Every macrostep row preserved; serialization/checkpointing is outside transition(z)."

    maxdn=max(maxdn,float(np.max(np.abs(block[:,8]-0.25))))
    maxdc=max(maxdc,float(np.max(np.abs(block[:,9]-c0))))
    maxdd=max(maxdd,float(np.max(np.abs(block[:,10]-d0))))
    if first_bad_step is None:
        bad=np.nonzero(~np.isfinite(block[:,11]))[0]
        if bad.size:
            k=int(bad[0]); first_bad_step=int(block[k,0]); first_bad_p=float(block[k,1])

    parts=[]
    if manifest_path.exists():
        parts=json.loads(manifest_path.read_text(encoding="utf-8")).get("parts",[])
    parts.append({"file":part.name,"start_step":step0,"row_count":count,"sha256":core.sha256(part)})

    if reached:
        last=block[-1].copy()
        prev=block[-2].copy() if count>=2 else prev_last.copy()
        nstep=int(last[0])
        frac=(float(prev[1])-20.0)/(float(prev[1])-float(last[1]))
        t0=float(prev[11]); t1=float(last[11])
        tcross=(t0+frac*(t1-t0)) if math.isfinite(t0) and math.isfinite(t1) else None
        phicross=(nstep-1+frac)*core.HSTEP
        np.save(out/"final_full_state.npy",z_after)
        summary={
          "case_id":case_id,"status":"completed","lambda_A":la,"lambda_B":lb,
          "initial_state":{"P":50.0,"N":0.25,"C":c0,"D":d0},
          "macro_steps":nstep,"microsteps":nstep*11,"final_P":float(last[1]),
          "t_cross_r20":tcross,"phi_cross_r20":phicross,
          "max_identity_drift":{"N":maxdn,"C":maxdc,"D":maxdd},
          "raw_macro_parts":[x["file"] for x in parts],"raw_macro_format":"HDF5 parts",
          "raw_micro_csv":"raw_first_macro_microsteps.csv",
          "first_nonfinite_time_step":first_bad_step,"first_nonfinite_time_P":first_bad_p
        }
        (out/"generator_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        if cp.exists(): cp.unlink()
        status="completed"
    else:
        last=block[-1].copy()
        np.savez(cp,z=z_after,next_step=np.int64(step0+count),part_index=np.int64(part_index+1),last_row=last,
                 first_bad_step=np.int64(-1 if first_bad_step is None else first_bad_step),
                 first_bad_p=np.float64(np.nan if first_bad_p is None else first_bad_p),
                 maxdn=np.float64(maxdn),maxdc=np.float64(maxdc),maxdd=np.float64(maxdd))
        status="in_progress"

    manifest={"case_id":case_id,"status":status,"parts":parts,
              "first_nonfinite_time_step":first_bad_step,"first_nonfinite_time_P":first_bad_p,
              "max_identity_drift":{"N":maxdn,"C":maxdc,"D":maxdd}}
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"case_id":case_id,"status":status,"part":part_index,"rows":count,"last_P":float(block[-1,1]),"first_bad_step":first_bad_step},ensure_ascii=False))
    return status

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--resume",action="store_true")
    args=ap.parse_args()
    save_part("n3_repulsive",0.90,0.90,0.19,0.0,resume=args.resume)

if __name__=="__main__": main()
