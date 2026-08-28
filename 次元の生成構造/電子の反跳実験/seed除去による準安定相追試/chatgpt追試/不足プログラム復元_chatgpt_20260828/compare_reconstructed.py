# -*- coding: utf-8 -*-
"""ChatGPT 復元プログラムの出力（verify_root）を、論文v1_全再現テスト の original と突合。対象＝comparison.json の NOT_REGENERATED 101 件。"""
import os, json, numpy as np, pandas as pd
from PIL import Image
H=os.path.dirname(os.path.abspath(__file__)); O=os.path.join(os.path.dirname(H),"論文v1_全再現テスト_20260828","original"); V=os.path.join(H,"verify_root"); T0=float(open(os.path.join(H,"results","run_start_epoch.txt")).read())
miss=[r for r in json.load(open(os.path.join(os.path.dirname(H),"論文v1_全再現テスト_20260828","results","comparison.json"))) if r["verdict"].startswith("NOT_REGENERATED")]
def cmp_csv(a,b):
    A=pd.read_csv(a); B=pd.read_csv(b)
    if A.shape!=B.shape or list(A.columns)!=list(B.columns): return "SHAPE_DIFF",f"{A.shape}->{B.shape}"
    mx=0.0; sd=0
    for c in A.columns:
        if pd.api.types.is_numeric_dtype(A[c]) and pd.api.types.is_numeric_dtype(B[c]):
            d=np.abs(A[c].to_numpy(float)-B[c].to_numpy(float)); d=d[~np.isnan(d)]; mx=max(mx,float(d.max()) if len(d) else 0.0)
        else: sd+=int((A[c].astype(str)!=B[c].astype(str)).sum())
    return ("IDENTICAL" if mx==0 and sd==0 else "ROUNDOFF" if mx<1e-8 and sd==0 else "DIFF"), f"max|Δ|={mx:.2e} str_mismatch={sd}"
def leaves(o,p=""):
    if isinstance(o,dict):
        for k,v in o.items(): yield from leaves(v,p+"/"+str(k))
    elif isinstance(o,list):
        for i,v in enumerate(o): yield from leaves(v,p+f"[{i}]")
    else: yield p,o
def cmp_json(a,b):
    A=dict(leaves(json.load(open(a)))); B=dict(leaves(json.load(open(b)))); mx=0.0; oth=0; miss=0
    for k,v in A.items():
        if k not in B: miss+=1; continue
        w=B[k]
        if isinstance(v,(int,float)) and not isinstance(v,bool) and isinstance(w,(int,float)): mx=max(mx,abs(float(v)-float(w)))
        elif v!=w: oth+=1
    return ("IDENTICAL" if mx==0 and oth==0 and miss==0 else "ROUNDOFF" if mx<1e-8 and oth==0 and miss==0 else "DIFF"), f"max|Δ|={mx:.2e} other={oth} missing={miss}"
def cmp_png(a,b):
    A=np.asarray(Image.open(a).convert("RGB")); B=np.asarray(Image.open(b).convert("RGB"))
    if A.shape!=B.shape: return "SHAPE_DIFF",f"{A.shape}->{B.shape}"
    fr=float((np.abs(A.astype(int)-B.astype(int)).sum(axis=2)>0).mean()); return ("IDENTICAL" if fr==0 else "NEAR" if fr<0.02 else "DIFF"), f"pixel_mismatch={fr:.3f}"
rows=[]
for r in miss:
    o=os.path.join(O,r["package"],r["file"]); v=os.path.join(V,r["package"],r["file"])
    if not os.path.exists(v): rows.append((r["package"],r["file"],"MISSING","")); continue
    if os.path.getmtime(v)<=T0: rows.append((r["package"],r["file"],"NOT_REGENERATED","（復元プログラムはこのファイルを書かなかった）")); continue
    try:
        f=r["file"]; res=cmp_csv(o,v) if f.endswith(".csv") else cmp_json(o,v) if f.endswith(".json") else cmp_png(o,v) if f.endswith(".png") else (("IDENTICAL" if open(o,"rb").read()==open(v,"rb").read() else "DIFF"),"bytes")
    except Exception as e: res=("ERROR",str(e)[:80])
    rows.append((r["package"],r["file"],res[0],res[1]))
from collections import Counter; c=Counter(x[2] for x in rows)
L=["# ChatGPT 復元プログラムの検証結果（original との突合）","",f"対象 {len(rows)} 件: "+", ".join(f"{k} {v}" for k,v in sorted(c.items())),"","| package | file | verdict | detail |","|---|---|---|---|"]+[f"| {p} | {f} | {v} | {d} |" for p,f,v,d in rows]
open(os.path.join(H,"results","compare_reconstructed.md"),"w",encoding="utf-8").write("\n".join(L)+"\n"); print("\n".join(L[:3]))
for p,f,v,d in rows:
    if v not in ("IDENTICAL","ROUNDOFF","NEAR"): print(f"  {v:15s} {p[:38]:38s} {f:55s} {d}")
