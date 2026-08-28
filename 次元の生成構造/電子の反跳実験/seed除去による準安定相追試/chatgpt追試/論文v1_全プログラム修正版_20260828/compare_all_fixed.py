# -*- coding: utf-8 -*-
"""original/<pkg> の全データ・図ファイルを rerun/<pkg> の同名ファイルと突合する。
 csv: 数値列は max|Δ|・max 相対差、文字列列は不一致セル数、形状不一致は記録
 json: 数値葉の max|Δ|、その他は等値
 png: ピクセル不一致率（matplotlib のバージョン差で 0 にはならないので参考値）
 生成判定: rerun 側の mtime が results/run_start_epoch.txt より新しいものだけ「再生成」。古いものは「入力コピー／未生成」
出力: results/comparison.json, results/comparison_table.md"""
import os, json, math, csv, io, numpy as np, pandas as pd
from PIL import Image
HERE=os.path.dirname(os.path.abspath(__file__)); O=os.path.join(os.path.dirname(HERE),"論文v1_全再現テスト_20260828","original"); R=os.path.join(HERE,"fixed")
T0=float(open(os.path.join(HERE,"results","run_start_epoch.txt")).read())
SKIP=(".py",".md",".txt",".pyc",".cpp",".patch")
ALT={ "N3_N16_partial_zero_closure_analysis_20260826": "out" }  # 出力先を変えたパッケージ（rmtree 回避）
def cmp_csv(a,b):
    A=pd.read_csv(a); B=pd.read_csv(b); rec={"shape_o":list(A.shape),"shape_r":list(B.shape)}
    if A.shape!=B.shape: rec["verdict"]="SHAPE_DIFF"; return rec
    cols=[c for c in A.columns if c in B.columns]; rec["cols_missing"]=[c for c in A.columns if c not in B.columns]
    mx=0.0; rel=0.0; sdiff=0; ncells=0
    for c in cols:
        x,y=A[c],B[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            d=np.abs(x.to_numpy(float)-y.to_numpy(float)); d=d[~np.isnan(d)]
            if len(d): mx=max(mx,float(d.max())); s=np.abs(x.to_numpy(float)); m=s>1e-300; 
            if len(d) and m.any(): rel=max(rel,float((d[m[:len(d)]]/s[m]).max()) if len(d)==len(s) else rel)
        else:
            sdiff+=int((x.astype(str)!=y.astype(str)).sum()); ncells+=len(x)
    rec.update(max_abs_diff=mx,max_rel_diff=rel,str_mismatch=sdiff,str_cells=ncells)
    rec["verdict"]="IDENTICAL" if mx==0 and sdiff==0 else ("ROUNDOFF(<1e-8)" if mx<1e-8 and sdiff==0 else "DIFF")
    return rec
def leaves(o,p=""):
    if isinstance(o,dict):
        for k,v in o.items(): yield from leaves(v,p+"/"+str(k))
    elif isinstance(o,list):
        for i,v in enumerate(o): yield from leaves(v,p+f"[{i}]")
    else: yield p,o
def cmp_json(a,b):
    A=dict(leaves(json.load(open(a)))); B=dict(leaves(json.load(open(b)))); mx=0.0; other=0; miss=0
    for k,v in A.items():
        if k not in B: miss+=1; continue
        w=B[k]
        if isinstance(v,(int,float)) and isinstance(w,(int,float)) and not isinstance(v,bool): mx=max(mx,abs(float(v)-float(w)))
        elif v!=w: other+=1
    rec={"max_abs_diff":mx,"nonnumeric_mismatch":other,"keys_missing":miss}
    rec["verdict"]="IDENTICAL" if mx==0 and other==0 and miss==0 else ("ROUNDOFF(<1e-8)" if mx<1e-8 and other==0 and miss==0 else "DIFF"); return rec
def cmp_png(a,b):
    A=np.asarray(Image.open(a).convert("RGB")); B=np.asarray(Image.open(b).convert("RGB"))
    if A.shape!=B.shape: return {"shape_o":list(A.shape),"shape_r":list(B.shape),"verdict":"SHAPE_DIFF"}
    fr=float((np.abs(A.astype(int)-B.astype(int)).sum(axis=2)>0).mean()); return {"pixel_mismatch_fraction":fr,"verdict":"IDENTICAL" if fr==0 else ("NEAR(<2%)" if fr<0.02 else "DIFF")}
out=[]
for pkg in sorted(os.listdir(O)):
    po=os.path.join(O,pkg)
    if not os.path.isdir(po): continue
    for root,_,files in os.walk(po):
        for f in sorted(files):
            if f.endswith(SKIP) or "__pycache__" in root: continue
            relp=os.path.relpath(os.path.join(root,f),po)
            cand=[os.path.join(R,pkg,relp)]
            if pkg in ALT: cand.insert(0,os.path.join(R,pkg,ALT[pkg],relp))
            rr=next((c for c in cand if os.path.exists(c)),None)
            rec={"package":pkg,"file":relp}
            if rr is None: rec["verdict"]="MISSING_IN_RERUN"; out.append(rec); continue
            regenerated=os.path.getmtime(rr)>T0
            if not regenerated: rec["verdict"]="NOT_REGENERATED(no program in package)"; out.append(rec); continue
            try:
                if f.endswith(".csv"): rec.update(cmp_csv(os.path.join(root,f),rr))
                elif f.endswith(".json"): rec.update(cmp_json(os.path.join(root,f),rr))
                elif f.endswith(".png"): rec.update(cmp_png(os.path.join(root,f),rr))
                else: rec["verdict"]="IDENTICAL" if open(os.path.join(root,f),"rb").read()==open(rr,"rb").read() else "DIFF"
            except Exception as e: rec["verdict"]=f"ERROR: {e}"
            out.append(rec)
json.dump(out,open(os.path.join(HERE,"results","comparison.json"),"w"),indent=1,ensure_ascii=False)
from collections import Counter
lines=["# 再現テスト突合表（original（旧エンジン） vs fixed（修正版 FIX1-4））","",f"ファイル数 {len(out)}: "+", ".join(f"{k} {v}" for k,v in sorted(Counter(r['verdict'].split('(')[0] for r in out).items())),""]
lines+=["| package | file | verdict | max|Δ| / pixel mismatch | 備考 |","|---|---|---|---|---|"]
for r in out:
    v=r.get("max_abs_diff", r.get("pixel_mismatch_fraction","")); v=f"{v:.3g}" if isinstance(v,float) else v
    note=[]; 
    if r.get("shape_o") and r.get("shape_o")!=r.get("shape_r"): note.append(f"shape {r['shape_o']}→{r['shape_r']}")
    if r.get("str_mismatch"): note.append(f"str mismatch {r['str_mismatch']}/{r['str_cells']}")
    if r.get("nonnumeric_mismatch"): note.append(f"nonnumeric {r['nonnumeric_mismatch']}")
    if r.get("keys_missing"): note.append(f"keys missing {r['keys_missing']}")
    lines.append(f"| {r['package']} | {r['file']} | {r['verdict']} | {v} | {'; '.join(note)} |")
open(os.path.join(HERE,"results","comparison_table.md"),"w",encoding="utf-8").write("\n".join(lines)+"\n"); print("\n".join(lines[:3]))
