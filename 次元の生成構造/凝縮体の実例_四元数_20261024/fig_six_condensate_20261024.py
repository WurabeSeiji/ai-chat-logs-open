"""図化（SVG と PNG）：six_condensate_results_20261024.json を読む
python fig_six_condensate_20261024.py → six_condensate_20261024.svg / .png
"""
import json, numpy as np

res=json.load(open("six_condensate_results_20261024.json")); r=res["six_condensate"]
names=["A1","A2","A3","B1","B2","B3"]; A={}
for p in r["pairs"]:
    a,b=p["pair"]; A[(a,b)]=A[(b,a)]=p["alpha_deg"]

# ---- SVG（軽量・手組み）
def col(a):
    t=a/180
    if t<.5: c=(59+(221-59)*2*t,76+(221-76)*2*t,192+(221-192)*2*t)
    else:    c=(221+(180-221)*(2*t-1),221+(4-221)*(2*t-1),221+(38-221)*(2*t-1))
    return "rgb(%d,%d,%d)"%c
cs,x0,y0=50,50,70
o=['<svg xmlns="http://www.w3.org/2000/svg" width="820" height="440" font-family="Noto Sans CJK JP,sans-serif" font-size="12">','<rect width="100%" height="100%" fill="white"/>',
   f'<text x="{x0+3*cs}" y="30" text-anchor="middle" font-size="14">四元数 6 個の相対位相の角 α</text>',
   f'<text x="{x0+3*cs}" y="50" text-anchor="middle">灰＝90°（|k|=1）、赤ほど結合が深い、青寄りは自由</text>']
for i,a in enumerate(names):
    o.append(f'<text x="{x0-8}" y="{y0+i*cs+cs/2+4}" text-anchor="end">{a}</text>')
    o.append(f'<text x="{x0+i*cs+cs/2}" y="{y0+6*cs+18}" text-anchor="middle">{a}</text>')
    for j,b in enumerate(names):
        fill="rgb(255,255,255)" if i==j else col(A[(a,b)])
        o.append(f'<rect x="{x0+j*cs}" y="{y0+i*cs}" width="{cs}" height="{cs}" fill="{fill}" stroke="#888"/>')
        if i!=j: o.append(f'<text x="{x0+j*cs+cs/2}" y="{y0+i*cs+cs/2+4}" text-anchor="middle">{A[(a,b)]:.1f}°</text>')
o.append(f'<rect x="{x0}" y="{y0}" width="{3*cs}" height="{3*cs}" fill="none" stroke="black" stroke-width="2"/>')
o.append(f'<rect x="{x0+3*cs}" y="{y0+3*cs}" width="{3*cs}" height="{3*cs}" fill="none" stroke="black" stroke-width="2"/>')
mm=res["maxmin_angle_under_closures"]; tx=420
lines=["閉塞 Σq=0, Σq²=0：成立",f"虚部の特異値：{', '.join(f'{v:.3f}' for v in r['vector_part_singular_values'])}（3 方向）",
       "A 内・B 内（6 組）：α=75.52°, |k|=0.775, D=1.291（自由）","A–B 間（9 組）：α=93.8°〜158.9°（すべて結合）",
       f"Δ₃：20 組すべて正（最小 {min(t['Delta3'] for t in r['triples']):.4f}）","",
       "閉塞のもとでの最小角の最大値（全組の結合の上限）："]+[f"  N={k}：{v}°" for k,v in mm.items()]+["  → N=6 で全組 |k|≥1 は不可能。本例は上限 75.52° に一致"]
for k,l in enumerate(lines): o.append(f'<text x="{tx}" y="{90+k*22}">{l}</text>')
o.append('</svg>')
open("six_condensate_20261024.svg","w").write("\n".join(o))

# ---- PNG（matplotlib）
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt; from matplotlib import rcParams
rcParams["font.family"]=["Noto Sans CJK JP","DejaVu Sans"]
M=np.zeros((6,6))
for (a,b),v in A.items(): M[names.index(a),names.index(b)]=v
fig,ax=plt.subplots(figsize=(6.4,5.6))
im=ax.pcolormesh(np.arange(7)-.5,np.arange(7)-.5,M,vmin=0,vmax=180,cmap="coolwarm"); ax.invert_yaxis(); ax.set_aspect("equal")
for i in range(6):
    for j in range(6):
        if i!=j: ax.text(j,i,f"{M[i,j]:.1f}°",ha="center",va="center",fontsize=9)
ax.set_xticks(range(6)); ax.set_yticks(range(6)); ax.set_xticklabels(names); ax.set_yticklabels(names)
for o_ in [-.5,2.5]: ax.add_patch(plt.Rectangle((o_,o_),3,3,fill=False,lw=2))
ax.set_title("四元数 6 個（複素数 12 個）の相対位相の角 α\nA 内・B 内は 75.52°（自由）、A–B 間はすべて 90° 超（結合）")
fig.colorbar(im,label="α（90°＝|k|=1）"); fig.savefig("six_condensate_20261024.png",dpi=150,bbox_inches="tight")
