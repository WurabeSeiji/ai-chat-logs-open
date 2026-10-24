"""図化（SVG と PNG）：octonion_condensate_results_20261024.json を読む
python fig_octonion_condensate_20261024.py → octonion_condensate_20261024.svg / .png
"""
import json, numpy as np
res=json.load(open("octonion_condensate_results_20261024.json")); r=res["octonion_condensate_N10"]
names=[f"A{i}" for i in range(1,6)]+[f"B{i}" for i in range(1,6)]; M=np.zeros((10,10))
for p in r["pairs"]:
    a,b=[names.index(x) for x in p["pair"]]; M[a,b]=M[b,a]=p["alpha_deg"]
def col(a):
    t=a/180
    if t<.5: c=(59+(221-59)*2*t,76+(221-76)*2*t,192+(221-192)*2*t)
    else:    c=(221+(180-221)*(2*t-1),221+(4-221)*(2*t-1),221+(38-221)*(2*t-1))
    return "rgb(%d,%d,%d)"%c
cs,x0,y0=38,45,70
o=['<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520" font-family="Noto Sans CJK JP,sans-serif" font-size="11">','<rect width="100%" height="100%" fill="white"/>',
   f'<text x="{x0+5*cs}" y="30" text-anchor="middle" font-size="14">八元数 10 個（複素数 40 個）の相対位相の角 α</text>',
   f'<text x="{x0+5*cs}" y="50" text-anchor="middle">灰＝90°（|k|=1）、赤寄りは結合、青寄りは自由</text>']
for i,a in enumerate(names):
    o.append(f'<text x="{x0-6}" y="{y0+i*cs+cs/2+4}" text-anchor="end">{a}</text>')
    o.append(f'<text x="{x0+i*cs+cs/2}" y="{y0+10*cs+15}" text-anchor="middle">{a}</text>')
    for j in range(10):
        fill="rgb(255,255,255)" if i==j else col(M[i,j])
        o.append(f'<rect x="{x0+j*cs}" y="{y0+i*cs}" width="{cs}" height="{cs}" fill="{fill}" stroke="#999"/>')
        if i!=j: o.append(f'<text x="{x0+j*cs+cs/2}" y="{y0+i*cs+cs/2+4}" text-anchor="middle" font-size="10">{M[i,j]:.0f}°</text>')
for s in [0,5]: o.append(f'<rect x="{x0+s*cs}" y="{y0+s*cs}" width="{5*cs}" height="{5*cs}" fill="none" stroke="black" stroke-width="2"/>')
mm=res.get("maxmin_angle_under_closures_8D",{}); d3=r["Delta3_min_max_count_pos"]; asc=r["associator_min_max_count_pos"]
lines=["閉塞 Σo=0, Σo²=0：成立","虚部の特異値：7 個すべて正（7 方向に張る）",
       f"群内（20 組）：α={r['within_alpha_range'][0]:.2f}°（自由）",
       f"群間（25 組）：α={r['cross_alpha_range'][0]:.1f}°〜{r['cross_alpha_range'][1]:.1f}°（すべて結合）",
       f"三体の非可換性 Δ₃：{d3[2]}/{d3[3]} 組が正（最小 {d3[0]:.3f}）",
       f"四体の非結合性 |(uv)w−u(vw)|：{asc[2]}/{asc[3]} 組が正（最小 {asc[0]:.3f}）","",
       "閉塞のもとでの最小角の最大値（8 次元）："]+[f"  N={k}：{v}°" for k,v in mm.items()]+["  → 本例は N=10 の上限に一致"]
for k,l in enumerate(lines): o.append(f'<text x="470" y="{90+k*22}" font-size="12">{l}</text>')
o.append('</svg>'); open("octonion_condensate_20261024.svg","w").write("\n".join(o))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib import rcParams
rcParams["font.family"]=["Noto Sans CJK JP","DejaVu Sans"]
fig,ax=plt.subplots(figsize=(8,7))
im=ax.pcolormesh(np.arange(11)-.5,np.arange(11)-.5,M,vmin=0,vmax=180,cmap="coolwarm"); ax.invert_yaxis(); ax.set_aspect("equal")
for i in range(10):
    for j in range(10):
        if i!=j: ax.text(j,i,f"{M[i,j]:.0f}°",ha="center",va="center",fontsize=8)
ax.set_xticks(range(10)); ax.set_yticks(range(10)); ax.set_xticklabels(names); ax.set_yticklabels(names)
for s in [-.5,4.5]: ax.add_patch(plt.Rectangle((s,s),5,5,fill=False,lw=2))
ax.set_title("八元数 10 個（複素数 40 個）の相対位相の角 α\n群内 67.98°（自由）、群間はすべて 90° 超（結合）")
fig.colorbar(im,label="α（90°＝|k|=1）"); fig.savefig("octonion_condensate_20261024.png",dpi=150,bbox_inches="tight")
print("ok")
