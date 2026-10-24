"""図化：condensate_example_results_20261024.json と Δ3(θ) 曲線"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
import importlib.util, sys
rcParams["font.family"] = ["Noto Sans CJK JP","DejaVu Sans"]; rcParams["svg.fonttype"] = "none"
spec = importlib.util.spec_from_file_location("ce","condensate_example_20261024.py")
ce = importlib.util.module_from_spec(spec); sys.argv=["x","/tmp/_r.json"]; spec.loader.exec_module(ce)
res = json.load(open("condensate_example_results_20261024.json"))
fig, axes = plt.subplots(1,4, figsize=(16,4))
for ax,key,title in zip(axes[:3],
    ["C_degenerate_theta0","intermediate_theta45","H_condensate_theta90"],
    ["θ=0°（m∥n：複素数）","θ=45°","θ=90°（m⊥n：四元数）"]):
    M = np.zeros((4,4))
    for p in res[key]["pairs"]:
        a,b = p["pair"]; M[a-1,b-1]=M[b-1,a-1]=p["alpha_deg"]
    im = ax.pcolormesh(np.arange(5)-0.5, np.arange(5)-0.5, M, vmin=0, vmax=180, cmap="coolwarm"); ax.invert_yaxis(); ax.set_aspect("equal")
    for i in range(4):
        for j in range(4):
            if i!=j: ax.text(j,i,f"{M[i,j]:.0f}°",ha="center",va="center",fontsize=9)
    d3 = res[key]["triples"][0]["Delta3"]
    ax.set_title(f"{title}\nΔ₃={d3}"); ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels([f"q{i+1}" for i in range(4)]); ax.set_yticklabels([f"q{i+1}" for i in range(4)])
fig.colorbar(im, ax=axes[:3], shrink=0.8, label="相対位相の角 α（90°＝|k|=1）")
th = np.linspace(0, np.pi/2, 91)
d3 = []
for t in th:
    Q = ce.make_set(t); u = ce.qmul(ce.qinv(Q[0]),Q[1]); v = ce.qmul(ce.qinv(Q[0]),Q[2])
    c = ce.qmul(u,v)-ce.qmul(v,u); d3.append(0.25*np.dot(c,c))
axes[3].plot(np.degrees(th), d3); axes[3].set_xlabel("m と n のなす角 θ [°]")
axes[3].set_ylabel("Δ₃"); axes[3].set_title("三体判別式 Δ₃（Σq=0, Σq²=0 は全θで成立）")
axes[3].grid(alpha=.3)
for ext in ["svg","png"]:
    fig.savefig(f"condensate_example_20261024.{ext}", dpi=150, bbox_inches="tight")
print("ok", max(d3))
# SVG を軽くする（数値を小数1桁に丸め、メタデータを除く）
import re
s = open("condensate_example_20261024.svg").read()
s = re.sub(r'<metadata>.*?</metadata>', '', s, flags=re.S)
s = re.sub(r'(\d+\.\d)\d+', r'\1', s)
s = re.sub(r'\n\s*', '\n', s)
open("condensate_example_20261024.svg", "w").write(s)
