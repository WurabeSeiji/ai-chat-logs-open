"""軽量 SVG の図化（Drive 保存用）：condensate_example_results_20261024.json を読む"""
import json, numpy as np, importlib.util, sys
spec = importlib.util.spec_from_file_location("ce","condensate_example_20261024.py")
ce = importlib.util.module_from_spec(spec); sys.argv=["x","/tmp/_r.json"]; spec.loader.exec_module(ce)
res = json.load(open("condensate_example_results_20261024.json"))
def col(a):  # 0°青 → 90°灰 → 180°赤
    t = a/180
    if t < .5: r,g,b = 59+(221-59)*2*t, 76+(221-76)*2*t, 192+(221-192)*2*t
    else:      r,g,b = 221+(180-221)*(2*t-1), 221+(4-221)*(2*t-1), 221+(38-221)*(2*t-1)
    return f"rgb({r:.0f},{g:.0f},{b:.0f})"
W,H,cs = 1000, 300, 42
o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" font-family="Noto Sans CJK JP,sans-serif" font-size="12">',
     '<rect width="100%" height="100%" fill="white"/>']
panels = [("C_degenerate_theta0","θ=0°（m∥n：複素数）"),("intermediate_theta45","θ=45°"),("H_condensate_theta90","θ=90°（m⊥n：四元数）")]
for p,(key,title) in enumerate(panels):
    x0, y0 = 30+p*230, 60
    M = np.zeros((4,4))
    for e in res[key]["pairs"]:
        a,b = e["pair"]; M[a-1,b-1]=M[b-1,a-1]=e["alpha_deg"]
    o.append(f'<text x="{x0+2*cs}" y="25" text-anchor="middle" font-size="13">{title}</text>')
    o.append(f'<text x="{x0+2*cs}" y="43" text-anchor="middle">Δ₃ = {res[key]["triples"][0]["Delta3"]}</text>')
    for i in range(4):
        o.append(f'<text x="{x0-6}" y="{y0+i*cs+cs/2+4}" text-anchor="end">q{i+1}</text>')
        o.append(f'<text x="{x0+i*cs+cs/2}" y="{y0+4*cs+16}" text-anchor="middle">q{i+1}</text>')
        for j in range(4):
            fill = "rgb(59,76,192)" if i==j else col(M[i,j])
            o.append(f'<rect x="{x0+j*cs}" y="{y0+i*cs}" width="{cs}" height="{cs}" fill="{fill}" stroke="white"/>')
            if i!=j: o.append(f'<text x="{x0+j*cs+cs/2}" y="{y0+i*cs+cs/2+4}" text-anchor="middle">{M[i,j]:.0f}°</text>')
o.append(f'<text x="375" y="{H-12}" text-anchor="middle">相対位相の角 α（青 0° ／ 灰 90° ＝ |k|=1 ／ 赤 180°）</text>')
# Δ3(θ) 曲線
gx, gy, gw, gh = 760, 60, 210, 170
th = np.linspace(0, np.pi/2, 46); pts=[]
for t in th:
    Q = ce.make_set(t); u = ce.qmul(ce.qinv(Q[0]),Q[1]); v = ce.qmul(ce.qinv(Q[0]),Q[2])
    c = ce.qmul(u,v)-ce.qmul(v,u); d = 0.25*np.dot(c,c)
    pts.append(f"{gx+gw*t/(np.pi/2):.1f},{gy+gh-gh*d/0.5:.1f}")
o.append(f'<rect x="{gx}" y="{gy}" width="{gw}" height="{gh}" fill="none" stroke="black"/>')
o.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="rgb(31,119,180)" stroke-width="2"/>')
o.append(f'<text x="{gx+gw/2}" y="25" text-anchor="middle" font-size="13">三体判別式 Δ₃</text>')
o.append(f'<text x="{gx+gw/2}" y="43" text-anchor="middle">Σq=0, Σq²=0 は全θで成立</text>')
for v,lab in [(0,"0"),(0.25,"0.25"),(0.5,"0.5")]:
    o.append(f'<text x="{gx-5}" y="{gy+gh-gh*v/0.5+4}" text-anchor="end">{lab}</text>')
for d in [0,45,90]:
    o.append(f'<text x="{gx+gw*d/90}" y="{gy+gh+16}" text-anchor="middle">{d}°</text>')
o.append(f'<text x="{gx+gw/2}" y="{gy+gh+34}" text-anchor="middle">m と n のなす角 θ</text>')
o.append('</svg>')
open("condensate_example_20261024.svg","w").write("\n".join(o))
print(len("\n".join(o)))
