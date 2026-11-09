# -*- coding: utf-8 -*-
# index_v12q4i.html の初期値を「斜めの衝突・振幅ばらつき」に（2026-11-09 木原氏「正面衝突ではない角度を持った衝突、パワーも 0.5 1 1.5 2 程度に」）。
# 核の主峰が狭い（半値幅 1.2°）ので、軸を少し違えた二波が主峰の中で出会うには出会う場所が原点に近い必要がある（d_min ≈ Θ × 軸の角）。
#  対 1–2：n₁=(1,0,0)、n₂=(1,0.05,0)（軸の角 2.86°）、θ₁=0 v₁=+1 / θ₂=30 v₂=−1 → tick 15 に Θ=15、最接近 d=0.75°、F≈0.75。δ₂=12 で tick 15 に Δ_F≈180°。A₁=0.5、A₂=1.0
#  対 3–4：n₃=(0.5,0.5,0.7)、n₄=(0.5,0.55,0.7)（軸の角 2.43°）、θ₃=5 v₃=+2 / θ₄=35 v₄=−1 → tick 10 に Θ=25、最接近 d≈1.06°、F≈0.55。δ₄=9 で tick 10 に Δ_F≈90°。A₃=1.5、A₄=2.0
# 走行系は一切変えない（スライダーの value と表示だけ）。冪等。
import io, re
p = "index_v12q4i.html"; s = io.open(p, encoding="utf-8").read()
def slider(i, key, val, label):
    global s
    pat = re.compile(r'(<input id="%s%s" type="range"[^>]*?value=")[^"]*("><span class="v" id="%s%sv">)[^<]*(</span>)' % (key, i, key, i))
    s, n = pat.subn(lambda m: m.group(1) + str(val) + m.group(2) + label + m.group(3), s); assert n == 1, (key, i, n)
cfg = {1: dict(th=0, vv=1, dl=0, aa=0.5, n=(1, 0, 0)), 2: dict(th=30, vv=-1, dl=12, aa=1.0, n=(1, 0.05, 0)),
       3: dict(th=5, vv=2, dl=0, aa=1.5, n=(0.5, 0.5, 0.7)), 4: dict(th=35, vv=-1, dl=9, aa=2.0, n=(0.5, 0.55, 0.7))}
for i, c in cfg.items():
    slider(i, "th", c["th"], "%d°" % c["th"]); slider(i, "vv", c["vv"], "%d°/t" % c["vv"]); slider(i, "dl", c["dl"], "%d°/t" % c["dl"]); slider(i, "aa", c["aa"], "%.1f" % c["aa"])
    for key, val in zip(("ni", "nj", "nk"), c["n"]): slider(i, key, val, str(val))
io.open(p, "w", encoding="utf-8").write(s); print("defaults set")
