# -*- coding: utf-8 -*-
# index_v12q4i.html の初期値：斜めの衝突・振幅ばらつき・幅広のパルス（2026-11-9 木原氏「角度を持った衝突、パワー 0.5 1 1.5 2、半値幅を広げて衝突しやすく見やすく」）。
# 核の距離を位相の差 |Θ_i − Θ_j| に直したので、出会いは Θ ≈ 110 に置ける（軸の角は向きの因子 (1+n_i·n_j)/2 にだけ入る）。
#  幅 φ = 20°（全パルス）：核の半値幅 10.2°、15° 以遠の裾は 10⁻³ 以下。
#  対 1–2：n₁=(1,0,0)、n₂=(1,0.3,0)（軸の角 16.7°、因子 0.979）、θ₁=20 v₁=+1 / θ₂=200 v₂=−1 → tick 90 に Θ=110。δ₂=2 で Δ_F≈180°。A₁=0.5、A₂=1.0
#  対 3–4：n₃=(0.5,0.5,0.7)、n₄=(0.5,0.8,0.7)（軸の角 12.8°、因子 0.988）、θ₃=40 v₃=+2 / θ₄=280 v₄=−1 → tick 80 に Θ=200。δ₄=1 で Δ_F≈80°。A₃=1.5、A₄=2.0
# 走行系は一切変えない（スライダーの value と表示だけ）。冪等。
import io, re
p = "index_v12q4i.html"; s = io.open(p, encoding="utf-8").read()
def slider(i, key, val, label):
    global s
    pat = re.compile(r'(<input id="%s%s" type="range"[^>]*?value=")[^"]*("><span class="v" id="%s%sv">)[^<]*(</span>)' % (key, i, key, i))
    s, n = pat.subn(lambda m: m.group(1) + str(val) + m.group(2) + label + m.group(3), s); assert n == 1, (key, i, n)
cfg = {1: dict(th=20, vv=1, dl=0, aa=0.5, ph=20, n=(1, 0, 0)), 2: dict(th=200, vv=-1, dl=2, aa=1.0, ph=20, n=(1, 0.3, 0)),
       3: dict(th=40, vv=2, dl=0, aa=1.5, ph=20, n=(0.5, 0.5, 0.7)), 4: dict(th=280, vv=-1, dl=1, aa=2.0, ph=20, n=(0.5, 0.8, 0.7))}
for i, c in cfg.items():
    slider(i, "th", c["th"], "%d°" % c["th"]); slider(i, "vv", c["vv"], "%d°/t" % c["vv"]); slider(i, "dl", c["dl"], "%d°/t" % c["dl"]); slider(i, "aa", c["aa"], "%.1f" % c["aa"]); slider(i, "ph", c["ph"], "%d°" % c["ph"])
    for key, val in zip(("ni", "nj", "nk"), c["n"]): slider(i, key, val, str(val))
io.open(p, "w", encoding="utf-8").write(s); print("defaults set")
