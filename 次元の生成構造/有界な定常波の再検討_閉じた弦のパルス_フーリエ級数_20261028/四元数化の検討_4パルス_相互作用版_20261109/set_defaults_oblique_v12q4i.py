# -*- coding: utf-8 -*-
# index_v12q4i.html の初期値（2026-11-09 木原氏「2 体は X 軸、2 体は Y 軸に直交させ、わずかに Z だけ変えて必ずオフセットして衝突するように」）。
#  幅 φ = 20°（全パルス、核の半値幅 10.2°、裾なし）。振幅 0.5 1 1.5 2。
#  対 1–2（X 軸）：n₁=(1,0,0)、n₂=(1,0,0.05)、θ₁=20 v₁=+1 / θ₂=200 v₂=−1 → tick 90 に Θ=110 で位相一致、図のオフセット 5.5°、向きの因子 0.999。δ₂=2 で Δ_F=180°
#  対 3–4（Y 軸）：n₃=(0,1,0)、n₄=(0,1,0.05)、θ₃=220 v₃=+1 / θ₄=340 v₄=−1 → tick 60 に Θ=280 で位相一致、図のオフセット 14.0°、向きの因子 0.999。δ₄=2 で Δ_F=120°
#  組をまたぐ対（因子 0.5）は 150 tick の中で位相が一致しない（1–4 は tick 160、2–3 は負、1–3 と 2–4 は同じ速さ）。球（半径 11.5°）は中心距離 22.9° 未満で重なる。
# 走行系は一切変えない（スライダーの value と表示だけ）。冪等。
import io, re
p = "index_v12q4i.html"; s = io.open(p, encoding="utf-8").read()
def slider(i, key, val, label):
    global s
    pat = re.compile(r'(<input id="%s%s" type="range"[^>]*?value=")[^"]*("><span class="v" id="%s%sv">)[^<]*(</span>)' % (key, i, key, i))
    s, n = pat.subn(lambda m: m.group(1) + str(val) + m.group(2) + label + m.group(3), s); assert n == 1, (key, i, n)
cfg = {1: dict(th=20, vv=1, dl=0, aa=0.5, ph=20, n=(1, 0, 0)), 2: dict(th=200, vv=-1, dl=2, aa=1.0, ph=20, n=(1, 0, 0.05)),
       3: dict(th=220, vv=1, dl=0, aa=1.5, ph=20, n=(0, 1, 0)), 4: dict(th=340, vv=-1, dl=2, aa=2.0, ph=20, n=(0, 1, 0.05))}
for i, c in cfg.items():
    slider(i, "th", c["th"], "%d°" % c["th"]); slider(i, "vv", c["vv"], "%d°/t" % c["vv"]); slider(i, "dl", c["dl"], "%d°/t" % c["dl"]); slider(i, "aa", c["aa"], "%.1f" % c["aa"]); slider(i, "ph", c["ph"], "%d°" % c["ph"])
    for key, val in zip(("ni", "nj", "nk"), c["n"]): slider(i, key, val, str(val))
io.open(p, "w", encoding="utf-8").write(s); print("defaults set")
