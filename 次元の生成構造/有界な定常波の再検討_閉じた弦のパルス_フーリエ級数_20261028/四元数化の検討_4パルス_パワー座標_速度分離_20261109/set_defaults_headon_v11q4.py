# -*- coding: utf-8 -*-
# index_v11q4.html の初期値を正面衝突の配置にする（2026-11-09 木原氏「初期値のままで実行すればよいようにして」）。
# 対 1–2：軸 x、θ₁=20 v₁=+1 / θ₂=200 v₂=−1（tick 90、Θ=110° で衝突）。対 3–4：軸 (0.5, 0.5, 0.7)（正規化で (0.503, 0.503, 0.704)）、θ₃=40 v₃=+2 / θ₄=280 v₄=−1（tick 80、Θ=200° で衝突）。
# 11/9 木原氏「二つの組が同じ平面に乗っているのがつまらない、違う面に」。原点を通る二本の直線は必ず一つの平面を張るので、張る平面を座標平面から外す（x 軸と (0.5,0.5,0.7) が張る斜めの面）。
# φ=1、M=100、δ=0、A=1 は共通。右パネルの既定を「座標」にする。走行系（計算・描画）は一切変えない。
import io, re
p = "index_v11q4.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new, n=1):
    global s
    c = s.count(old); assert c == n, (old[:80], c, n); s = s.replace(old, new)
def slider(i, key, val, label):
    global s
    pat = re.compile(r'(<input id="%s%s" type="range"[^>]*?value=")[^"]*("><span class="v" id="%s%sv">)[^<]*(</span>)' % (key, i, key, i))
    s, n = pat.subn(lambda m: m.group(1) + str(val) + m.group(2) + label + m.group(3), s); assert n == 1, (key, i, n)
cfg = {1: dict(th=20, vv=1, dl=0, n=(1, 0, 0)), 2: dict(th=200, vv=-1, dl=0, n=(1, 0, 0)), 3: dict(th=40, vv=2, dl=0, n=(0.5, 0.5, 0.7)), 4: dict(th=280, vv=-1, dl=0, n=(0.5, 0.5, 0.7))}
for i, c in cfg.items():
    slider(i, "th", c["th"], "%d°" % c["th"]); slider(i, "vv", c["vv"], "%d°/t" % c["vv"]); slider(i, "dl", c["dl"], "%d°/t" % c["dl"])
    for key, val in zip(("ni", "nj", "nk"), c["n"]): slider(i, key, val, str(val))
rep('<input type="radio" name="rc" value="0" checked>円周上', '<input type="radio" name="rc" value="0">円周上')
rep('<input type="radio" name="rc" value="1">座標', '<input type="radio" name="rc" value="1" checked>座標')
io.open(p, "w", encoding="utf-8").write(s); print("ok")
