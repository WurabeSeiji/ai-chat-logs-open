# -*- coding: utf-8 -*-
# index_v10q4.html → index_v11q4.html：位置の速度 v_i（°/tick、負も可）を幅 φ_i から切り離す。
# 2026-11-09 木原氏「幅を兼ねているための制限で、設定項目が分かれればレンジの制限は無くなる。物理ではない制限」。
# φ_i は幅だけ（c_m = sin(mπφ/360)/m）。v_i は θ ← θ + v と累積 Θ（右パネルの位置）と対の履歴に入る。既定 v = 1 で v10q4 と同じ動き。
# δ_i（内部時計）・読み・左パネル・球の半径 2/√M は変更なし。
import io, re
p = "index_v11q4.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new, n=1):
    global s
    c = s.count(old); assert c == n, (old[:80], c, n); s = s.replace(old, new)
# スライダー v_i（φ_i の隣）
for i, sub in zip("1234", "₁₂₃₄"):
    rep('<span class="v" id="ph%sv">1°</span></div>\n' % i,
        '<span class="v" id="ph%sv">1°</span></div>\n  <div class="sl"><span>v%s</span><input id="vv%s" type="range" min="-360" max="360" step="1" value="1"><span class="v" id="vv%sv">1°/t</span></div>\n' % (i, sub, i, i))
# P
rep("ph: $('ph' + i), mm: $('mm' + i),", "ph: $('ph' + i), vv: $('vv' + i), mm: $('mm' + i),")
rep("phv: $('ph' + i + 'v'), mmv: $('mm' + i + 'v'),", "phv: $('ph' + i + 'v'), vvv: $('vv' + i + 'v'), mmv: $('mm' + i + 'v'),")
# パルス生成
rep("phi = Number(u.ph.value), M = Number(u.mm.value), delta = Number(u.dl.value), phase0 = u.phase0;",
    "phi = Number(u.ph.value), v = Number(u.vv.value), M = Number(u.mm.value), delta = Number(u.dl.value), phase0 = u.phase0;   // φ は幅、v は速度（°/tick、負も可）")
rep("return {theta, phi, M, delta, phase0, co, A,", "return {theta, phi, v, M, delta, phase0, co, A,")
# 数値表示
rep("u.phv.textContent = q.phi + '°'; u.mmv.textContent = q.M;", "u.phv.textContent = q.phi + '°'; u.vvv.textContent = q.v + '°/t'; u.mmv.textContent = q.M;")
rep("φ = ${u.phi}°, δ = ${u.delta}°/t,", "φ = ${u.phi}°, v = ${u.v}°/t, δ = ${u.delta}°/t,")
# 位置の進み・対の履歴・累積 Θ
rep("let t = Number(u.th.value) + Number(u.ph.value); if (t > 360) t -= 360; u.th.value = String(t);   // 位置 θ ← θ + φ",
    "let t = Number(u.th.value) + Number(u.vv.value); t = ((t % 360) + 360) % 360; u.th.value = String(t);   // 位置 θ ← θ + v（負も折り返す）")
rep("hist[k].push([last[0] + Number(P[i].ph.value) - Number(P[j].ph.value),", "hist[k].push([last[0] + Number(P[i].vv.value) - Number(P[j].vv.value),")
rep("cum[i].push([last[0] + Number(u.ph.value), last[1] + Number(u.dl.value)]);", "cum[i].push([last[0] + Number(u.vv.value), last[1] + Number(u.dl.value)]);")
# resetAll の監視・URL
rep("for (const s of [u.th, u.ph, u.mm, u.dl, u.aa, u.ni, u.nj, u.nk])", "for (const s of [u.th, u.ph, u.vv, u.mm, u.dl, u.aa, u.ni, u.nj, u.nk])")
rep("if (q.has('phi' + j)) u.ph.value = q.get('phi' + j);", "if (q.has('phi' + j)) u.ph.value = q.get('phi' + j); if (q.has('v' + j)) u.vv.value = q.get('v' + j);")
assert s.count("u.ph.value") == 2, s.count("u.ph.value")   # 残るのは幅（180 行）と URL（447 行）だけ
io.open(p, "w", encoding="utf-8").write(s); print("ok")
