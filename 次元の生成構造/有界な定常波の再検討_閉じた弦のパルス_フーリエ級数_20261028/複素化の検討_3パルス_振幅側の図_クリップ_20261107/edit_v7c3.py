# -*- coding: utf-8 -*-
# index_v6c3.html → index_v7c3.html：右パネルの各パルスを、そのパルスの山（|Z_i|² が最大の 5% 以上の範囲）だけ描く
import io, re
p = "index_v7c3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の複素パルス 3 波（index_v5c3.html の振幅側の図版：右パネルは各パルスを自分の中心位相が原点に来るよう平行移動）</title>",
    "<title>閉じた弦の複素パルス 3 波（index_v6c3.html のクリップ版：右パネルは各パルスの山だけ描く）</title>")
start = s.index("  p.pulses.forEach((u, i) => {\n    const [ddx, ddy] = shifts[i], Zi = Zs[i];\n")
end = s.index("  traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x'")
new_body = (
"  p.pulses.forEach((u, i) => {\n"
"    const [ddx, ddy] = shifts[i], Zi = Zs[i];\n"
"    // クリップ：そのパルスの山だけ描く。中心位相 θ_i の格子点から両側へ、自分のパワー |Z_i|² が下り続ける間、または最大の 5% 以上である間は広げ、\n"
"    //          5% を下回ったあと最初に上がり始める点（根元の谷）で止める。谷の点を含むので山は根元（ほぼ 0）まで描かれ、Gibbs の副山（4.3% 以下）は含まない。\n"
"    //          スライダーの φ_i で切ると M の打ち切りで広がった山の根元が切れるので、山の幅はパワーから読む。360° の継ぎ目をまたいでも一本につなぐ\n"
"    const N = n - 1, at = k => ((k % N) + N) % N, Pi = Zi.re.map((v, k) => v * v + Zi.im[k] * Zi.im[k]);\n"
"    const k0 = at(Math.round(u.theta / 360 * N)); let pk = 0; for (let k = 0; k < N; k++) if (Pi[k] > pk) pk = Pi[k];\n"
"    let a = k0, b = k0;\n"
"    while (k0 - a < N / 2 && (Pi[at(a - 1)] <= Pi[at(a)] || Pi[at(a - 1)] >= 0.05 * pk)) a--;\n"
"    while (b - k0 < N / 2 && (Pi[at(b + 1)] <= Pi[at(b)] || Pi[at(b + 1)] >= 0.05 * pk)) b++;\n"
"    const idx = []; for (let k = a; k <= b; k++) idx.push(at(k));\n"
"    u.clipA = wrap180(g.psi[at(a)] - u.theta); u.clipB = wrap180(g.psi[at(b)] - u.theta);   // 山の範囲（θ_i からの相対角、度）\n"
"    const C = idx.map(k => g.c[k]), S = idx.map(k => g.s[k]), PSI = idx.map(k => g.psi[k]);\n"
"    { const m = 73, xs = [], ys = [], z0 = [], z1 = [];   // 同じ範囲の部分円筒（高さ ±1）\n"
"      for (let q = 0; q < m; q++) { const a2 = Math.PI * (u.theta + u.clipA + (u.clipB - u.clipA) * q / (m - 1)) / 180; xs.push(p.Rd * Math.cos(a2) + ddx); ys.push(p.Rd * Math.sin(a2) + ddy); z0.push(-1); z1.push(1); }\n"
"      traces2.push({type: 'surface', x: [xs, xs], y: [ys, ys], z: [z0, z1], scene: 'scene2', opacity: 0.12, showscale: false, colorscale: [[0, '#8a8a8f'], [1, '#8a8a8f']], hoverinfo: 'skip', name: '円筒', showlegend: false}); }\n"
"    const Xc = C.map(c => p.Rd * c + ddx), Yc = S.map(s => p.Rd * s + ddy);   // 同じ範囲の基準円の弧\n"
"    traces2.push({type: 'scatter3d', mode: 'lines', x: Xc, y: Yc, z: Xc.map(() => 0), scene: 'scene2', line: {color: '#8a8a8f', width: 2}, hoverinfo: 'skip', showlegend: false});\n"
"    const radI = pmode === 1 ? Pi : pmode === 2 ? Zi.re.map(v => v * v) : Zi.re;\n"
"    const verI = pmode === 1 ? null : pmode === 2 ? Zi.im.map(v => v * v) : Zi.im;\n"
"    traces2.push({type: 'scatter3d', mode: 'lines', x: C.map((c, j) => (p.Rd + radI[idx[j]]) * c + ddx), y: S.map((s, j) => (p.Rd + radI[idx[j]]) * s + ddy), z: Xc.map(() => 0), scene: 'scene2', line: {color: pmode === 1 ? '#2ca02c' : '#2f6fdf', width: 5},\n"
"      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>` + (pmode === 1 ? '|Z|²' : pmode === 2 ? 'Re²' : 'Re Z') + '=%{customdata[1]:.4f}（半径方向）<extra></extra>', customdata: PSI.map((v, j) => [v, radI[idx[j]]])});\n"
"    if (verI) traces2.push({type: 'scatter3d', mode: 'lines', x: Xc, y: Yc, z: idx.map(k => verI[k]), scene: 'scene2', line: {color: '#e8862a', width: 5},\n"
"      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata:.2f}°<br>` + (pmode === 2 ? 'Im²' : 'Im Z') + '=%{z:.4f}（垂直）<extra></extra>', customdata: PSI});\n"
"    traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: ['diamond', 'square', 'cross'][i], size: 3, color: ['#1f77b4', '#9467bd', '#2ca02c'][i]}, name: '中心 ' + (i + 1), hovertemplate: `パルス ${i + 1} の中心 θ = ${u.theta}°（原点に）<extra></extra>`});\n"
"  });\n")
s = s[:start] + new_body + s[end:]
rep("、各円の中心は反対側の円周上）`;",
    "、各円の中心は反対側の円周上）　クリップ：各パルスの山（頂点から根元の谷まで、θ_i からの相対角）${p.pulses.map((u, i) => `${i + 1}: ${u.clipA.toFixed(1)}°〜+${u.clipB.toFixed(1)}°`).join('、')}`;")
rep("位相 ψ・振幅の向き・縮尺はそのままで、各パルスの位相幅は原点の周りで自分の位相の向きと直交する接線の上に乗る。",
    "位相 ψ・振幅の向き・縮尺はそのままで、各パルスの位相幅は原点の周りで自分の位相の向きと直交する接線の上に乗る。描くのは各パルスの山だけ：中心位相 θ_i から両側へ、そのパルス自身のパワー |Z_i|² が下り続ける間（または最大の 5% 以上の間）は広げ、5% を下回ったあと最初に上がり始める根元の谷で止める（谷を含む）。基準円と円筒も同じ範囲の弧にクリップする。スライダーの φ_i で切ると M の打ち切りで広がった山の根元が切れるので、幅はパワーから読む。")
rep('<input id="th1" type="range" min="0" max="360" step="1" value="180"><span class="v" id="th1v">180°</span>', '<input id="th1" type="range" min="0" max="360" step="1" value="97"><span class="v" id="th1v">97°</span>')
rep('<input id="th2" type="range" min="0" max="360" step="1" value="180"><span class="v" id="th2v">180°</span>', '<input id="th2" type="range" min="0" max="360" step="1" value="194"><span class="v" id="th2v">194°</span>')
rep('<input id="th3" type="range" min="0" max="360" step="1" value="180"><span class="v" id="th3v">180°</span>', '<input id="th3" type="range" min="0" max="360" step="1" value="58"><span class="v" id="th3v">58°</span>')
rep('<input id="ph1" type="range" min="0" max="360" step="1" value="20"><span class="v" id="ph1v">20°</span>', '<input id="ph1" type="range" min="0" max="360" step="1" value="1"><span class="v" id="ph1v">1°</span>')
rep('<input id="ph2" type="range" min="0" max="360" step="1" value="10"><span class="v" id="ph2v">10°</span>', '<input id="ph2" type="range" min="0" max="360" step="1" value="1"><span class="v" id="ph2v">1°</span>')
rep('<input id="ph3" type="range" min="0" max="360" step="1" value="10"><span class="v" id="ph3v">10°</span>', '<input id="ph3" type="range" min="0" max="360" step="1" value="1"><span class="v" id="ph3v">1°</span>')
for j in (1, 2, 3): rep(f'<input id="mm{j}" type="range" min="1" max="100" step="1" value="10"><span class="v" id="mm{j}v">10</span>', f'<input id="mm{j}" type="range" min="1" max="100" step="1" value="100"><span class="v" id="mm{j}v">100</span>')
rep('<input id="rr" type="range" min="1" max="20" step="1" value="10"><span class="v" id="rrv">10</span>', '<input id="rr" type="range" min="1" max="20" step="1" value="1"><span class="v" id="rrv">1</span>')
rep('<input type="radio" name="pw" value="0" checked>Re・Im</label><label class="ck"><input type="radio" name="pw" value="1">パワー合計', '<input type="radio" name="pw" value="0">Re・Im</label><label class="ck"><input type="radio" name="pw" value="1" checked>パワー合計')
rep("初期値は φ ＝ (20, 10, 10)、δ ＝ (20, 0, 10)、観測者 3、標準 23：パルス 1 は傾き 1 で光円錐の上、パルス 2 は縦軸の上（位置は動かず内部時計だけ進む）、辺 1–2 は傾き 2。",
    "初期値は θ ＝ (97, 194, 58)、φ ＝ (1, 1, 1)、M ＝ 100、δ ＝ (20, 0, 10)、R× ＝ 1、表示はパワー合計、観測者 3、標準 23（2026-11-07 木原氏の指定。index_v6c3 までの初期値 φ ＝ (20, 10, 10)、M ＝ 10、R× ＝ 10 は URL 引数で再現できる）。")
rep("（既定は φ₁=20, δ₁=20, φ₂=10, δ₂=0, φ₃=10, δ₃=10）", "（既定は θ=97,194,58、φ=1、M=100、δ₁=20, δ₂=0, δ₃=10、k=1、pow=1）")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v7c3_check.js", "w", encoding="utf-8").write(js)
