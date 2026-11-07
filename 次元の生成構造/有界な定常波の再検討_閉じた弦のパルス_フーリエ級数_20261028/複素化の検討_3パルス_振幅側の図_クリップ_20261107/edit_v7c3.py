# -*- coding: utf-8 -*-
import io, re
p = "index_v7c3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の複素パルス 3 波（index_v5c3.html の振幅側の図版：右パネルは各パルスを自分の中心位相が原点に来るよう平行移動）</title>",
    "<title>閉じた弦の複素パルス 3 波（index_v6c3.html のクリップ版：右パネルは各パルスを中心位相 ± 位相幅の範囲だけ描く）</title>")
start = s.index("  p.pulses.forEach((u, i) => {\n    const [ddx, ddy] = shifts[i], Zi = Zs[i];\n")
end = s.index("  traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x'")
body = s[start:end]
assert body.count("  });\n") == 1
new_body = (
"  p.pulses.forEach((u, i) => {\n"
"    const [ddx, ddy] = shifts[i], Zi = Zs[i];\n"
"    // クリップ：ψ が θ_i − φ_i 〜 θ_i ＋ φ_i の範囲だけ描く（φ_i は位相幅のスライダー値）。相対角の順に並べ、360° の継ぎ目をまたいでも一本につなぐ\n"
"    const idx = []; for (let k = 0; k < n - 1; k++) if (Math.abs(wrap180(g.psi[k] - u.theta)) <= u.phi) idx.push(k);\n"
"    idx.sort((a, b) => wrap180(g.psi[a] - u.theta) - wrap180(g.psi[b] - u.theta));\n"
"    const C = idx.map(k => g.c[k]), S = idx.map(k => g.s[k]), PSI = idx.map(k => g.psi[k]);\n"
"    { const m = 73, xs = [], ys = [], z0 = [], z1 = [];   // 同じ範囲の部分円筒（高さ ±1）\n"
"      for (let q = 0; q < m; q++) { const a = Math.PI * (u.theta - u.phi + 2 * u.phi * q / (m - 1)) / 180; xs.push(p.Rd * Math.cos(a) + ddx); ys.push(p.Rd * Math.sin(a) + ddy); z0.push(-1); z1.push(1); }\n"
"      traces2.push({type: 'surface', x: [xs, xs], y: [ys, ys], z: [z0, z1], scene: 'scene2', opacity: 0.12, showscale: false, colorscale: [[0, '#8a8a8f'], [1, '#8a8a8f']], hoverinfo: 'skip', name: '円筒', showlegend: false}); }\n"
"    const Xc = C.map(c => p.Rd * c + ddx), Yc = S.map(s => p.Rd * s + ddy);   // 同じ範囲の基準円の弧\n"
"    traces2.push({type: 'scatter3d', mode: 'lines', x: Xc, y: Yc, z: Xc.map(() => 0), scene: 'scene2', line: {color: '#8a8a8f', width: 2}, hoverinfo: 'skip', showlegend: false});\n"
"    const radI = pmode === 1 ? Zi.re.map((v, k) => v * v + Zi.im[k] * Zi.im[k]) : pmode === 2 ? Zi.re.map(v => v * v) : Zi.re;\n"
"    const verI = pmode === 1 ? null : pmode === 2 ? Zi.im.map(v => v * v) : Zi.im;\n"
"    traces2.push({type: 'scatter3d', mode: 'lines', x: C.map((c, j) => (p.Rd + radI[idx[j]]) * c + ddx), y: S.map((s, j) => (p.Rd + radI[idx[j]]) * s + ddy), z: Xc.map(() => 0), scene: 'scene2', line: {color: pmode === 1 ? '#2ca02c' : '#2f6fdf', width: 5},\n"
"      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>` + (pmode === 1 ? '|Z|²' : pmode === 2 ? 'Re²' : 'Re Z') + '=%{customdata[1]:.4f}（半径方向）<extra></extra>', customdata: PSI.map((v, j) => [v, radI[idx[j]]])});\n"
"    if (verI) traces2.push({type: 'scatter3d', mode: 'lines', x: Xc, y: Yc, z: idx.map(k => verI[k]), scene: 'scene2', line: {color: '#e8862a', width: 5},\n"
"      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata:.2f}°<br>` + (pmode === 2 ? 'Im²' : 'Im Z') + '=%{z:.4f}（垂直）<extra></extra>', customdata: PSI});\n"
"    traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: ['diamond', 'square', 'cross'][i], size: 3, color: ['#1f77b4', '#9467bd', '#2ca02c'][i]}, name: '中心 ' + (i + 1), hovertemplate: `パルス ${i + 1} の中心 θ = ${u.theta}°（原点に）<extra></extra>`});\n"
"  });\n")
s = s[:start] + new_body + s[end:]
rep("、各円の中心は反対側の円周上）`;", "、各円の中心は反対側の円周上）　クリップ：各パルス θ_i ± φ_i（${p.pulses.map((u, i) => `${i + 1}: ${u.theta}° ± ${u.phi}°`).join('、')}）`;")
rep("位相 ψ・振幅の向き・縮尺はそのままで、各パルスの位相幅は原点の周りで自分の位相の向きと直交する接線の上に乗る。",
    "位相 ψ・振幅の向き・縮尺はそのままで、各パルスの位相幅は原点の周りで自分の位相の向きと直交する接線の上に乗る。描くのは各パルスの中心位相 ± 位相幅（θ_i − φ_i 〜 θ_i ＋ φ_i）の範囲だけで、基準円と円筒も同じ範囲の弧にクリップする。")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v7c3_check.js", "w", encoding="utf-8").write(js)
