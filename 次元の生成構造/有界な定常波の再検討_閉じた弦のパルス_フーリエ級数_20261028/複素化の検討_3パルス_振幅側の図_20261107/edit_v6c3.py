# -*- coding: utf-8 -*-
import io, re
p = "index_v6c3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の複素パルス 3 波（index_v4c3.html のパワー表示 3 択版）</title>",
    "<title>閉じた弦の複素パルス 3 波（index_v5c3.html の振幅側の図版：右パネルは各パルスを自分の中心位相が原点に来るよう平行移動）</title>")
rep("  // ---- 右パネル：観測者の時空。観測者 o を原点に置き、他の二つのパルス i を点 P_i = (Δ_io, Φ_io) で描く。\n",
    "  // ---- 右パネル：複素数は全て原点にある。パルス i ごとに、自分の中心位相の点 C_i = kR(cos θ_i, sin θ_i) が原点に来るように −C_i だけ平行移動して、左と同じ描き方で Z_i を描く\n"
    "  //      （各パルスの円は原点を通り、中心は反対側の円周上 −C_i。位相 ψ・振幅の向き・縮尺・表示の 3 択はそのまま。計算は変えない。重心と三角形は和の量なので描かない）\n"
    "  const traces2 = [];\n"
    "  const shifts = p.pulses.map(u => [-p.Rd * Math.cos(Math.PI * u.theta / 180), -p.Rd * Math.sin(Math.PI * u.theta / 180)]);\n"
    "  p.pulses.forEach((u, i) => {\n"
    "    const [ddx, ddy] = shifts[i], Zi = Zs[i];\n"
    "    const cyl = cylinder(p.Rd, -1, 1); cyl.x = cyl.x.map(r => r.map(v => v + ddx)); cyl.y = cyl.y.map(r => r.map(v => v + ddy)); cyl.scene = 'scene2'; traces2.push(cyl);\n"
    "    const Xc = X0.map(v => v + ddx), Yc = Y0.map(v => v + ddy);\n"
    "    traces2.push({type: 'scatter3d', mode: 'lines', x: Xc, y: Yc, z: Xc.map(() => 0), scene: 'scene2', line: {color: '#8a8a8f', width: 2}, hoverinfo: 'skip', showlegend: false});\n"
    "    const radI = pmode === 1 ? Zi.re.map((v, k) => v * v + Zi.im[k] * Zi.im[k]) : pmode === 2 ? Zi.re.map(v => v * v) : Zi.re;\n"
    "    const verI = pmode === 1 ? null : pmode === 2 ? Zi.im.map(v => v * v) : Zi.im;\n"
    "    traces2.push({type: 'scatter3d', mode: 'lines', x: Array.from(g.c, (c, k) => (p.Rd + radI[k]) * c + ddx), y: Array.from(g.s, (s, k) => (p.Rd + radI[k]) * s + ddy), z: Xc.map(() => 0), scene: 'scene2', line: {color: pmode === 1 ? '#2ca02c' : '#2f6fdf', width: 5},\n"
    "      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>` + (pmode === 1 ? '|Z|²' : pmode === 2 ? 'Re²' : 'Re Z') + '=%{customdata[1]:.4f}（半径方向）<extra></extra>', customdata: Array.from(g.psi, (v, k) => [v, radI[k]])});\n"
    "    if (verI) traces2.push({type: 'scatter3d', mode: 'lines', x: Xc, y: Yc, z: Array.from(verI), scene: 'scene2', line: {color: '#e8862a', width: 5},\n"
    "      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata:.2f}°<br>` + (pmode === 2 ? 'Im²' : 'Im Z') + '=%{z:.4f}（垂直）<extra></extra>', customdata: Array.from(g.psi)});\n"
    "    traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: ['diamond', 'square', 'cross'][i], size: 3, color: ['#1f77b4', '#9467bd', '#2ca02c'][i]}, name: '中心 ' + (i + 1), hovertemplate: `パルス ${i + 1} の中心 θ = ${u.theta}°（原点に）<extra></extra>`});\n"
    "  });\n"
    "  traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点<extra></extra>'});\n"
    "  // ---- 旧右パネル（観測者の時空）の計算は読める量（h2）のために残す。曲線は描かない\n"
    "  // ---- 旧：観測者の時空。観測者 o を原点に置き、他の二つのパルス i を点 P_i = (Δ_io, Φ_io) で描く。\n")
rep("${Lnow > L ? '、はみ出し' : ''}）`;",
    "${Lnow > L ? '、はみ出し' : ''}）　右：各パルスの中心位相を原点に（平行移動 ${shifts.map((d, i) => `${i + 1}: (${d[0].toFixed(3)}, ${d[1].toFixed(3)})`).join('、')}、各円の中心は反対側の円周上）`;")
rep("    scene: {domain: {x: [0, 0.6], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // x・y・z を同じ範囲・同じ縮尺に固定\n",
    "    scene: {domain: {x: [0, 0.5], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // x・y・z を同じ範囲・同じ縮尺に固定\n")
old_axes = re.search(r"    xaxis2: \{domain: \[0\.66, 1\].*?\n    yaxis2: \{.*?\n", s, re.S).group(0)
rep(old_axes,
    "    scene2: {domain: {x: [0.5, 1], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // 右パネル：左と同じ枠・同じ向き・同じ縮尺\n"
    "            xaxis: ax('x', [-L, L]), yaxis: ax('y', [-L, L]), zaxis: ax('z', zrange),\n"
    "            camera: {eye: {x: 1.2, y: -1.45, z: 0.72}, projection: {type: 'orthographic'}}}\n")
rep("  Plotly.react('plot', traces, layout, {responsive: true});",
    "  Plotly.react('plot', traces.filter(t => t.xaxis !== 'x2').concat(traces2), layout, {responsive: true});   // 左の図＋右の図。観測者の時空の曲線は描かない")
rep("「パワー分離」は軸ごとのパワー Re²（半径方向、青）と Im²（垂直、橙）を Re・Im と同じ置き方で描く（合計が |Z|²）。",
    "「パワー分離」は軸ごとのパワー Re²（半径方向、青）と Im²（垂直、橙）を Re・Im と同じ置き方で描く（合計が |Z|²）。右パネルでは複素数は全て原点にある：パルス i ごとに、円の中心と自分の中心位相の点 C_i ＝ kR(cos θ_i, sin θ_i) を結ぶ線分を半径の分だけ中心の方向へ動かし、C_i を原点に、円の中心を反対側の円周上 −C_i に置いて、左と同じ描き方で Z_i を描く。位相 ψ・振幅の向き・縮尺はそのままで、各パルスの位相幅は原点の周りで自分の位相の向きと直交する接線の上に乗る。観測者の時空の図は外し、その読みは文字のまま残す。")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v6c3_check.js", "w", encoding="utf-8").write(js)
