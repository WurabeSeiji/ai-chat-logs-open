# -*- coding: utf-8 -*-
# index_v5q3.html → index_v7q3.html：複素版 index_v7c3 と同じ構成。右パネル＝各パルスを自分の中心位相が原点に来るよう平行移動し、山だけ描く。初期値を変更
import io, re
p = "index_v7q3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の四元数パルス 3 波（index_v4q3.html のパワー表示 3 択版）</title>",
    "<title>閉じた弦の四元数パルス 3 波（index_v5q3.html の振幅側の図版：右パネルは各パルスを自分の中心位相が原点に来るよう平行移動し、山だけ描く）</title>")
rep("let Lfix = null;                              // 3D の枠 ±L。null なら次の draw で合わせ直す（スケール調整ボタン）。それ以外では R が変わっても枠は動かさない\n",
    "let Lfix = null;                              // 3D の枠 ±L。null なら次の draw で合わせ直す（スケール調整ボタン）。それ以外では R が変わっても枠は動かさない\n")
rep("  // ---- 右パネル：観測者の時空。観測者 o を原点に置き、他の二つのパルス i を点 P_i = (Δ_io, Φ_io) で描く。\n",
    "  // ---- 右パネル：複素数（四元数）は全て原点にある。パルス i ごとに、自分の中心位相の点 C_i = kR(cos θ_i, sin θ_i) が原点に来るように −C_i だけ平行移動して、左と同じ描き方で Q_i を描く\n"
    "  //      （各パルスの円は原点を通り、中心は反対側の円周上 −C_i。位相 ψ・振幅の向き・縮尺・表示の 3 択はそのまま。計算は変えない。重心と三角形は和の量なので描かない）\n"
    "  //      描くのはそのパルスの山だけ：θ_i の格子点から両側へ、自分のパワー |Q_i|² が下り続ける間または最大の 5% 以上の間は広げ、5% を下回ったあと最初に上がる根元の谷で止める（谷を含む）\n"
    "  const traces2 = [];\n"
    "  const shifts = p.pulses.map(u => [-p.Rd * Math.cos(Math.PI * u.theta / 180), -p.Rd * Math.sin(Math.PI * u.theta / 180)]);\n"
    "  p.pulses.forEach((u, i) => {\n"
    "    const [ddx, ddy] = shifts[i], Zi = Zs[i], nv = u.nv;\n"
    "    const N = n - 1, at = k => ((k % N) + N) % N, Pi = Zi.re.map((v, k) => v * v + Zi.im[k] * Zi.im[k]);   // |Q_i|² = W_i² + |V_i|²（|n| = 1）\n"
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
    "    // 各パルスの 4 成分：W_i = Re、V_i の i・j・k 成分 = Im × n_i の成分（左の線形合成の各パルスの寄与）\n"
    "    const Wi = Zi.re, Vii = Zi.im.map(v => v * nv[0]), Vji = Zi.im.map(v => v * nv[1]), Vki = Zi.im.map(v => v * nv[2]);\n"
    "    const radial = (A, color, width, label) => ({type: 'scatter3d', mode: 'lines', x: C.map((c, j) => (p.Rd + A[idx[j]]) * c + ddx), y: S.map((s, j) => (p.Rd + A[idx[j]]) * s + ddy), z: Xc.map(() => 0), scene: 'scene2', line: {color, width},\n"
    "      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>${label}=%{customdata[1]:.4f}（半径方向）<extra></extra>`, customdata: PSI.map((v, j) => [v, A[idx[j]]])});\n"
    "    const vertical = (A, color, width, label, sgn) => ({type: 'scatter3d', mode: 'lines', x: Xc, y: Yc, z: idx.map(k => sgn * A[k]), scene: 'scene2', line: {color, width},\n"
    "      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata:.2f}°<br>${label}=%{z:.4f}（垂直${sgn < 0 ? '、反対符号で描画' : ''}）<extra></extra>`, customdata: PSI});\n"
    "    if (pmode === 1) {                                           // パワー合計：|Q_i|² を半径方向に 1 本（緑）\n"
    "      traces2.push(radial(Pi, '#2ca02c', 5, '|Q|²'));\n"
    "    } else if (pmode === 2) {                                    // パワー分離：軸ごと W²・V_j²（半径方向、黒・赤）、V_i²・V_k²（垂直、緑・青）\n"
    "      traces2.push(radial(Wi.map(v => v * v), fg, 5, 'W²'), radial(Vji.map(v => v * v), '#d62728', 4, 'V_j²'), vertical(Vii.map(v => v * v), '#2ca02c', 4, 'V_i²', 1), vertical(Vki.map(v => v * v), '#2f6fdf', 5, 'V_k²', 1));\n"
    "    } else {                                                     // 4 成分：黒 W・赤 V_j（半径方向）、緑 −V_i・青 V_k（垂直）。左と同じ置き方\n"
    "      traces2.push(radial(Wi, fg, 5, 'Re Q'), radial(Vji, '#d62728', 4, 'j 成分'), vertical(Vii, '#2ca02c', 4, 'i 成分', -1), vertical(Vki, '#2f6fdf', 5, 'k 成分', 1));\n"
    "    }\n"
    "    traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: ['diamond', 'square', 'cross'][i], size: 3, color: ['#1f77b4', '#9467bd', '#2ca02c'][i]}, name: '中心 ' + (i + 1), hovertemplate: `パルス ${i + 1} の中心 θ = ${u.theta}°（原点に）<extra></extra>`});\n"
    "  });\n"
    "  traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点<extra></extra>'});\n"
    "  // ---- 旧右パネル（観測者の時空）の計算は読める量（h2）のために残す。曲線は描かない\n"
    "  // ---- 旧：観測者の時空。観測者 o を原点に置き、他の二つのパルス i を点 P_i = (Δ_io, Φ_io) で描く。\n")
rep("${Lnow > L ? '、はみ出し' : ''}）`;",
    "${Lnow > L ? '、はみ出し' : ''}）　右：各パルスの中心位相を原点に（平行移動 ${shifts.map((d, i) => `${i + 1}: (${d[0].toFixed(3)}, ${d[1].toFixed(3)})`).join('、')}、各円の中心は反対側の円周上）　クリップ：各パルスの山（頂点から根元の谷まで、θ_i からの相対角）${p.pulses.map((u, i) => `${i + 1}: ${u.clipA.toFixed(1)}°〜+${u.clipB.toFixed(1)}°`).join('、')}`;")
rep("    scene: {domain: {x: [0, 0.6], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // x・y・z を同じ範囲・同じ縮尺に固定\n",
    "    scene: {domain: {x: [0, 0.5], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // x・y・z を同じ範囲・同じ縮尺に固定\n")
old_axes = re.search(r"    xaxis2: \{domain: \[0\.66, 1\].*?\n    yaxis2: \{.*?\n", s, re.S).group(0)
rep(old_axes,
    "    scene2: {domain: {x: [0.5, 1], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // 右パネル：左と同じ枠・同じ向き・同じ縮尺\n"
    "            xaxis: ax('x', [-L, L]), yaxis: ax('y', [-L, L]), zaxis: ax('z', zrange),\n"
    "            camera: {eye: {x: 1.2, y: -1.45, z: 0.72}, projection: {type: 'orthographic'}}}\n")
rep("  Plotly.react('plot', traces, layout, {responsive: true});",
    "  Plotly.react('plot', traces.filter(t => t.xaxis !== 'x2').concat(traces2), layout, {responsive: true});   // 左の図＋右の図。観測者の時空の曲線は描かない")
rep("「パワー分離」は、4 成分（黒 W・赤 V_j・緑 −V_i・青 V_k）の代わりに、軸ごとのパワー W²・V_j²（半径方向、黒・赤）と V_i²・V_k²（垂直、緑・青）を同じ置き方で描く。",
    "「パワー分離」は、4 成分（黒 W・赤 V_j・緑 −V_i・青 V_k）の代わりに、軸ごとのパワー W²・V_j²（半径方向、黒・赤）と V_i²・V_k²（垂直、緑・青）を同じ置き方で描く。右パネルでは四元数は全て原点にある：パルス i ごとに、円の中心と自分の中心位相の点 C_i ＝ kR(cos θ_i, sin θ_i) を結ぶ線分を半径の分だけ中心の方向へ動かし、C_i を原点に、円の中心を反対側の円周上 −C_i に置いて、左と同じ描き方で Q_i（各パルスの寄与 W_i ＝ Re、V_i ＝ Im·n_i）を描く。位相 ψ・振幅の向き・縮尺はそのままで、各パルスの位相幅は原点の周りで自分の位相の向きと直交する接線の上に乗る。描くのは各パルスの山だけ：中心位相 θ_i から両側へ、そのパルス自身のパワー |Q_i|² が下り続ける間（または最大の 5% 以上の間）は広げ、5% を下回ったあと最初に上がり始める根元の谷で止める（谷を含む）。基準円と円筒も同じ範囲の弧にクリップする。観測者の時空の図は外し、その読みは文字のまま残す。")
# 初期値：θ = 97, 194, 58、φ = 1、M = 100、δ = 20, 0, 10、R× = 1、表示はパワー合計（複素版 index_v7c3 と同じ。A は 1、n は k のまま）
rep('<input id="th1" type="range" min="0" max="360" step="1" value="180"><span class="v" id="th1v">180°</span>', '<input id="th1" type="range" min="0" max="360" step="1" value="97"><span class="v" id="th1v">97°</span>')
rep('<input id="th2" type="range" min="0" max="360" step="1" value="180"><span class="v" id="th2v">180°</span>', '<input id="th2" type="range" min="0" max="360" step="1" value="194"><span class="v" id="th2v">194°</span>')
rep('<input id="th3" type="range" min="0" max="360" step="1" value="180"><span class="v" id="th3v">180°</span>', '<input id="th3" type="range" min="0" max="360" step="1" value="58"><span class="v" id="th3v">58°</span>')
rep('<input id="ph1" type="range" min="0" max="360" step="1" value="20"><span class="v" id="ph1v">20°</span>', '<input id="ph1" type="range" min="0" max="360" step="1" value="1"><span class="v" id="ph1v">1°</span>')
rep('<input id="ph2" type="range" min="0" max="360" step="1" value="10"><span class="v" id="ph2v">10°</span>', '<input id="ph2" type="range" min="0" max="360" step="1" value="1"><span class="v" id="ph2v">1°</span>')
rep('<input id="ph3" type="range" min="0" max="360" step="1" value="10"><span class="v" id="ph3v">10°</span>', '<input id="ph3" type="range" min="0" max="360" step="1" value="1"><span class="v" id="ph3v">1°</span>')
for j in (1, 2, 3): rep(f'<input id="mm{j}" type="range" min="1" max="100" step="1" value="10"><span class="v" id="mm{j}v">10</span>', f'<input id="mm{j}" type="range" min="1" max="100" step="1" value="100"><span class="v" id="mm{j}v">100</span>')
rep('<input id="rr" type="range" min="1" max="20" step="1" value="10"><span class="v" id="rrv">10</span>', '<input id="rr" type="range" min="1" max="20" step="1" value="1"><span class="v" id="rrv">1</span>')
rep('<input type="radio" name="pw" value="0" checked>4 成分</label><label class="ck"><input type="radio" name="pw" value="1">パワー合計', '<input type="radio" name="pw" value="0">4 成分</label><label class="ck"><input type="radio" name="pw" value="1" checked>パワー合計')
rep("初期値は φ ＝ (20, 10, 10)、δ ＝ (20, 0, 10)、観測者 3、標準 23：パルス 1 は傾き 1 で光円錐の上、パルス 2 は縦軸の上（位置は動かず内部時計だけ進む）、辺 1–2 は傾き 2。",
    "初期値は θ ＝ (97, 194, 58)、φ ＝ (1, 1, 1)、M ＝ 100、δ ＝ (20, 0, 10)、R× ＝ 1、表示はパワー合計、A ＝ 1、n ＝ k、観測者 3、標準 23（2026-11-07 木原氏の指定。index_v5q3 までの初期値 φ ＝ (20, 10, 10)、M ＝ 10、R× ＝ 10 は URL 引数で再現できる）。")
rep("（既定は φ₁=20, δ₁=20, φ₂=10, δ₂=0, φ₃=10, δ₃=10、A は全部 1、n は全部 k）", "（既定は θ=97,194,58、φ=1、M=100、δ₁=20, δ₂=0, δ₃=10、k=1、pow=1、A は全部 1、n は全部 k）")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v7q3_check.js", "w", encoding="utf-8").write(js)
