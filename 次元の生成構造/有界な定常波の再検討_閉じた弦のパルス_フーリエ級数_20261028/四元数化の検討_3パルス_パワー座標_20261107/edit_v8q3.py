# -*- coding: utf-8 -*-
# index_v7q3.html → index_v8q3.html：右パネルのパワー表示に「円周上／座標」の切替。座標のときは右に軌道だけ描く：
# 波 i の点 (i_i Θ₁, j_i Θ₂, k_i Θ₃)、Θ = (Θ₁, Θ₂, Θ₃) は三つの中心位相の累積（度、折り返さない）、(i_i, j_i, k_i) = n_i。Auto のティックごとの点をつないだ曲線。原点と結ぶ線分は描かない。左は不変。右の枠は初期 ±360°、Auto で軌道に合わせて計算し直す
import io, re
p = "index_v8q3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の四元数パルス 3 波（index_v5q3.html の振幅側の図版：右パネルは各パルスを自分の中心位相が原点に来るよう平行移動し、山だけ描く）</title>",
    "<title>閉じた弦の四元数パルス 3 波（index_v7q3.html のパワー座標版：右パネルに (i_i Θ₁, j_i Θ₂, k_i Θ₃) の軌道、Θ ＝ 中心位相の累積）</title>")
rep('<label class="ck"><input type="radio" name="pw" value="2">パワー分離 軸ごと W²・V_j²（半径方向）・V_i²・V_k²（垂直）</label></span>\n',
    '<label class="ck"><input type="radio" name="pw" value="2">パワー分離 軸ごと W²・V_j²（半径方向）・V_i²・V_k²（垂直）</label></span>\n'
    '  <span>右のパワー</span>\n'
    '  <span class="grp"><label class="ck"><input type="radio" name="rc" value="0" checked>円周上</label><label class="ck"><input type="radio" name="rc" value="1">座標 波 i ごとに x ＝ i_i Θ₁, y ＝ j_i Θ₂, z ＝ k_i Θ₃（Θ ＝ (Θ₁, Θ₂, Θ₃) は三つの中心位相の累積、度、360° で折り返さない。(i_i, j_i, k_i) ＝ n_i）。Auto のティックごとの点をつないだ軌道、波ごとの色。右の枠は初期 ±360°、Auto で軌道に合わせて計算し直す</label></span>\n')
rep("  const pmode = Number(document.querySelector('input[name=pw]:checked').value), usePow = pmode > 0;   // 表示：0 = 4 成分、1 = パワー合計 |Q|²、2 = パワー分離 軸ごと（計算は変えない）\n",
    "  const pmode = Number(document.querySelector('input[name=pw]:checked').value), usePow = pmode > 0;   // 表示：0 = 4 成分、1 = パワー合計 |Q|²、2 = パワー分離 軸ごと（計算は変えない）\n"
    "  const rcoord = Number(document.querySelector('input[name=rc]:checked').value) === 1;              // 右のパワー：false = 円周上（従来）、true = 座標（軌道）。左には影響しない\n")
rep("let Lfix = null;                              // 3D の枠 ±L。null なら次の draw で合わせ直す（スケール調整ボタン）。それ以外では R が変わっても枠は動かさない\n",
    "let Lfix = null;                              // 3D の枠 ±L。null なら次の draw で合わせ直す（スケール調整ボタン）。それ以外では R が変わっても枠は動かさない\n"
    "let orb = [[], [], []], orbTick = -1;         // 右の座標表示の軌道：波ごとに、ティックごとの点 (i_i Θ₁, j_i Θ₂, k_i Θ₃)。スライダーを動かすと消える\n")
rep("  const traces2 = [];\n",
    "  const traces2 = [];\n"
    "  // 右の座標表示の軌道：波 i の点 = (i_i Θ₁, j_i Θ₂, k_i Θ₃)。Θ = (Θ₁, Θ₂, Θ₃) は三つの中心位相の累積（度、360 で折り返さない。base + cum）、(i_i, j_i, k_i) = n_i。点 (Θ₁, Θ₂, Θ₃) は 1 ティックに (φ₁, φ₂, φ₃) 進む（光速）。ティックごとに 1 点、波ごとにつないで曲線\n  const Th = [0, 1, 2].map(i => base[i][0] + cum[i][cum[i].length - 1][0]);\n"
    "  const orbPts = p.pulses.map(u => { const nv = u.nv; return [nv[0] * Th[0], nv[1] * Th[1], nv[2] * Th[2]]; });\n"
    "  if (orb[0].length === 0 || ticks !== orbTick) { orbPts.forEach((q, i) => { orb[i].push(q); if (orb[i].length > 4000) orb[i].shift(); }); orbTick = ticks; } else orbPts.forEach((q, i) => { orb[i][orb[i].length - 1] = q; });\n  const L2 = Math.max(360, ...orb.map(o => o.reduce((m, q) => Math.max(m, Math.abs(q[0]), Math.abs(q[1]), Math.abs(q[2])), 0)));   // 座標表示の右の枠（度）：初期 ±360、軌道がはみ出したら合わせ直す\n")
rep("  p.pulses.forEach((u, i) => {\n    const [ddx, ddy] = shifts[i], Zi = Zs[i]",
    "  if (rcoord) orb.forEach((o, i) => traces2.push({type: 'scatter3d', mode: 'lines', x: o.map(q => q[0]), y: o.map(q => q[1]), z: o.map(q => q[2]), scene: 'scene2', line: {color: ['#1f77b4', '#9467bd', '#2ca02c'][i], width: 4}, showlegend: false,   // 座標表示：右は波ごとの軌道だけ（原点と結ぶ線分・円筒・弧・印は描かない）\n"
    "    hovertemplate: `波 ${i + 1}：(i_i Θ₁, j_i Θ₂, k_i Θ₃) = (%{x:.1f}, %{y:.1f}, %{z:.1f})<extra></extra>`}));\n"
    "  if (rcoord) orb.forEach((o, i) => {                           // 現在位置に、その波の色の球。半径 = 振幅 A_i（q = A(a + n b) の A）。k は掛けない\n"
    "    const c = o[o.length - 1], r = p.pulses[i].A, col = ['#1f77b4', '#9467bd', '#2ca02c'][i], xs = [], ys = [], zs = [];\n"
    "    for (let a = 0; a <= 12; a++) { const th = Math.PI * a / 12, rx = [], ry = [], rz = []; for (let b = 0; b <= 24; b++) { const ph = 2 * Math.PI * b / 24; rx.push(c[0] + r * Math.sin(th) * Math.cos(ph)); ry.push(c[1] + r * Math.sin(th) * Math.sin(ph)); rz.push(c[2] + r * Math.cos(th)); } xs.push(rx); ys.push(ry); zs.push(rz); }\n"
    "    traces2.push({type: 'surface', x: xs, y: ys, z: zs, scene: 'scene2', opacity: 0.35, showscale: false, colorscale: [[0, col], [1, col]], hoverinfo: 'skip', showlegend: false, name: `波 ${i + 1} の球（半径 A = ${r}）`});\n"
    "  });\n"
    "  p.pulses.forEach((u, i) => {\n    const [ddx, ddy] = shifts[i], Zi = Zs[i]")
rep("    u.clipA = wrap180(g.psi[at(a)] - u.theta); u.clipB = wrap180(g.psi[at(b)] - u.theta);   // 山の範囲（θ_i からの相対角、度）\n",
    "    u.clipA = wrap180(g.psi[at(a)] - u.theta); u.clipB = wrap180(g.psi[at(b)] - u.theta);   // 山の範囲（θ_i からの相対角、度）\n"
    "    if (rcoord) return;                                          // 座標表示では右の各パルス（円筒・弧・曲線・印）は描かない。山の範囲は読みのために計算する\n")
rep("  traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点<extra></extra>'});\n",
    "  if (!rcoord) traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点<extra></extra>'});\n")
rep("${Lnow > L ? '、はみ出し' : ''}）　右：各パルスの中心位相を原点に",
    "${Lnow > L ? '、はみ出し' : ''}）${rcoord ? `　軌道（i_i Θ₁, j_i Θ₂, k_i Θ₃）：Θ = (${Th.map(v => v.toFixed(1)).join(', ')})、` + orbPts.map((q, i) => `${i + 1}: (${q[0].toFixed(1)}, ${q[1].toFixed(1)}, ${q[2].toFixed(1)})`).join('、') + `、各 ${orb[0].length} 点、右の枠 ±${L2.toFixed(0)}°` : ''}　右：各パルスの中心位相を原点に")
rep("for (const r of document.querySelectorAll('input[name=pw]')) r.addEventListener('change', draw);",
    "for (const r of document.querySelectorAll('input[name=pw], input[name=rc]')) r.addEventListener('change', draw);")
rep("const resetAll = () => { ticks = 0; hist = PAIRS.map(() => [[0, 0]]);", "const resetAll = () => { ticks = 0; orb = [[], [], []]; orbTick = -1; hist = PAIRS.map(() => [[0, 0]]);")
rep("  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n",
    "  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n"
    "  if (q.get('rc') === '1') document.querySelector('input[name=rc][value=\"1\"]').checked = true;   // rc=1 で右を座標表示（撮影用）\n")
rep("観測者の時空の図は外し、その読みは文字のまま残す。",
    "観測者の時空の図は外し、その読みは文字のまま残す。「右のパワー」を「座標」にすると、右パネルには波ごとの軌道だけを描く：三つの中心位相の累積 Θ ＝ (Θ₁, Θ₂, Θ₃)（度、360° で折り返さない）を、波 i の軸 n_i ＝ (i_i, j_i, k_i) で成分ごとに写した点 (i_i Θ₁, j_i Θ₂, k_i Θ₃)。点 (Θ₁, Θ₂, Θ₃) は 1 ティックに (φ₁, φ₂, φ₃) だけ進む（光速）。Auto のティックごとに 1 点を取り、波ごとに色を分けてつなぐ。現在位置には、その波の色の球を置く（半径 ＝ 振幅 A_i、q ＝ A(a ＋ n b) の A。表示倍率 k は掛けない）。原点と結ぶ線分・円筒・弧・印は描かない。スライダーを動かすと軌道は消える。右の枠は初期 ±360°、Auto で軌道に合わせて計算し直す。左パネル、読める量、左の枠は変わらない。")
# 追加（11/7 木原氏）：R× の横にスライダー s（1〜100、初期 10）。球の半径に s を掛ける。他は不変
rep('<div class="sl"><span>R×</span><input id="rr" type="range" min="1" max="20" step="1" value="1"><span class="v" id="rrv">1</span></div>\n',
    '<div class="sl"><span>R×</span><input id="rr" type="range" min="1" max="20" step="1" value="1"><span class="v" id="rrv">1</span></div>\n'
    '  <div class="sl"><span>s</span><input id="ss" type="range" min="1" max="100" step="1" value="10"><span class="v" id="ssv">10</span></div>\n')
rep("const rr = $('rr'), btn = $('auto');", "const rr = $('rr'), ss = $('ss'), btn = $('auto');   // ss = 球の倍率 s（1〜100、初期 10）")
rep("  $('rrv').textContent = p.k;\n", "  $('rrv').textContent = p.k; $('ssv').textContent = ss.value;\n")
rep("// 現在位置に、その波の色の球。半径 = 振幅 A_i（q = A(a + n b) の A）。k は掛けない", "// 現在位置に、その波の色の球。半径 = s × 振幅 A_i（q = A(a + n b) の A、s = スライダー s、1〜100）。k は掛けない")
rep('const c = o[o.length - 1], r = p.pulses[i].A, col = ', 'const c = o[o.length - 1], r = Number(ss.value) * p.pulses[i].A, col = ')
rep("name: `波 ${i + 1} の球（半径 A = ${r}）`", "name: `波 ${i + 1} の球（半径 s A = ${r.toFixed(3)}）`")
rep("rr.addEventListener('input', draw);\n", "rr.addEventListener('input', draw);\nss.addEventListener('input', draw);   // s は球の大きさだけ（軌道は消さない）\n")
rep("  if (q.has('k')) rr.value = q.get('k');\n", "  if (q.has('k')) rr.value = q.get('k');\n  if (q.has('s')) ss.value = q.get('s');   // s = 球の倍率（撮影用）\n")
rep("現在位置には、その波の色の球を置く（半径 ＝ 振幅 A_i、q ＝ A(a ＋ n b) の A。表示倍率 k は掛けない）。",
    "現在位置には、その波の色の球を置く（半径 ＝ s × 振幅 A_i、q ＝ A(a ＋ n b) の A、s ＝ R× の横のスライダー s、1〜100、初期 10。表示倍率 k は掛けない）。")
# 追加（11/7 木原氏）：座標表示では右の枠を度にする。初期 ±360、Auto で軌道に合わせて計算し直す（L2）。左の枠 ±L は不変
rep("aspectratio: {x: 1, y: 1, z: 1},   // 右パネル：左と同じ枠・同じ向き・同じ縮尺\n            xaxis: ax('x', [-L, L]), yaxis: ax('y', [-L, L]), zaxis: ax('z', zrange),\n",
    "aspectratio: {x: 1, y: 1, z: 1},   // 右パネル：左と同じ枠・同じ向き・同じ縮尺。座標表示では枠 ±L2（度、初期 ±360、Auto で軌道に合わせて計算し直す）\n            xaxis: ax('x', rcoord ? [-L2, L2] : [-L, L]), yaxis: ax('y', rcoord ? [-L2, L2] : [-L, L]), zaxis: ax('z', rcoord ? [-L2, L2] : zrange),\n")
# 追加（11/7 木原氏）：n の初期値を i, j, k（(1,0,0)・(0,1,0)・(0,0,1)）にする。それまでは全部 k
rep('<input id="ni1" type="range" min="-1" max="1" step="0.05" value="0"><span class="v" id="ni1v">0</span>', '<input id="ni1" type="range" min="-1" max="1" step="0.05" value="1"><span class="v" id="ni1v">1</span>')
rep('<input id="nk1" type="range" min="-1" max="1" step="0.05" value="1"><span class="v" id="nk1v">1</span>', '<input id="nk1" type="range" min="-1" max="1" step="0.05" value="0"><span class="v" id="nk1v">0</span>')
rep('<input id="nj2" type="range" min="-1" max="1" step="0.05" value="0"><span class="v" id="nj2v">0</span>', '<input id="nj2" type="range" min="-1" max="1" step="0.05" value="1"><span class="v" id="nj2v">1</span>')
rep('<input id="nk2" type="range" min="-1" max="1" step="0.05" value="1"><span class="v" id="nk2v">1</span>', '<input id="nk2" type="range" min="-1" max="1" step="0.05" value="0"><span class="v" id="nk2v">0</span>')
rep('A ＝ 1、n ＝ k、観測者 3', 'A ＝ 1、n ＝ (1, 0, 0)・(0, 1, 0)・(0, 0, 1)（i, j, k。index_v8q3 で変更、index_v7q3 までは全部 k）、観測者 3')
rep('A は全部 1、n は全部 k）', 'A は全部 1、n は i, j, k）')
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/vx_check.js", "w", encoding="utf-8").write(js)
