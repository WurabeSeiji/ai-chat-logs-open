# -*- coding: utf-8 -*-
# index_v7c3.html → index_v8c3.html：右パネルのパワー表示に「円周上／座標」の切替。座標のときは右に軌道だけ描く：
# 点 (P sin θ₁, P sin θ₂, P sin θ₃)、P = kR（二乗パワー R² の平方根に表示倍率 k）、θ_i は三つのパルスの中心位相。Auto のティックごとの点をつないだ曲線。原点と結ぶ線分は描かない。左は不変
import io, re
p = "index_v8c3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の複素パルス 3 波（index_v6c3.html のクリップ版：右パネルは各パルスの山だけ描く）</title>",
    "<title>閉じた弦の複素パルス 3 波（index_v7c3.html のパワー座標版：右パネルに (P sin θ₁, P sin θ₂, P sin θ₃) の軌道）</title>")
rep('<label class="ck"><input type="radio" name="pw" value="2">パワー分離 Re²（半径方向）・Im²（垂直）</label></span>\n',
    '<label class="ck"><input type="radio" name="pw" value="2">パワー分離 Re²（半径方向）・Im²（垂直）</label></span>\n'
    '  <span>右のパワー</span>\n'
    '  <span class="grp"><label class="ck"><input type="radio" name="rc" value="0" checked>円周上</label><label class="ck"><input type="radio" name="rc" value="1">座標 x ＝ P sin θ₁, y ＝ P sin θ₂, z ＝ P sin θ₃（P ＝ kR。Auto のティックごとの点をつないだ軌道）</label></span>\n')
rep("  const pmode = Number(document.querySelector('input[name=pw]:checked').value), usePow = pmode > 0;   // 表示：0 = Re・Im、1 = パワー合計 |Z|²、2 = パワー分離 Re²・Im²（計算は変えない）\n",
    "  const pmode = Number(document.querySelector('input[name=pw]:checked').value), usePow = pmode > 0;   // 表示：0 = Re・Im、1 = パワー合計 |Z|²、2 = パワー分離 Re²・Im²（計算は変えない）\n"
    "  const rcoord = Number(document.querySelector('input[name=rc]:checked').value) === 1;              // 右のパワー：false = 円周上（従来）、true = 座標（軌道）。左には影響しない\n")
rep("let Lfix = null;                              // 3D の枠 ±L。null なら次の draw で合わせ直す（スケール調整ボタン）。それ以外では R が変わっても枠は動かさない\n",
    "let Lfix = null;                              // 3D の枠 ±L。null なら次の draw で合わせ直す（スケール調整ボタン）。それ以外では R が変わっても枠は動かさない\n"
    "let orb = [], orbTick = -1;                   // 右の座標表示の軌道：ティックごとの点 (kR sin θ₁, kR sin θ₂, kR sin θ₃)。スライダーを動かすと消える\n")
rep("  const traces2 = [];\n",
    "  const traces2 = [];\n"
    "  // 右の座標表示の軌道：点 (P sin θ₁, P sin θ₂, P sin θ₃)、P = kR（二乗パワー R² の平方根 × 表示倍率）、θ_i = 三つのパルスの中心位相。ティックごとに 1 点、つないで曲線\n"
    "  const orbPt = [p.Rd * Math.sin(Math.PI * p.pulses[0].theta / 180), p.Rd * Math.sin(Math.PI * p.pulses[1].theta / 180), p.Rd * Math.sin(Math.PI * p.pulses[2].theta / 180)];\n"
    "  if (orb.length === 0 || ticks !== orbTick) { orb.push(orbPt); orbTick = ticks; if (orb.length > 4000) orb.shift(); } else orb[orb.length - 1] = orbPt;\n")
rep("  p.pulses.forEach((u, i) => {\n    const [ddx, ddy] = shifts[i], Zi = Zs[i]",
    "  if (rcoord) traces2.push({type: 'scatter3d', mode: 'lines', x: orb.map(q => q[0]), y: orb.map(q => q[1]), z: orb.map(q => q[2]), scene: 'scene2', line: {color: '#2ca02c', width: 4}, showlegend: false,   // 座標表示：右は軌道だけ（原点と結ぶ線分・円筒・弧・印は描かない）\n"
    "    hovertemplate: '(kR sin θ₁, kR sin θ₂, kR sin θ₃) = (%{x:.4f}, %{y:.4f}, %{z:.4f})<extra></extra>'});\n"
    "  p.pulses.forEach((u, i) => {\n    const [ddx, ddy] = shifts[i], Zi = Zs[i]")
rep("    u.clipA = wrap180(g.psi[at(a)] - u.theta); u.clipB = wrap180(g.psi[at(b)] - u.theta);   // 山の範囲（θ_i からの相対角、度）\n",
    "    u.clipA = wrap180(g.psi[at(a)] - u.theta); u.clipB = wrap180(g.psi[at(b)] - u.theta);   // 山の範囲（θ_i からの相対角、度）\n"
    "    if (rcoord) return;                                          // 座標表示では右の各パルス（円筒・弧・曲線・印）は描かない。山の範囲は読みのために計算する\n")
rep("  traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点<extra></extra>'});\n",
    "  if (!rcoord) traces2.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点<extra></extra>'});\n")
rep("${Lnow > L ? '、はみ出し' : ''}）　右：各パルスの中心位相を原点に",
    "${Lnow > L ? '、はみ出し' : ''}）${rcoord ? `　軌道：(kR sin θ₁, kR sin θ₂, kR sin θ₃) = (${orbPt[0].toFixed(4)}, ${orbPt[1].toFixed(4)}, ${orbPt[2].toFixed(4)})、${orb.length} 点` : ''}　右：各パルスの中心位相を原点に")
rep("for (const r of document.querySelectorAll('input[name=pw]')) r.addEventListener('change', draw);",
    "for (const r of document.querySelectorAll('input[name=pw], input[name=rc]')) r.addEventListener('change', draw);")
rep("const resetAll = () => { ticks = 0; hist = PAIRS.map(() => [[0, 0]]);", "const resetAll = () => { ticks = 0; orb = []; orbTick = -1; hist = PAIRS.map(() => [[0, 0]]);")
rep("  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n",
    "  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n"
    "  if (q.get('rc') === '1') document.querySelector('input[name=rc][value=\"1\"]').checked = true;   // rc=1 で右を座標表示（撮影用）\n")
rep("観測者の時空の図は外し、その読みは文字のまま残す。",
    "観測者の時空の図は外し、その読みは文字のまま残す。「右のパワー」を「座標」にすると、右パネルには軌道だけを描く：点 (P sin θ₁, P sin θ₂, P sin θ₃)、P ＝ kR（二乗パワー R² の平方根に表示倍率 k を掛けたもの）、θ_i は三つのパルスの中心位相。Auto のティックごとに 1 点を取り、つないで一本の曲線にする。原点と結ぶ線分・円筒・弧・印は描かない。スライダーを動かすと軌道は消える。左パネル、読める量、枠は変わらない。")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/vx_check.js", "w", encoding="utf-8").write(js)
