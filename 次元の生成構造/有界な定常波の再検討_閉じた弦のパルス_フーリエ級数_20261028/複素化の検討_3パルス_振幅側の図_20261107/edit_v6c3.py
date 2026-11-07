# -*- coding: utf-8 -*-
import io, re
p = "index_v6c3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の複素パルス 3 波（index_v4c3.html のパワー表示 3 択版）</title>",
    "<title>閉じた弦の複素パルス 3 波（index_v5c3.html の振幅側の図版：右パネルは左の図を観測者の中心位相が原点に来るよう平行移動）</title>")
rep("  // ---- 右パネル：観測者の時空。観測者 o を原点に置き、他の二つのパルス i を点 P_i = (Δ_io, Φ_io) で描く。\n",
    "  // ---- 右パネル：左の図をそのまま平行移動したもの。観測者 o の中心位相の点 C_o = kR(cos θ_o, sin θ_o) が原点に来るように、全部の描画要素を −C_o だけ動かす\n"
    "  //      （円の中心は反対側の円周上 −C_o へ。位相 ψ・振幅の向き・縮尺はそのまま。計算は変えない）\n"
    "  const oo = obsOf(), dx = -p.Rd * Math.cos(Math.PI * p.pulses[oo].theta / 180), dy = -p.Rd * Math.sin(Math.PI * p.pulses[oo].theta / 180);\n"
    "  const shift = (arr, d) => Array.isArray(arr[0]) ? arr.map(row => Array.from(row, v => v + d)) : Array.from(arr, v => v + d);\n"
    "  const traces2 = traces.map(t => Object.assign({}, t, {scene: 'scene2', x: shift(t.x, dx), y: shift(t.y, dy)}));\n"
    "  // ---- 旧右パネル（観測者の時空）の計算は読める量（h2）のために残す。曲線は描かない\n")
rep("${Lnow > L ? '、はみ出し' : ''}）`;",
    "${Lnow > L ? '、はみ出し' : ''}）　右：観測者 ${oo + 1} の中心位相 θ = ${p.pulses[oo].theta}° を原点に（平行移動 (${dx.toFixed(3)}, ${dy.toFixed(3)})、円の中心は反対側の円周上）`;")
rep("    scene: {domain: {x: [0, 0.6], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // x・y・z を同じ範囲・同じ縮尺に固定\n",
    "    scene: {domain: {x: [0, 0.5], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // x・y・z を同じ範囲・同じ縮尺に固定\n")
old_axes = re.search(r"    xaxis2: \{domain: \[0\.66, 1\].*?\n    yaxis2: \{.*?\n", s, re.S).group(0)
rep(old_axes,
    "    scene2: {domain: {x: [0.5, 1], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // 右パネル：左と同じ枠・同じ向き・同じ縮尺（平行移動しただけ）\n"
    "            xaxis: ax('x', [-L, L]), yaxis: ax('y', [-L, L]), zaxis: ax('z', zrange),\n"
    "            camera: {eye: {x: 1.2, y: -1.45, z: 0.72}, projection: {type: 'orthographic'}}}\n")
rep("  Plotly.react('plot', traces, layout, {responsive: true});",
    "  Plotly.react('plot', traces.filter(t => t.xaxis !== 'x2').concat(traces2), layout, {responsive: true});   // 左の図＋平行移動した右の図。観測者の時空の曲線は描かない")
rep("「パワー分離」は軸ごとのパワー Re²（半径方向、青）と Im²（垂直、橙）を Re・Im と同じ置き方で描く（合計が |Z|²）。",
    "「パワー分離」は軸ごとのパワー Re²（半径方向、青）と Im²（垂直、橙）を Re・Im と同じ置き方で描く（合計が |Z|²）。右パネルは左の図をそのまま平行移動したもの：円の中心と、観測者（右：観測者で選ぶ）の中心位相の点 C_o ＝ kR(cos θ_o, sin θ_o) を結ぶ線分を、半径の分だけ中心の方向へ動かし、C_o を原点に、円の中心を反対側の円周上に置く。位相 ψ・振幅の向き・縮尺はそのままで、パルスの位相幅は原点の周りで半径方向と直交する接線の向きに乗る。観測者の時空の図は外し、その読みは文字のまま残す。")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v6c3_check.js", "w", encoding="utf-8").write(js)
