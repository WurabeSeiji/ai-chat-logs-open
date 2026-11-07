# -*- coding: utf-8 -*-
# index_v7c3.html → index_v8c3.html：右パネルのパワー表示だけに「円周上／座標」の切替を足す（四元数版 v8q3 と同じ）。複素数の軸は i の 1 本なので x = P i だけ（P_i = √|Z_i|²）
import io, re
p = "index_v8c3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の複素パルス 3 波（index_v6c3.html のクリップ版：右パネルは各パルスの山だけ描く）</title>",
    "<title>閉じた弦の複素パルス 3 波（index_v7c3.html のパワー座標版：右パネルのパワーを P_i i の軌道で描く切替）</title>")
rep('<label class="ck"><input type="radio" name="pw" value="2">パワー分離 Re²（半径方向）・Im²（垂直）</label></span>\n',
    '<label class="ck"><input type="radio" name="pw" value="2">パワー分離 Re²（半径方向）・Im²（垂直）</label></span>\n'
    '  <span>右のパワー</span>\n'
    '  <span class="grp"><label class="ck"><input type="radio" name="rc" value="0" checked>円周上</label><label class="ck"><input type="radio" name="rc" value="1">座標 x ＝ P i（P ＝ √|Z_i|²、複素数の軸は i の 1 本）</label></span>\n')
rep("  const pmode = Number(document.querySelector('input[name=pw]:checked').value), usePow = pmode > 0;   // 表示：0 = Re・Im、1 = パワー合計 |Z|²、2 = パワー分離 Re²・Im²（計算は変えない）\n",
    "  const pmode = Number(document.querySelector('input[name=pw]:checked').value), usePow = pmode > 0;   // 表示：0 = Re・Im、1 = パワー合計 |Z|²、2 = パワー分離 Re²・Im²（計算は変えない）\n"
    "  const rcoord = Number(document.querySelector('input[name=rc]:checked').value) === 1;              // 右のパワー：false = 円周上（従来）、true = 座標 P_i i の軌道。左には影響しない\n")
old_block = ("    traces2.push({type: 'scatter3d', mode: 'lines', x: C.map((c, j) => (p.Rd + radI[idx[j]]) * c + ddx), y: S.map((s, j) => (p.Rd + radI[idx[j]]) * s + ddy), z: Xc.map(() => 0), scene: 'scene2', line: {color: pmode === 1 ? '#2ca02c' : '#2f6fdf', width: 5},\n"
             "      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>` + (pmode === 1 ? '|Z|²' : pmode === 2 ? 'Re²' : 'Re Z') + '=%{customdata[1]:.4f}（半径方向）<extra></extra>', customdata: PSI.map((v, j) => [v, radI[idx[j]]])});\n"
             "    if (verI) traces2.push({type: 'scatter3d', mode: 'lines', x: Xc, y: Yc, z: idx.map(k => verI[k]), scene: 'scene2', line: {color: '#e8862a', width: 5},\n"
             "      showlegend: false, hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata:.2f}°<br>` + (pmode === 2 ? 'Im²' : 'Im Z') + '=%{z:.4f}（垂直）<extra></extra>', customdata: PSI});\n")
new_block = ("    if (rcoord && pmode > 0) {                                   // 右のパワーを座標表示：パワーの曲線だけを P_i i の軌道（P_i = √|Z_i|²、x = P i）に置き換える。山の範囲だけ。他は同じ\n"
             "      traces2.push({type: 'scatter3d', mode: 'lines', x: idx.map(k => Math.sqrt(Pi[k])), y: idx.map(() => 0), z: idx.map(() => 0), scene: 'scene2', line: {color: '#2ca02c', width: 5}, showlegend: false,\n"
             "        hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>P=%{customdata[1]:.4f}（x 軸 = i）<extra></extra>`, customdata: PSI.map((v, j) => [v, Math.sqrt(Pi[idx[j]])])});\n"
             "    } else {\n" + old_block + "    }\n")
rep(old_block, new_block)
rep("for (const r of document.querySelectorAll('input[name=pw]')) r.addEventListener('change', draw);",
    "for (const r of document.querySelectorAll('input[name=pw], input[name=rc]')) r.addEventListener('change', draw);")
rep("  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n",
    "  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n"
    "  if (q.get('rc') === '1') document.querySelector('input[name=rc][value=\"1\"]').checked = true;   // rc=1 で右のパワーを座標表示（撮影用）\n")
rep("観測者の時空の図は外し、その読みは文字のまま残す。",
    "観測者の時空の図は外し、その読みは文字のまま残す。「右のパワー」を「座標」にすると、右パネルのパワーの曲線だけを、各パルスのパワー P_i ＝ √|Z_i|² を軸 i に置いた点 x ＝ P i の軌道（ψ に沿って一本の曲線、山の範囲だけ。複素数の軸は 1 本）に置き換える。左パネル、右の円筒・弧・平行移動・印、読める量、枠は変わらない。")
io.open(p, "w", encoding="utf-8").write(s); print("c html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v8c3_check.js", "w", encoding="utf-8").write(js)
