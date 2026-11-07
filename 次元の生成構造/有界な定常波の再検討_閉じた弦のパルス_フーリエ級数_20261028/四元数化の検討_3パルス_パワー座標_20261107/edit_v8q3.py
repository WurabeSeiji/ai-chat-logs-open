# -*- coding: utf-8 -*-
# index_v7q3.html → index_v8q3.html：右パネルのパワー表示だけに「円周上／座標」の切替を足す。座標のときは、右のパワーの曲線だけを
# P_i n_i（P_i = √|Q_i|²、x = P i, y = P j, z = P k）の軌道に置き換える。左パネル・右の他の要素（円筒・弧・平行移動・クリップ・印）・読める量・枠は一切変えない
import io, re
p = "index_v8q3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の四元数パルス 3 波（index_v5q3.html の振幅側の図版：右パネルは各パルスを自分の中心位相が原点に来るよう平行移動し、山だけ描く）</title>",
    "<title>閉じた弦の四元数パルス 3 波（index_v7q3.html のパワー座標版：右パネルのパワーを P_i n_i の軌道で描く切替）</title>")
rep('<label class="ck"><input type="radio" name="pw" value="2">パワー分離 軸ごと W²・V_j²（半径方向）・V_i²・V_k²（垂直）</label></span>\n',
    '<label class="ck"><input type="radio" name="pw" value="2">パワー分離 軸ごと W²・V_j²（半径方向）・V_i²・V_k²（垂直）</label></span>\n'
    '  <span>右のパワー</span>\n'
    '  <span class="grp"><label class="ck"><input type="radio" name="rc" value="0" checked>円周上</label><label class="ck"><input type="radio" name="rc" value="1">座標 x ＝ P i, y ＝ P j, z ＝ P k（P ＝ √|Q_i|²、n_i ＝ (i, j, k)）</label></span>\n')
rep("  const pmode = Number(document.querySelector('input[name=pw]:checked').value), usePow = pmode > 0;   // 表示：0 = 4 成分、1 = パワー合計 |Q|²、2 = パワー分離 軸ごと（計算は変えない）\n",
    "  const pmode = Number(document.querySelector('input[name=pw]:checked').value), usePow = pmode > 0;   // 表示：0 = 4 成分、1 = パワー合計 |Q|²、2 = パワー分離 軸ごと（計算は変えない）\n"
    "  const rcoord = Number(document.querySelector('input[name=rc]:checked').value) === 1;              // 右のパワー：false = 円周上（従来）、true = 座標 P_i n_i の軌道。左には影響しない\n")
rep("    if (pmode === 1) {                                           // パワー合計：|Q_i|² を半径方向に 1 本（緑）\n",
    "    if (rcoord && pmode > 0) {                                   // 右のパワーを座標表示：パワーの曲線だけを P_i n_i の軌道（P_i = √|Q_i|²、x = P i, y = P j, z = P k）に置き換える。山の範囲だけ。他は同じ\n"
    "      traces2.push({type: 'scatter3d', mode: 'lines', x: idx.map(k => Math.sqrt(Pi[k]) * nv[0]), y: idx.map(k => Math.sqrt(Pi[k]) * nv[1]), z: idx.map(k => Math.sqrt(Pi[k]) * nv[2]), scene: 'scene2', line: {color: '#2ca02c', width: 5}, showlegend: false,\n"
    "        hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>P=%{customdata[1]:.4f}（n = ${v3(nv, 3)}）<extra></extra>`, customdata: PSI.map((v, j) => [v, Math.sqrt(Pi[idx[j]])])});\n"
    "    } else if (pmode === 1) {                                    // パワー合計：|Q_i|² を半径方向に 1 本（緑）\n")
rep("for (const r of document.querySelectorAll('input[name=pw]')) r.addEventListener('change', draw);",
    "for (const r of document.querySelectorAll('input[name=pw], input[name=rc]')) r.addEventListener('change', draw);")
rep("  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n",
    "  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n"
    "  if (q.get('rc') === '1') document.querySelector('input[name=rc][value=\"1\"]').checked = true;   // rc=1 で右のパワーを座標表示（撮影用）\n")
rep("観測者の時空の図は外し、その読みは文字のまま残す。",
    "観測者の時空の図は外し、その読みは文字のまま残す。「右のパワー」を「座標」にすると、右パネルのパワーの曲線だけを、各パルスのパワー P_i ＝ √|Q_i|² を軸 n_i ＝ (i, j, k) に掛けた点 x ＝ P i、y ＝ P j、z ＝ P k の軌道（ψ に沿って一本の曲線、山の範囲だけ）に置き換える。左パネル、右の円筒・弧・平行移動・印、読める量、枠は変わらない。")
io.open(p, "w", encoding="utf-8").write(s); print("q html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v8q3_check.js", "w", encoding="utf-8").write(js)
