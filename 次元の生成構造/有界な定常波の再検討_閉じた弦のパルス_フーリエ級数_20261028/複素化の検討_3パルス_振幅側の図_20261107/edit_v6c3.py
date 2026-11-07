# -*- coding: utf-8 -*-
import io, re
p = "index_v6c3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の複素パルス 3 波（index_v4c3.html のパワー表示 3 択版）</title>",
    "<title>閉じた弦の複素パルス 3 波（index_v5c3.html の振幅側の図版：右パネルを円周の中心を原点にした複素波に）</title>")
rep('<label class="ck"><input type="radio" name="pw" value="2">パワー分離 Re²（半径方向）・Im²（垂直）</label></span>\n',
    '<label class="ck"><input type="radio" name="pw" value="2">パワー分離 Re²（半径方向）・Im²（垂直）</label></span>\n'
    '  <span>右</span>\n'
    '  <span class="grp"><label class="ck"><input type="radio" name="rm" value="0" checked>波（円周の中心を原点に）</label><label class="ck"><input type="radio" name="rm" value="1">線分（原点から |Z|² の山の頂点へ、先端に半値全幅を直径とする円）</label></span>\n')
rep("let Lfix = null;                              // 3D の枠 ±L。null なら次の draw で合わせ直す（スケール調整ボタン）。それ以外では R が変わっても枠は動かさない\n",
    "let Lfix = null;                              // 3D の枠 ±L。null なら次の draw で合わせ直す（スケール調整ボタン）。それ以外では R が変わっても枠は動かさない\n"
    "let Lfix2 = null;                             // 右パネルの枠 ±L2。同じくスケール調整で合わせ直す\n")
rep("  const zrange = [-L, L];                                     // x・y・z を同じ範囲にして目盛の間隔をそろえる\n",
    "  const zrange = [-L, L];                                     // x・y・z を同じ範囲にして目盛の間隔をそろえる\n"
    "  // ---- 右パネル：振幅側から見た図。左と同じ波を、R の円周上でなく円周の中心を原点にして描く（kR = 0。振幅の向きは左と同じ：半径方向に Re／Re²／|Z|²、垂直に Im／Im²）。表示の 3 択は左と共通\n"
    "  const rmode = Number(document.querySelector('input[name=rm]:checked').value);   // 右：0 = 波、1 = 線分\n"
    "  const rad = pmode === 1 ? PW : pmode === 2 ? RE2 : Z.re, ver = pmode === 1 ? null : pmode === 2 ? IM2 : Z.im;\n"
    "  const Xr = Array.from(g.c, (c, i) => rad[i] * c), Yr = Array.from(g.s, (s, i) => rad[i] * s), Zr = ver ? Array.from(ver) : Xr.map(() => 0);\n"
    "  const traces2 = [{type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], scene: 'scene2', marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点（円周の中心）<extra></extra>'}];\n"
    "  const lobes = [];\n"
    "  let L2now = 0.5;\n"
    "  if (rmode === 0) {\n"
    "    traces2.push({type: 'scatter3d', mode: 'lines', x: Xr, y: Yr, z: Zr, scene: 'scene2', line: {color: pmode === 1 ? '#2ca02c' : '#2f6fdf', width: 5}, showlegend: false,\n"
    "      hovertemplate: 'ψ=%{customdata[0]:.2f}°<br>' + (pmode === 1 ? '|Z|²' : pmode === 2 ? 'Re²' : 'Re Z') + '=%{customdata[1]:.4f}（半径方向）' + (ver ? '<br>' + (pmode === 2 ? 'Im²' : 'Im Z') + '=%{z:.4f}（垂直）' : '') + '<extra></extra>', customdata: Array.from(g.psi, (v, i) => [v, rad[i]])});\n"
    "    for (let i = 0; i < n; i++) L2now = Math.max(L2now, Math.abs(Xr[i]), Math.abs(Yr[i]), Math.abs(Zr[i]));\n"
    "  } else {\n"
    "    // 線分：パワー |Z|² の山ごとに、原点から頂点 (P cos ψ*, P sin ψ*, 0) への線分（緑）と、半値全幅の両端（極座標の曲線上の 2 点、橙の点）を結ぶ長さを直径とする円（橙、z = 0 の面、頂点が中心）\n"
    "    const N = n - 1, at = i => ((i % N) + N) % N; let top = 0; for (let i = 0; i < N; i++) if (PW[i] > top) top = PW[i];\n"
    "    const cand = [];\n"
    "    for (let i = 0; i < N; i++) if (PW[i] >= PW[at(i - 1)] && PW[i] > PW[at(i + 1)] && PW[i] >= 0.1 * top) {   // 極大（最大の 1 割以上。単独パルスの Gibbs の副山は 4.3% 以下）\n"
    "      let a = i, b = i; while (PW[at(a - 1)] > PW[i] / 2 && i - a < N) a--; while (PW[at(b + 1)] > PW[i] / 2 && b - i < N) b++;\n"
    "      cand.push({i, a, b});\n"
    "    }\n"
    "    cand.sort((u, v) => PW[v.i] - PW[u.i]);\n"
    "    for (const c of cand) if (!lobes.some(k => ((c.i - k.a) % N + N) % N <= k.b - k.a)) lobes.push(c);   // 高い順に。既に選んだ山の半値幅の中にある極大は捨てる\n"
    "    for (const l of lobes) {\n"
    "      const P = PW[l.i], ps = g.psi[l.i], ia = at(l.a), ib = at(l.b);\n"
    "      const pa = [PW[ia] * g.c[ia], PW[ia] * g.s[ia]], pb = [PW[ib] * g.c[ib], PW[ib] * g.s[ib]];\n"
    "      const d = Math.hypot(pa[0] - pb[0], pa[1] - pb[1]), cx = P * g.c[l.i], cy = P * g.s[l.i];\n"
    "      l.P = P; l.psi = ps; l.d = d; l.w = (l.b - l.a) * 360 / N;\n"
    "      traces2.push({type: 'scatter3d', mode: 'lines+markers', x: [0, cx], y: [0, cy], z: [0, 0], scene: 'scene2', line: {color: '#2ca02c', width: 5}, marker: {size: [0, 4], color: '#2ca02c'}, showlegend: false, hovertemplate: `山 ψ* = ${ps.toFixed(1)}°、P = ${P.toFixed(4)}<extra></extra>`});\n"
    "      const cxs = [], cys = []; for (let k = 0; k <= 72; k++) { const a = 2 * Math.PI * k / 72; cxs.push(cx + d / 2 * Math.cos(a)); cys.push(cy + d / 2 * Math.sin(a)); }\n"
    "      traces2.push({type: 'scatter3d', mode: 'lines', x: cxs, y: cys, z: cxs.map(() => 0), scene: 'scene2', line: {color: '#e8862a', width: 3}, showlegend: false, hovertemplate: `幅（半値全幅 ${l.w.toFixed(1)}°）直径 ${d.toFixed(4)}<extra></extra>`});\n"
    "      traces2.push({type: 'scatter3d', mode: 'markers', x: [pa[0], pb[0]], y: [pa[1], pb[1]], z: [0, 0], scene: 'scene2', marker: {size: 3, color: '#e8862a'}, showlegend: false, hoverinfo: 'skip'});\n"
    "      L2now = Math.max(L2now, P + d / 2);\n"
    "    }\n"
    "  }\n"
    "  L2now *= 1.15;\n"
    "  if (Lfix2 === null) Lfix2 = L2now;\n"
    "  const L2 = Lfix2;\n"
    "  const lobeText = lobes.map(l => `ψ* ${l.psi.toFixed(1)}°, P ${l.P.toFixed(4)}, 幅 ${l.w.toFixed(1)}°, 直径 ${l.d.toFixed(4)}`).join('；') || 'なし';\n")
rep("${Lnow > L ? '、はみ出し' : ''}）`;",
    "${Lnow > L ? '、はみ出し' : ''}）　右枠 ±${L2.toFixed(3)}（必要 ±${L2now.toFixed(3)}${L2now > L2 ? '、はみ出し' : ''}）${rmode === 1 ? '　線分：' + lobeText : ''}`;")
rep("    scene: {domain: {x: [0, 0.6], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // x・y・z を同じ範囲・同じ縮尺に固定\n",
    "    scene: {domain: {x: [0, 0.5], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // x・y・z を同じ範囲・同じ縮尺に固定\n")
old_axes = re.search(r"    xaxis2: \{domain: \[0\.66, 1\].*?\n    yaxis2: \{.*?\n", s, re.S).group(0)
rep(old_axes,
    "    scene2: {domain: {x: [0.5, 1], y: [0, 1]}, aspectmode: 'manual', aspectratio: {x: 1, y: 1, z: 1},   // 右パネル：振幅側から見た図（原点は円周の中心）。向きと縮尺の決め方は左と同じ\n"
    "            xaxis: ax('x', [-L2, L2]), yaxis: ax('y', [-L2, L2]), zaxis: ax('z', [-L2, L2]),\n"
    "            camera: {eye: {x: 1.2, y: -1.45, z: 0.72}, projection: {type: 'orthographic'}}}\n")
rep("  Plotly.react('plot', traces, layout, {responsive: true});",
    "  Plotly.react('plot', traces.filter(t => t.xaxis !== 'x2').concat(traces2), layout, {responsive: true});   // 旧右パネル（観測者の時空）の曲線は描かない。読める量の文字（h2）はそのまま")
rep("$('fit').addEventListener('click', () => { Lfix = null; draw(); });",
    "$('fit').addEventListener('click', () => { Lfix = null; Lfix2 = null; draw(); });")
rep("for (const r of document.querySelectorAll('input[name=pw]')) r.addEventListener('change', draw);",
    "for (const r of document.querySelectorAll('input[name=pw], input[name=rm]')) r.addEventListener('change', draw);")
rep("  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n",
    "  if (q.has('pow')) document.querySelector(`input[name=pw][value=\"${q.get('pow')}\"]`).checked = true;   // pow=1 合計、pow=2 分離（撮影用）\n"
    "  if (q.has('right')) document.querySelector(`input[name=rm][value=\"${q.get('right')}\"]`).checked = true;   // right=1 で右を線分（撮影用）\n")
rep("「パワー分離」は軸ごとのパワー Re²（半径方向、青）と Im²（垂直、橙）を Re・Im と同じ置き方で描く（合計が |Z|²）。",
    "「パワー分離」は軸ごとのパワー Re²（半径方向、青）と Im²（垂直、橙）を Re・Im と同じ置き方で描く（合計が |Z|²）。右パネルは振幅側から見た図：半径 R の方向は振幅側からは見えない方向なので、左と同じ波を、R の円周上でなく円周の中心を原点にして描く（kR ＝ 0。振幅の向きは同じで、半径方向に Re／Re²／|Z|²、垂直に Im／Im²）。右の「線分」は、|Z|² の山ごとに原点から頂点までの線分（緑）と、半値全幅の両端を結ぶ長さを直径とする円（橙）を頂点を中心に置く。観測者の時空の図は外し、その読み（三辺・傾き・合成則など）は文字のまま残す。")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v6c3_check.js", "w", encoding="utf-8").write(js)
