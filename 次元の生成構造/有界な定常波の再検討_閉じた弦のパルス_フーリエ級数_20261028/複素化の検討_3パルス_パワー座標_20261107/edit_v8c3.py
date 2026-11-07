# -*- coding: utf-8 -*-
# index_v7c3.html → index_v8c3.html：表示に 4 つ目「パワー座標」を足す（四元数版 index_v8q3 と同じ）。複素数の軸は i の 1 本なので x = P i だけ：
# P_i = √|Z_i|² を x 軸に置き、左は和 Σ_i P_i(ψ) の軌道（x 軸上の往復）、右は各パルスの P_i（山の範囲）
import io, re
p = "index_v8c3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の複素パルス 3 波（index_v6c3.html のクリップ版：右パネルは各パルスの山だけ描く）</title>",
    "<title>閉じた弦の複素パルス 3 波（index_v7c3.html のパワー座標版：P_i ＝ √|Z_i|² を軸 i に置いた軌道）</title>")
rep('<label class="ck"><input type="radio" name="pw" value="2">パワー分離 Re²（半径方向）・Im²（垂直）</label></span>\n',
    '<label class="ck"><input type="radio" name="pw" value="2">パワー分離 Re²（半径方向）・Im²（垂直）</label><label class="ck"><input type="radio" name="pw" value="3">パワー座標 x ＝ P i（P ＝ √|Z_i|²、複素数の軸は i の 1 本。左は Σ_i P_i の軌道、右は各パルスの P_i）</label></span>\n')
rep("  const [rmaxP, rmeanP] = stat(RE2), [imaxP, imeanP] = stat(IM2);\n",
    "  const [rmaxP, rmeanP] = stat(RE2), [imaxP, imeanP] = stat(IM2);\n"
    "  // パワー座標：各パルスの P_i = √|Z_i|² を軸 i の向き（x 軸）に置く（x = P i。複素数の軸は 1 本）。左は和 S(ψ) = Σ_i P_i(ψ) の軌道、右は各パルスの P_i（山の範囲だけ）。計算は変えない\n"
    "  const Pn = Zs.map(Zi => Zi.re.map((v, k) => Math.hypot(v, Zi.im[k])));\n"
    "  const SX = new Float64Array(n); Pn.forEach(P => { for (let k = 0; k < n; k++) SX[k] += P[k]; });\n"
    "  let smax = 0; for (let k = 0; k < n; k++) smax = Math.max(smax, SX[k]);\n")
rep("  const Lnow = (pmode === 1 ? Math.max(p.Rd + pmax, p.Rd + 0.5) : ",
    "  const Lnow = (pmode === 3 ? Math.max(smax, 0.5) : pmode === 1 ? Math.max(p.Rd + pmax, p.Rd + 0.5) : ")
rep("  ];\n  // ---- 右パネル：複素数は全て原点にある。",
    "  ];\n"
    "  if (pmode === 3) {                                           // パワー座標（左）：和の軌道 S(ψ) = Σ_i P_i を x 軸上に 1 本（緑）と原点。円周・円筒・印は描かない\n"
    "    traces.length = 0;\n"
    "    traces.push({type: 'scatter3d', mode: 'lines', x: Array.from(SX), y: SX.map(() => 0), z: SX.map(() => 0), line: {color: '#2ca02c', width: 5}, showlegend: false, hovertemplate: 'ψ=%{customdata:.2f}°<br>Σ P_i = %{x:.4f}（x 軸 = i）<extra></extra>', customdata: Array.from(g.psi)});\n"
    "    traces.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点<extra></extra>'});\n"
    "  }\n"
    "  // ---- 右パネル：複素数は全て原点にある。")
rep("    const C = idx.map(k => g.c[k]), S = idx.map(k => g.s[k]), PSI = idx.map(k => g.psi[k]);\n",
    "    const C = idx.map(k => g.c[k]), S = idx.map(k => g.s[k]), PSI = idx.map(k => g.psi[k]);\n"
    "    if (pmode === 3) {                                           // パワー座標（右）：各パルスの P_i を x 軸上に山の範囲だけ（緑）。円周・円筒・平行移動なし\n"
    "      const P = Pn[i];\n"
    "      traces2.push({type: 'scatter3d', mode: 'lines', x: idx.map(k => P[k]), y: idx.map(() => 0), z: idx.map(() => 0), scene: 'scene2', line: {color: '#2ca02c', width: 5}, showlegend: false,\n"
    "        hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>P=%{customdata[1]:.4f}（x 軸 = i）<extra></extra>`, customdata: PSI.map((v, j) => [v, P[idx[j]]])});\n"
    "      return;\n"
    "    }\n")
rep("` : ''}　枠 ±", "` : pmode === 3 ? `　パワー座標：P_i = √|Z_i|² を軸 i（x）に、max Σ P_i = ${smax.toFixed(4)}` : ''}　枠 ±")
rep("観測者の時空の図は外し、その読みは文字のまま残す。",
    "観測者の時空の図は外し、その読みは文字のまま残す。表示の 4 つ目「パワー座標」は、ψ を描画の座標に使わず、各パルスのパワー P_i ＝ √|Z_i|² を軸 i の向き（x 軸）に置いた点を座標にする（四元数版 x ＝ P i, y ＝ P j, z ＝ P k の、軸が i の 1 本の場合）：左は和 S(ψ) ＝ Σ_i P_i(ψ) の軌道を ψ に沿って一本の曲線で結び（x 軸上の往復）、右は各パルスの P_i を山の範囲だけ描く。円周・半径 R・円筒は出てこない。")
rep("// pow=1 合計、pow=2 分離（撮影用）", "// pow=1 合計、pow=2 分離、pow=3 座標（撮影用）")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v8c3_check.js", "w", encoding="utf-8").write(js)
