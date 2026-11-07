# -*- coding: utf-8 -*-
# index_v7q3.html → index_v8q3.html：表示に 4 つ目「パワー座標」を足す。P_i = √|Q_i|² を軸 n_i の向きに置き x = P i, y = P j, z = P k。左は和 Σ_i P_i n_i の軌道、右は各パルスの P_i n_i（山の範囲）
import io, re
p = "index_v8q3.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new):
    global s
    assert s.count(old) == 1, (old[:70], s.count(old)); s = s.replace(old, new)
rep("<title>閉じた弦の四元数パルス 3 波（index_v5q3.html の振幅側の図版：右パネルは各パルスを自分の中心位相が原点に来るよう平行移動し、山だけ描く）</title>",
    "<title>閉じた弦の四元数パルス 3 波（index_v7q3.html のパワー座標版：P_i ＝ √|Q_i|² を軸 n_i の向きに置いた軌道）</title>")
rep('<label class="ck"><input type="radio" name="pw" value="2">パワー分離 軸ごと W²・V_j²（半径方向）・V_i²・V_k²（垂直）</label></span>\n',
    '<label class="ck"><input type="radio" name="pw" value="2">パワー分離 軸ごと W²・V_j²（半径方向）・V_i²・V_k²（垂直）</label><label class="ck"><input type="radio" name="pw" value="3">パワー座標 x ＝ P i, y ＝ P j, z ＝ P k（P ＝ √|Q_i|²、軸 n_i の向きに。左は Σ_i P_i n_i の軌道、右は各パルスの P_i n_i）</label></span>\n')
rep("  const PWt = Array.from(PWw, (v, i) => v + PWi[i] + PWj[i] + PWk[i]); const [tmaxP, tmeanP] = stat(PWt);   // 合計 |Q|²\n",
    "  const PWt = Array.from(PWw, (v, i) => v + PWi[i] + PWj[i] + PWk[i]); const [tmaxP, tmeanP] = stat(PWt);   // 合計 |Q|²\n"
    "  // パワー座標：各パルスの P_i = √|Q_i|² を軸 n_i の向きに置く（x = P i, y = P j, z = P k）。左は和 S(ψ) = Σ_i P_i(ψ) n_i の軌道、右は各パルスの P_i n_i（山の範囲だけ）。計算は変えない\n"
    "  const Pn = Zs.map(Zi => Zi.re.map((v, k) => Math.hypot(v, Zi.im[k])));\n"
    "  const SX = new Float64Array(n), SY = new Float64Array(n), SZ = new Float64Array(n);\n"
    "  Pn.forEach((P, q) => { const nv = p.pulses[q].nv; for (let k = 0; k < n; k++) { SX[k] += P[k] * nv[0]; SY[k] += P[k] * nv[1]; SZ[k] += P[k] * nv[2]; } });\n"
    "  let smax = 0; for (let k = 0; k < n; k++) smax = Math.max(smax, Math.hypot(SX[k], SY[k], SZ[k]));\n")
rep("  const Lnow = (pmode === 1 ? Math.max(p.Rd + tmaxP, p.Rd + 0.5) : ",
    "  const Lnow = (pmode === 3 ? Math.max(smax, 0.5) : pmode === 1 ? Math.max(p.Rd + tmaxP, p.Rd + 0.5) : ")
rep("  ];\n  // ---- 右パネル：複素数（四元数）は全て原点にある。",
    "  ];\n"
    "  if (pmode === 3) {                                           // パワー座標（左）：和の軌道 S(ψ) = Σ_i P_i n_i を 1 本（緑）と原点。円周・円筒・印は描かない\n"
    "    traces.length = 0;\n"
    "    traces.push({type: 'scatter3d', mode: 'lines', x: Array.from(SX), y: Array.from(SY), z: Array.from(SZ), line: {color: '#2ca02c', width: 5}, showlegend: false, hovertemplate: 'ψ=%{customdata:.2f}°<br>Σ P_i n_i = (%{x:.4f}, %{y:.4f}, %{z:.4f})<extra></extra>', customdata: Array.from(g.psi)});\n"
    "    traces.push({type: 'scatter3d', mode: 'markers', x: [0], y: [0], z: [0], marker: {symbol: 'x', size: 2, color: fg}, name: '原点', hovertemplate: '原点<extra></extra>'});\n"
    "  }\n"
    "  // ---- 右パネル：複素数（四元数）は全て原点にある。")
rep("    const C = idx.map(k => g.c[k]), S = idx.map(k => g.s[k]), PSI = idx.map(k => g.psi[k]);\n",
    "    const C = idx.map(k => g.c[k]), S = idx.map(k => g.s[k]), PSI = idx.map(k => g.psi[k]);\n"
    "    if (pmode === 3) {                                           // パワー座標（右）：各パルスの P_i n_i を山の範囲だけ（緑）。円周・円筒・平行移動なし（原点が共通）\n"
    "      const P = Pn[i];\n"
    "      traces2.push({type: 'scatter3d', mode: 'lines', x: idx.map(k => P[k] * nv[0]), y: idx.map(k => P[k] * nv[1]), z: idx.map(k => P[k] * nv[2]), scene: 'scene2', line: {color: '#2ca02c', width: 5}, showlegend: false,\n"
    "        hovertemplate: `パルス ${i + 1}<br>ψ=%{customdata[0]:.2f}°<br>P=%{customdata[1]:.4f}（n = ${v3(nv, 3)}）<extra></extra>`, customdata: PSI.map((v, j) => [v, P[idx[j]]])});\n"
    "      return;\n"
    "    }\n")
rep("` : ''}　枠 ±", "` : pmode === 3 ? `　パワー座標：P_i = √|Q_i|² を n_i の向きに、max |Σ P_i n_i| = ${smax.toFixed(4)}` : ''}　枠 ±")
rep("観測者の時空の図は外し、その読みは文字のまま残す。",
    "観測者の時空の図は外し、その読みは文字のまま残す。表示の 4 つ目「パワー座標」は、ψ を描画の座標に使わず、各パルスのパワー P_i ＝ √|Q_i|² を軸 n_i の向きに置いた点 (P_i n_i) を座標にする：左は和 S(ψ) ＝ Σ_i P_i(ψ) n_i の軌道を ψ に沿って一本の曲線で結び、右は各パルスの P_i n_i を山の範囲だけ描く。円周・半径 R・円筒は出てこない。軸が全部同じなら軌道は一本の軸の上を往復し、軸が違えば空間の曲線になる。")
rep("// pow=1 合計、pow=2 分離（撮影用）", "// pow=1 合計、pow=2 分離、pow=3 座標（撮影用）")
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/v8q3_check.js", "w", encoding="utf-8").write(js)
