# -*- coding: utf-8 -*-
# index_v8q3.html → index_v9q4.html：パルスを 4 つにする（2026-11-07 木原氏「3 体実験系を拡張して 4 体に。四元数版だけ」）。
# 行 3 の複製で行 4（θ₄=300, φ₄=1, M₄=100, δ₄=−10, A₄=1, n₄=(1,1,1)/√3、色 #8c564b、印 circle）。対は 12, 23, 13 の後ろに 14, 24, 34。
# 観測者は 1〜4（既定 3 のまま）、標準の対は 6 組（既定 23 のまま）。3 固有の読み（三角形の向き・合成則・比・鏡映対）は 1〜3 のまま残す。
# 座標表示の点は x_i = n_i ∘ (Θ₁, Θ₂, Θ₃) のまま（Θ₄ は軌道に使わない。4 体目の軸の写像は未定）。
import io, re
p = "index_v9q4.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new, n=1):
    global s
    c = s.count(old); assert c == n, (old[:70], c, n); s = s.replace(old, new)
# 題
rep("<title>閉じた弦の四元数パルス 3 波（index_v7q3.html のパワー座標版：右パネルに (i_i Θ₁, j_i Θ₂, k_i Θ₃) の軌道、Θ ＝ 中心位相の累積）</title>",
    "<title>閉じた弦の四元数パルス 4 波（index_v8q3.html の 4 パルス版）</title>")
# 行 4（行 3 の複製）
m = re.search(r'<div class="bar">\n  <span>3:</span>\n.*?<span class="nn" id="nn3"></span>\n</div>\n', s, re.S); assert m
r = m.group(0).replace("<span>3:</span>", "<span>4:</span>").replace("₃", "₄")
for key in ["th", "ph", "mm", "odd", "even", "dl", "aa", "ni", "nj", "nk", "nn"]: r = r.replace('id="%s3' % key, 'id="%s4' % key)
r = r.replace('id="th4" type="range" min="0" max="360" step="1" value="58"><span class="v" id="th4v">58°</span>', 'id="th4" type="range" min="0" max="360" step="1" value="300"><span class="v" id="th4v">300°</span>')
r = r.replace('id="dl4" type="range" min="-30" max="30" step="1" value="10"><span class="v" id="dl4v">10°/t</span>', 'id="dl4" type="range" min="-30" max="30" step="1" value="-10"><span class="v" id="dl4v">-10°/t</span>')
r = r.replace('id="ni4" type="range" min="-1" max="1" step="0.05" value="0"><span class="v" id="ni4v">0</span>', 'id="ni4" type="range" min="-1" max="1" step="0.05" value="1"><span class="v" id="ni4v">1</span>')
r = r.replace('id="nj4" type="range" min="-1" max="1" step="0.05" value="0"><span class="v" id="nj4v">0</span>', 'id="nj4" type="range" min="-1" max="1" step="0.05" value="1"><span class="v" id="nj4v">1</span>')
assert r.count('value="300"') == 1 and r.count('value="-10"') == 1 and r.count('id="ni4v">1<') == 1 and r.count('id="nj4v">1<') == 1 and r.count('id="nk4v">1<') == 1 and not re.search(r'id="\w+3v?"', r) and "₃" not in r and "<span>3:</span>" not in r
s = s[:m.end()] + r + s[m.end():]
# 観測者・標準の対
rep('<label class="ck"><input type="radio" name="obs" value="2" checked>3</label></span>', '<label class="ck"><input type="radio" name="obs" value="2" checked>3</label><label class="ck"><input type="radio" name="obs" value="3">4</label></span>')
rep('<label class="ck"><input type="radio" name="std" value="2">13</label></span>', '<label class="ck"><input type="radio" name="std" value="2">13</label><label class="ck"><input type="radio" name="std" value="3">14</label><label class="ck"><input type="radio" name="std" value="4">24</label><label class="ck"><input type="radio" name="std" value="5">34</label></span>')
# 配列の長さ
rep("const PAIRS = [[0, 1], [1, 2], [0, 2]], PAIRC = ['#e8862a', '#9467bd', '#2ca02c'], PAIRL = ['12', '23', '13'];",
    "const PAIRS = [[0, 1], [1, 2], [0, 2], [0, 3], [1, 3], [2, 3]], PAIRC = ['#e8862a', '#9467bd', '#2ca02c', '#8c564b', '#e377c2', '#17becf'], PAIRL = ['12', '23', '13', '14', '24', '34'];   // 先頭 3 対は 3 パルス版と同じ順（q12, q23, q13 の分解をそのまま使う）")
rep("let cum = [0, 1, 2].map(() => [[0, 0]]);", "let cum = [0, 1, 2, 3].map(() => [[0, 0]]);")
rep("let base = [0, 1, 2].map(() => [0, 0]);", "let base = [0, 1, 2, 3].map(() => [0, 0]);")
rep("const P = [1, 2, 3].map(i => ({", "const P = [1, 2, 3, 4].map(i => ({")
rep("orb = [[], [], []]", "orb = [[], [], [], []]", 2)
rep("const Th = [0, 1, 2].map(i => base[i][0] + cum[i][cum[i].length - 1][0]);", "const Th = [0, 1, 2, 3].map(i => base[i][0] + cum[i][cum[i].length - 1][0]);   // Θ₄ も累積するが、軌道の点は n_i ∘ (Θ₁, Θ₂, Θ₃) のまま（4 体目の軸の写像は未定）")
rep("const others = [0, 1, 2].filter(i => i !== o);", "const others = [0, 1, 2, 3].filter(i => i !== o);")
rep("cum = [0, 1, 2].map(() => [[0, 0]]); win0 = null;", "cum = [0, 1, 2, 3].map(() => [[0, 0]]); win0 = null;")
rep("['#1f77b4', '#9467bd', '#2ca02c']", "['#1f77b4', '#9467bd', '#2ca02c', '#8c564b']", 5)
rep("['diamond', 'square', 'cross']", "['diamond', 'square', 'cross', 'circle']", 4)
rep("z: [0, 0, 0, 0], line: {color: '#8a8a8f', width: 1.5, dash: 'dot'}", "z: [0, 0, 0, 0, 0], line: {color: '#8a8a8f', width: 1.5, dash: 'dot'}")
# 観測者の時空：他の体は 3 つ、辺は 6 本
m2 = re.search(r"^  const \[a, b\] = others, Pa = curP\[0\], Pb = curP\[1\], Pab = .*$", s, re.M); assert m2
s = s[:m2.start()] + "  const edges = others.map((i, q) => [`${o + 1}–${i + 1}`, curP[q]]); for (let q1 = 0; q1 < others.length; q1++) for (let q2 = q1 + 1; q2 < others.length; q2++) { const a = others[q1], b = others[q2]; edges.push([`${a + 1}–${b + 1}`, [rlen(relv(a, b, nT - 1)) / (stdOK ? Lstd : 1), curP[q1][1] - curP[q2][1]]]); }   // 辺：観測者–他の各体と、他の体どうし。距離は |vec x_a − vec x_b|（(r, t) 面の差ではない）、時計のずれは t の差" + s[m2.end():]
rep("  // 三辺：観測者–a、観測者–b、a–b。各辺に距離（x の差）と時計のずれ（t の差）", "  // 辺：観測者–他の各体（3 本）と、他の体どうし（3 本）。各辺に距離（x の差）と時計のずれ（t の差）")
rep("　三辺　${side(`${o + 1}–${a + 1}`, Pa)}　${side(`${o + 1}–${b + 1}`, Pb)}　${side(`${a + 1}–${b + 1}`, Pab)}　｜", "　辺 ${edges.length} 本　${edges.map(e => side(e[0], e[1])).join('　')}　｜")
rep("x: [0, curP[0][0], curP[1][0], 0], y: [0, curP[0][1], curP[1][1], 0]", "x: [0, ...curP.map(v => v[0]), 0], y: [0, ...curP.map(v => v[1]), 0]")
# 説明
rep("（i, j, k。index_v8q3 で変更、index_v7q3 までは全部 k）、観測者 3、標準 23（2026-11-07 木原氏の指定。",
    "（i, j, k。index_v8q3 で変更、index_v7q3 までは全部 k）、パルス 4（index_v9q4 で追加）：θ₄ ＝ 300、φ₄ ＝ 1、M₄ ＝ 100、δ₄ ＝ −10、A₄ ＝ 1、n₄ ＝ (1, 1, 1)/√3、色 赤茶・印 ●。観測者 3、標準 23（2026-11-07 木原氏の指定。")
rep("Auto で軌道に合わせて計算し直す</label></span>", "Auto で軌道に合わせて計算し直す。4 波でも点は x_i ＝ n_i ∘ (Θ₁, Θ₂, Θ₃) のまま（Θ₄ は軌道に使わない。4 体目の軸の写像は未定）</label></span>")
rep("'n' + '₁₂₃'[i]", "'n' + '₁₂₃₄'[i]")   # n の添字（3 文字だと 4 本目が undefined になる）
io.open(p, "w", encoding="utf-8").write(s); print("html ok")
js = re.search(r"<script>\n(.*?)</script>", s, re.S).group(1); open("/tmp/vx_check.js", "w", encoding="utf-8").write(js)
