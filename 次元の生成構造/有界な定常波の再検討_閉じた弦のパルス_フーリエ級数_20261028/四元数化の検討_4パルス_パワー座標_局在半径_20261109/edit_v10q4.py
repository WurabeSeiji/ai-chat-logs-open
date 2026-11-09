# -*- coding: utf-8 -*-
# index_v9q4.html → index_v10q4.html：座標表示の球の半径を s × A_i から局在の幅 2/√M_i（度）にする。
# 2026-11-09 木原氏「現在粒子の位置としてプロットしている球体の大きさを 2/√k にしてみては。相互作用を考える前に衝突可能性が理解できる」。
# k ＝ その波の倍音数 M_i。表示座標が度なので 2/√k ラジアンを度に直す（M=100 で 11.46°、M=10 で 36.2°）。
# 変えるのは球の半径の行・注釈・名札だけ。スライダー s は残すが球には効かない。A_i は波形だけ。左パネル・読み・軌道・枠は変更なし。
import io
p = "index_v10q4.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new, n=1):
    global s
    c = s.count(old); assert c == n, (old[:70], c, n); s = s.replace(old, new)
rep("// 現在位置に、その波の色の球。半径 = s × 振幅 A_i（q = A(a + n b) の A、s = スライダー s、1〜100）。k は掛けない",
    "// 現在位置に、その波の色の球。半径 = 局在の幅 2/√M_i（ラジアンを度に、M_i = その波の倍音数。11/9 木原氏：cos² 型分布として読んだ大きさ、球の重なり = 衝突可能性）。s も A_i も掛けない")
rep("r = Number(ss.value) * p.pulses[i].A,", "r = (2 / Math.sqrt(p.pulses[i].M)) * 180 / Math.PI,")
rep("name: `波 ${i + 1} の球（半径 s A = ${r.toFixed(3)}）`", "name: `波 ${i + 1} の球（半径 2/√M = ${r.toFixed(2)}°）`")
io.open(p, "w", encoding="utf-8").write(s); print("ok")
