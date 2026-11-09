# -*- coding: utf-8 -*-
# index_v11q4.html → index_v12q4i.html：相互作用版（原本からの一回の再生成。2026-11-09 夜、木原氏の叱責のあと作り直し）。
# 入れるもの（算出根拠 v0.2 §8 の照合表）：κ スライダー、振幅四元数 a_i = A_i e^{n_i φ₀i} の動的状態（nDyn, aDyn）、核 F_ij（原本の波形どうしの重なり、直流を含む）、
#   二チャネル交換の全対同時指数、交換の二方式（率／全量、木原氏「スイッチ付けて両方」）、読み戻し、記録、自己検査（?test=1）。
# 触らないもの：パルスの定義（倍音の係数式・幅 φ・帯域 M・偶奇）、全スライダーの既定値、位置 θ・Θ、速さ v、左右の描画、読み。
# 外したもの：cos² 型の係数（推測を実装にした誤り）、包絡スイッチ、既定の δ の変更。
import io
p = "index_v12q4i.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new, n=1):
    global s
    c = s.count(old); assert c == n, (old[:90], c, n); s = s.replace(old, new)
rep("<title>閉じた弦の四元数パルス 4 波（index_v8q3.html の 4 パルス版）</title>", "<title>閉じた弦の四元数パルス 4 波 相互作用版（index_v11q4.html に重なりの核と二チャネル交換）</title>")
rep('  <div class="sl"><span>s</span><input id="ss" type="range" min="1" max="100" step="1" value="10"><span class="v" id="ssv">10</span></div>\n',
    '  <div class="sl"><span>s</span><input id="ss" type="range" min="1" max="100" step="1" value="10"><span class="v" id="ssv">10</span></div>\n'
    '  <div class="sl"><span>κ 相互作用</span><input id="kap" type="range" min="0" max="1" step="0.01" value="1"><span class="v" id="kapv">1.00</span></div>\n'
    '  <span class="grp">交換 <label class="ck"><input type="radio" name="xm" value="rate" checked>率（通過で合計 Δ_F）</label><label class="ck"><input type="radio" name="xm" value="full">全量（毎 tick F Δ_F）</label></span>\n')
rep("<script>\nconst $ = id => document.getElementById(id);", '<pre id="ktest" style="font-size:12px;white-space:pre-wrap;margin:4px 8px"></pre>\n<script>\nconst $ = id => document.getElementById(id);')
rep("const rr = $('rr'), ss = $('ss'), btn = $('auto');   // ss = 球の倍率 s（1〜100、初期 10）\n",
    "const rr = $('rr'), ss = $('ss'), btn = $('auto'), kap = $('kap');   // ss = 球の倍率 s（1〜100、初期 10）、kap = 相互作用の強さ κ（0〜1、0 で相互作用なし版と同一）\n"
    "let nDyn = null, aDyn = null, ihist = [], lastX = null;   // 交換で変わる軸 n_i と振幅 A_i（スライダーは初期値）、交換の記録\n"
    "const xmOf = () => document.querySelector('input[name=xm]:checked').value;   // 交換の方式：rate（率）| full（全量）\n")
rep("  const pulses = P.map(u => {\n", "  const pulses = P.map((u, i) => {\n")
rep("    const A = Number(u.aa.value); co.cs = co.cs.map(c => A * c); co.c0 *= A;", "    const A = aDyn ? aDyn[i] : Number(u.aa.value); co.cs = co.cs.map(c => A * c); co.c0 *= A;   // A_i は交換で変わる（aDyn）")
rep("    return {theta, phi, v, M, delta, phase0, co, A, nv: nvec(u), R: Math.sqrt(R2), R2,", "    return {theta, phi, v, M, delta, phase0, co, A, nv: nDyn ? nDyn[i] : nvec(u), R: Math.sqrt(R2), R2,")
rep("  $('rrv').textContent = p.k; $('ssv').textContent = ss.value;\n", "  $('rrv').textContent = p.k; $('ssv').textContent = ss.value; $('kapv').textContent = Number(kap.value).toFixed(2);\n")
rep("δ = ${u.delta}°/t, M = ${u.M}（奇", "δ = ${u.delta}°/t, φ₀ = ${u.phase0.toFixed(1)}°, M = ${u.M}（奇")
rep("  $('h1').textContent = `${desc}　　和：R = ", "  $('h1').textContent = `${desc}　　${xdesc()}　　和：R = ")
FUNCS = r'''
// ---------- 相互作用（算出根拠_相互作用版 v0.2 §2〜§4）。位置 θ・Θ・速さ v・包絡は触らず、振幅四元数 a_i = A_i e^{n_i φ₀i} だけを変える ----------
const mmul = (A, B) => A.map(r => B[0].map((_, j) => r.reduce((s, x, k) => s + x * B[k][j], 0)));
const Rmat = h => { const E = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]; const cols = E.map(e => qmul(e, h)); return [0, 1, 2, 3].map(r => [0, 1, 2, 3].map(c => cols[c][r])); };   // vec(q h) = Rmat(h) vec(q)（右乗）
function expmAnti(G) {   // 反対称行列の指数：スケーリングと二乗法＋テイラー 16 項。結果は直交行列
  const N = G.length; let nrm = 0; for (const row of G) { let t = 0; for (const x of row) t += Math.abs(x); nrm = Math.max(nrm, t); }
  let sq = 0; while (nrm > 0.25) { nrm /= 2; sq++; }
  const sc = Math.pow(2, -sq), A = G.map(r => r.map(x => x * sc)), I = G.map((r, i) => r.map((_, j) => i === j ? 1 : 0));
  let U = I.map(r => r.slice()), T = I.map(r => r.slice());
  for (let k = 1; k <= 16; k++) { T = mmul(T, A).map(r => r.map(x => x / k)); U = U.map((r, i) => r.map((x, j) => x + T[i][j])); }
  for (let t = 0; t < sq; t++) U = mmul(U, U);
  return U;
}
function xcorr(coI, coJ, dDeg) {   // 原本の波形どうしの重なり ⟨S_i S_j⟩(d) = c₀i c₀j + ½ Σ_m c_i,m c_j,m cos(m d)（共有する m、直流を含む。背景は対の核に入らない）
  const mapJ = new Map(); coJ.ms.forEach((m, k) => mapJ.set(m, coJ.cs[k])); const d = dDeg * Math.PI / 180; let t = coI.c0 * coJ.c0;
  coI.ms.forEach((m, k) => { const cj = mapJ.get(m); if (cj !== undefined) t += 0.5 * coI.cs[k] * cj * Math.cos(m * d); }); return t;
}
function exchangeStep(S) {   // S = {Th, n, v, A, phi0（度）, co, mode}。純関数。返り値 {A, n, phi0, F, dF, G, U, sumA2before, sumA2after}
  const N = S.n.length, k = Number(kap.value);
  const a = S.A.map((A, i) => qexp(S.n[i], S.phi0[i] * Math.PI / 180).map(x => A * x));                  // 振幅四元数 a_i = A_i e^{n_i φ₀i}
  const F = S.n.map(() => new Array(N).fill(0)), dF = S.n.map(() => new Array(N).fill(0)), G = Array.from({length: 4 * N}, () => new Array(4 * N).fill(0));
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
    const ri = S.n[i].map(x => x * S.Th[i]), rj = S.n[j].map(x => x * S.Th[j]);
    const d = Math.hypot(ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]);                                      // 図の距離 d_ij = |Θ_i n_i − Θ_j n_j|（度、巻数込み）
    const c = dot(S.n[i], S.n[j]), cii = xcorr(S.co[i], S.co[i], 0), cjj = xcorr(S.co[j], S.co[j], 0);
    const Ft = (cii > 0 && cjj > 0) ? xcorr(S.co[i], S.co[j], d) / Math.sqrt(cii * cjj) : 0;
    const Fij = 0.5 * (1 + c) * Ft; F[i][j] = F[j][i] = Fij;                                                 // 核 F_ij = (1 + n_i·n_j)/2 · ⟨S_i S_j⟩(d)/√(⟨S_i²⟩⟨S_j²⟩)
    const us = [S.n[i][0] + S.n[j][0], S.n[i][1] + S.n[j][1], S.n[i][2] + S.n[j][2]], L = Math.hypot(us[0], us[1], us[2]);
    if (L < 1e-9 || Math.abs(Fij) < 1e-300) continue;
    const u = us.map(x => x / L);                                                                            // 共有軸 û_ij ∝ n_i + n_j
    const rel = [0, 1, 2].map(q => S.v[i] * S.n[i][q] - S.v[j] * S.n[j][q]), sg = Math.sign(dot(rel, u));    // 向きの規約 s_ij = sign((v_i n_i − v_j n_j)·û)
    const dFij = sg * wrap180(S.phi0[i] - S.phi0[j]) * Math.PI / 180; dF[i][j] = dF[j][i] = dFij;            // Δ_F = s (φ₀i − φ₀j)
    let w = 1;                                                                                              // 全量：毎 tick κ F Δ_F
    if (S.mode === 'rate') {                                                                                 // 率：κ F Δ_F |Δd|/W、W = ∫ F⁺ dd（1° 刻み、円周一周）。通過で Σ ≈ Δ_F
      const rin = S.n[i].map(x => x * (S.Th[i] + S.v[i])), rjn = S.n[j].map(x => x * (S.Th[j] + S.v[j]));
      const dd = Math.abs(Math.hypot(rin[0] - rjn[0], rin[1] - rjn[1], rin[2] - rjn[2]) - d), nrm = 0.5 * (1 + c) / Math.sqrt(cii * cjj);
      let W = 0; for (let g = -180; g <= 180; g++) W += Math.max(0, nrm * xcorr(S.co[i], S.co[j], g));
      w = W > 0 ? dd / W : 0; }
    const coef = k * Fij * dFij * w, R = Rmat([0, u[0], u[1], u[2]]);                                        // 生成子 κ F Δ_F w · P_D ⊗ R(û)
    for (let r = 0; r < 4; r++) for (let cc = 0; cc < 4; cc++) { const g = coef * R[r][cc];
      G[4 * i + r][4 * i + cc] += 0.5 * g; G[4 * j + r][4 * j + cc] += 0.5 * g; G[4 * i + r][4 * j + cc] -= 0.5 * g; G[4 * j + r][4 * i + cc] -= 0.5 * g; }
  }
  const U = expmAnti(G), vec = [].concat(...a), out = U.map(row => row.reduce((t, x, cidx) => t + x * vec[cidx], 0));   // 全対同時：(a′) = exp(G) (a)
  const A2 = [], n2 = [], p2 = [];
  for (let i = 0; i < N; i++) { const q = out.slice(4 * i, 4 * i + 4), A = Math.hypot(q[0], q[1], q[2], q[3]), vn = Math.hypot(q[1], q[2], q[3]);
    let nn = S.n[i], ph = S.phi0[i];                                                                         // 読み戻し：A = |a′|、n = Vec/|Vec|（前の n に近い側）、φ₀ = atan2(|Vec|, Sc)
    if (A > 1e-9) { if (vn > 1e-9) { nn = [q[1] / vn, q[2] / vn, q[3] / vn]; const ang = Math.atan2(vn, q[0]) * 180 / Math.PI;   // 振幅が零を通る瞬間は軸を前の値に保つ（絶対しきい値）
        if (dot(nn, S.n[i]) < 0) { nn = nn.map(x => -x); ph = (360 - ang) % 360; } else ph = ang; } else ph = q[0] >= 0 ? 0 : 180; }
    A2.push(A); n2.push(nn); p2.push(ph); }
  const sb = vec.reduce((t, x) => t + x * x, 0), sa = out.reduce((t, x) => t + x * x, 0);
  return {A: A2, n: n2, phi0: p2, F, dF, G, U, sumA2before: sb, sumA2after: sa};
}
const xdesc = () => { const kv = Number(kap.value).toFixed(2); if (!lastX) return `κ = ${kv}、交換 ${xmOf() === 'rate' ? '率' : '全量'}（まだ交換なし）`;
  return `κ = ${kv}、交換 ${xmOf() === 'rate' ? '率' : '全量'}　核 F: ` + PAIRS.map(([i, j], k) => `${PAIRL[k]} ${lastX.F[i][j].toFixed(4)}`).join(', ') + `　Δ_F: ` + PAIRS.map(([i, j], k) => `${PAIRL[k]} ${(lastX.dF[i][j] * 180 / Math.PI).toFixed(1)}°`).join(', ') + `　ΣA² = ${lastX.sumA2after.toFixed(6)}`; };
function selfTest() {   // ?test=1：K1〜K7（算出根拠 §9）。結果は <pre id="ktest">。包絡は原本の式（φ, M, 偶奇はスライダーの値）
  const out = {}, x = [1, 0, 0], y = [0, 1, 0], mx = [-1, 0, 0];
  const mk = (Th, n, v, A, phi0, mode) => ({Th, n, v, A, phi0, co: P.map(u => coeffs(Number(u.ph.value), Number(u.mm.value), u.odd.checked, u.even.checked)), mode: mode || 'full'});
  const maxAbs = M => Math.max(...M.map(r => Math.max(...r.map(Math.abs)))), eye = (i, j) => (i === j ? 1 : 0);
  let S = mk([110, 110, 220, 190], [x, x, mx, mx], [1, -1, 2, -1], [1, 1, 1, 1], [0, 180, 30, 30]), X = exchangeStep(S);   // K3：同形同大・同軸・正面・重なり最大・Δ_F = 180。対 3–4 は軸 −x で 1, 2 から見えず、互いは Δ_F = 0
  out.K1_antisym = maxAbs(X.G.map((r, i) => r.map((g, j) => g + X.G[j][i])));
  out.K1_orth = maxAbs(mmul(X.U.map((r, i) => r.map((_, j) => X.U[j][i])), X.U).map((r, i) => r.map((g, j) => g - eye(i, j))));
  out.K2_sumA2_rel = Math.abs(X.sumA2after - X.sumA2before) / X.sumA2before;
  out.K3_F12 = X.F[0][1]; out.K3_dF12_deg = X.dF[0][1] * 180 / Math.PI; out.K3_phi0_after = X.phi0.slice(0, 2); out.K3_A_after = X.A.slice(0, 2);
  out.K3_swapped = Math.abs(X.phi0[0] - 180) < 1e-9 && Math.abs(X.phi0[1]) < 1e-9;
  S = mk([110, 110, 220, 190], [x, x, y, y], [1, -1, 2, -1], [1, 1, 1, 1], [40, 40, 40, 40]); X = exchangeStep(S);                 // K4：全対 Δ_F = 0 → 恒等
  out.K4_identity = maxAbs(X.U.map((r, i) => r.map((g, j) => g - eye(i, j))));
  const n5 = [Math.cos(5 * Math.PI / 180), Math.sin(5 * Math.PI / 180), 0];
  S = mk([110, 110, 500, 700], [x, n5, y, y], [1, -1, 2, -1], [1, 2, 1, 1], [0, 90, 0, 0]); X = exchangeStep(S);                   // K5：軸が 5° 違う対、A = 1 と 2、Θ 等しい
  out.K5_d12_deg = 2 * 110 * Math.sin(2.5 * Math.PI / 180); out.K5_F12 = X.F[0][1]; out.K5_A_after = X.A.slice(0, 2); out.K5_n1_after = X.n[0]; out.K5_n2_after = X.n[1];
  out.K5_turn_n1_deg = Math.acos(Math.max(-1, Math.min(1, dot(X.n[0], x)))) * 180 / Math.PI; out.K5_turn_n2_deg = Math.acos(Math.max(-1, Math.min(1, dot(X.n[1], n5)))) * 180 / Math.PI;
  const perm = [2, 0, 3, 1], Sp = mk(perm.map(i => S.Th[i]), perm.map(i => S.n[i]), perm.map(i => S.v[i]), perm.map(i => S.A[i]), perm.map(i => S.phi0[i])); Sp.co = perm.map(i => S.co[i]);   // K6：添字の置換
  const Xp = exchangeStep(Sp);
  out.K6_perm_max_diff = Math.max(...perm.map((i, k) => Math.max(Math.abs(Xp.A[k] - X.A[i]), Math.abs(wrap180(Xp.phi0[k] - X.phi0[i])), Math.hypot(...[0, 1, 2].map(q => Xp.n[k][q] - X.n[i][q])))));
  const passage = (mode, sp, nt) => { let T2 = {Th: [50, 170], n: [x, x], v: [sp, -sp], A: [1, 1], phi0: [0, 150], co: [S.co[0], S.co[1]], mode}, mF = 0;   // K7：2 波、同軸 x、正面、Δ_F = 150、相対速度 2sp
    for (let t = 0; t < nt; t++) { T2.Th = [T2.Th[0] + sp, T2.Th[1] - sp]; const Y = exchangeStep(T2); mF = Math.max(mF, Math.abs(Y.F[0][1])); T2 = {Th: T2.Th, n: Y.n, v: T2.v, A: Y.A, phi0: Y.phi0, co: T2.co, mode}; }
    return {phi0: T2.phi0, A: T2.A, maxF: mF}; };
  const r1 = passage('rate', 1, 120), r2 = passage('rate', 0.5, 240), f1 = passage('full', 1, 120), f2 = passage('full', 0.5, 240);   // 同じ道のりを相対速度 2 と 1 で
  out.K7_rate_v1 = r1; out.K7_rate_v05 = r2; out.K7_rate_speed_dep_deg = Math.max(Math.abs(wrap180(r1.phi0[0] - r2.phi0[0])), Math.abs(wrap180(r1.phi0[1] - r2.phi0[1])));
  out.K7_full_v1 = f1; out.K7_full_v05 = f2; out.K7_full_speed_dep_deg = Math.max(Math.abs(wrap180(f1.phi0[0] - f2.phi0[0])), Math.abs(wrap180(f1.phi0[1] - f2.phi0[1])));
  $('ktest').textContent = 'KTEST ' + JSON.stringify(out);
}

function params() {
'''
rep("\nfunction params() {\n", FUNCS)
rep("  P.forEach((u, i) => { const last = cum[i][cum[i].length - 1]; cum[i].push([last[0] + Number(u.vv.value), last[1] + Number(u.dl.value)]); if (cum[i].length > 4000) cum[i].shift(); });\n  draw();\n}\nfunction setAuto(on) {",
    "  P.forEach((u, i) => { const last = cum[i][cum[i].length - 1]; cum[i].push([last[0] + Number(u.vv.value), last[1] + Number(u.dl.value)]); if (cum[i].length > 4000) cum[i].shift(); });\n"
    "  if (Number(kap.value) > 0) {   // 相互作用：自由運動のあとに全対同時の交換\n"
    "    const Th = P.map((u, i) => base[i][0] + cum[i][cum[i].length - 1][0]);\n"
    "    const co = P.map(u => coeffs(Number(u.ph.value), Number(u.mm.value), u.odd.checked, u.even.checked));   // 包絡は原本の式そのまま\n"
    "    const S = {Th, n: nDyn.map(x => x.slice()), v: P.map(u => Number(u.vv.value)), A: aDyn.slice(), phi0: P.map(u => u.phase0), co, mode: xmOf()};\n"
    "    const X = exchangeStep(S); lastX = X;\n"
    "    const dph = P.map((u, i) => wrap180(X.phi0[i] - u.phase0));   // 交換による内部時計の回転（度）\n"
    "    P.forEach((u, i) => { aDyn[i] = X.A[i]; nDyn[i] = X.n[i]; u.phase0 = X.phi0[i]; cum[i][cum[i].length - 1][1] += dph[i]; });   // 時計の累積にも足す（読みを整合させる）\n"
    "    PAIRS.forEach(([i, j], k) => { hist[k][hist[k].length - 1][1] += dph[i] - dph[j]; });\n"
    "    ihist.push({t: ticks, F: X.F, dF: X.dF, s2: X.sumA2after, A: X.A.slice(), n: X.n.map(v => v.slice()), phi0: X.phi0.slice()}); if (ihist.length > 4000) ihist.shift();\n"
    "  }\n  draw();\n}\nfunction setAuto(on) {")
rep("const resetAll = () => { ticks = 0; orb = [[], [], [], []]; orbTick = -1; hist = PAIRS.map(() => [[0, 0]]); cum = [0, 1, 2, 3].map(() => [[0, 0]]); win0 = null; for (const w of P) w.phase0 = 0; base = P.map(u => [Number(u.th.value), 0]); draw(); };",
    "const resetAll = () => { ticks = 0; orb = [[], [], [], []]; orbTick = -1; hist = PAIRS.map(() => [[0, 0]]); cum = [0, 1, 2, 3].map(() => [[0, 0]]); win0 = null; for (const w of P) w.phase0 = 0; base = P.map(u => [Number(u.th.value), 0]); nDyn = P.map(u => nvec(u)); aDyn = P.map(u => Number(u.aa.value)); ihist = []; lastX = null; draw(); };\n"
    "kap.addEventListener('input', resetAll);   // κ を変えたら最初から（resetAll の定義のあと）\n"
    "for (const r of document.querySelectorAll('input[name=xm]')) r.addEventListener('change', resetAll);")
rep("  if (q.has('s')) ss.value = q.get('s');   // s = 球の倍率（撮影用）\n",
    "  if (q.has('s')) ss.value = q.get('s');   // s = 球の倍率（撮影用）\n"
    "  if (q.has('kap')) kap.value = q.get('kap');   // kap = 相互作用の強さ（0 で相互作用なし版と同一）\n"
    "  if (q.has('xm')) document.querySelector(`input[name=xm][value=\"${q.get('xm')}\"]`).checked = true;   // xm = rate | full\n")
rep("  base = P.map(u => [Number(u.th.value), u.phase0]);\n  draw();\n  setAuto(q.get('auto') !== '0');\n",
    "  base = P.map(u => [Number(u.th.value), u.phase0]); nDyn = P.map(u => nvec(u)); aDyn = P.map(u => Number(u.aa.value));\n  draw();\n  if (q.get('test') === '1') selfTest();\n  setAuto(q.get('auto') !== '0');\n")
assert "coeffsCos2" not in s and 'name="env"' not in s and s.count('id="dl2" type="range" min="-30" max="30" step="1" value="0"') == 1 and s.count('id="dl4" type="range" min="-30" max="30" step="1" value="0"') == 1
io.open(p, "w", encoding="utf-8").write(s); print("ok")
