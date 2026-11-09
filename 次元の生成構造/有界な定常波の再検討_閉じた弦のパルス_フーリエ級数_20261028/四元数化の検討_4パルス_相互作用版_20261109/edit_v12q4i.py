# -*- coding: utf-8 -*-
# index_v11q4.html → index_v12q4i.html：相互作用版 v0.6（2026-11-10）。
# 木原氏 11/10：相互作用は繰り返し起こる。一回で反転する写像を毎 tick 掛ければ反転を繰り返す。相互作用とは、自己無撞着のそれ自身が U^n = I になる反応。
# 一歩は 1 tick。重なっている（台が交わる：|Θ_i − Θ_j| < (φ_i + φ_j)/2）tick に、対の相対運動量を、一歩ぶんの位相のずれ δ_i − δ_j（内部時計の差の 1 tick 分）だけ回す。
# 回転軸は対の関係の Vec 部 n_i × n_j（添字の入替えで角も軸も反転して写像は同じ）。質量相当 m = A²、p = m v n、重心系で q を回し Σp は厳密に保存。
# 変わるのは向き n と速さ v だけ。A・φ₀・位相 Θ の進み方・包絡は触らない。位置 r は向きの積分（右パネルの点）。
# 核 F と Δ_F は表示だけ（写像には使わない。v0.5 の α = κ F Δ_F w は有限位数を壊す誤り）。
import io
p = "index_v12q4i.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new, n=1):
    global s
    c = s.count(old); assert c == n, (old[:90], c, n); s = s.replace(old, new)
rep("<title>閉じた弦の四元数パルス 4 波（index_v8q3.html の 4 パルス版）</title>", "<title>閉じた弦の四元数パルス 4 波 相互作用版（index_v11q4.html に、重なっている tick の有限位数回転）</title>")
rep('  <div class="sl"><span>s</span><input id="ss" type="range" min="1" max="100" step="1" value="10"><span class="v" id="ssv">10</span></div>\n',
    '  <div class="sl"><span>s</span><input id="ss" type="range" min="1" max="100" step="1" value="10"><span class="v" id="ssv">10</span></div>\n'
    '  <label class="ck"><input type="checkbox" id="kap" checked>相互作用（重なっている tick に δ_i − δ_j だけ回す）</label>\n')
rep("<script>\nconst $ = id => document.getElementById(id);", '<pre id="ktest" style="font-size:12px;white-space:pre-wrap;margin:4px 8px"></pre>\n<script>\nconst $ = id => document.getElementById(id);')
rep("const rr = $('rr'), ss = $('ss'), btn = $('auto');   // ss = 球の倍率 s（1〜100、初期 10）\n",
    "const rr = $('rr'), ss = $('ss'), btn = $('auto'), kap = $('kap');   // ss = 球の倍率 s（1〜100、初期 10）、kap = 相互作用の有無（外すと相互作用なし版と同一）\n"
    "let nDyn = null, vDyn = null, rDyn = null, ihist = [], lastX = null;   // 相互作用で変わる向き n_i と速さ v_i（スライダーは初期値）、位置 r_i（向きの積分、図の点）、記録\n"
    "const kapOn = () => kap.checked ? 1 : 0;\n"
    "const fmtv = v => Number.isInteger(v) ? String(v) : v.toFixed(2);   // 速さの表示：整数はそのまま（相互作用なしで原本と同じ文字列）、変われば小数 2 桁\n")
rep("  const pulses = P.map(u => {\n", "  const pulses = P.map((u, i) => {\n")
rep("    const theta = Number(u.th.value), phi = Number(u.ph.value), v = Number(u.vv.value), M = Number(u.mm.value), delta = Number(u.dl.value), phase0 = u.phase0;",
    "    const theta = Number(u.th.value), phi = Number(u.ph.value), v = vDyn ? vDyn[i] : Number(u.vv.value), M = Number(u.mm.value), delta = Number(u.dl.value), phase0 = u.phase0;   // v は相互作用で変わる（vDyn）")
rep("    return {theta, phi, v, M, delta, phase0, co, A, nv: nvec(u), R: Math.sqrt(R2), R2,", "    return {theta, phi, v, M, delta, phase0, co, A, nv: nDyn ? nDyn[i] : nvec(u), R: Math.sqrt(R2), R2,")
rep("u.vvv.textContent = q.v + '°/t';", "u.vvv.textContent = fmtv(q.v) + '°/t';")
rep("δ = ${u.delta}°/t, M = ${u.M}（奇", "δ = ${u.delta}°/t, φ₀ = ${u.phase0.toFixed(1)}°, M = ${u.M}（奇")
rep("v = ${u.v}°/t, δ", "v = ${fmtv(u.v)}°/t, δ")
rep("  $('h1').textContent = `${desc}　　和：R = ", "  $('h1').textContent = `${desc}　　${xdesc()}　　和：R = ")
rep("  const orbPts = p.pulses.map((u, i) => { const nv = u.nv; return [nv[0] * Th[i], nv[1] * Th[i], nv[2] * Th[i]]; });   // 波 i の点 = Θ_i n_i = (i_i Θ_i, j_i Θ_i, k_i Θ_i)",
    "  const orbPts = p.pulses.map((u, i) => rDyn ? rDyn[i].slice() : [u.nv[0] * Th[i], u.nv[1] * Th[i], u.nv[2] * Th[i]]);   // 波 i の点 = 位置 r_i（向きの積分 r ← r + v n。n が変わらなければ Θ_i n_i）")
FUNCS = r'''
// ---------- 相互作用（算出根拠_相互作用版 v0.6 §2〜§3）。重なっている tick に、対の相対運動量を δ_i − δ_j だけ n_i × n_j まわりに回す。A・φ₀・Θ は触らない ----------
const v3add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]], v3sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], v3mul = (a, k) => [a[0] * k, a[1] * k, a[2] * k], v3len = a => Math.hypot(a[0], a[1], a[2]);
const rotAxis = (x, e, a) => { const c = Math.cos(a), s = Math.sin(a), k = dot(e, x), cx = crs(e, x); return [x[0] * c + cx[0] * s + e[0] * k * (1 - c), x[1] * c + cx[1] * s + e[1] * k * (1 - c), x[2] * c + cx[2] * s + e[2] * k * (1 - c)]; };   // Rodrigues：x を単位軸 e まわりに a 回す
function xcorr(coI, coJ, dDeg) {   // 包絡どうしの重なり ⟨S_i S_j⟩(d) = c₀i c₀j + ½ Σ_m c_i,m c_j,m cos(m d)（表示用。写像には使わない）
  const mapJ = new Map(); coJ.ms.forEach((m, k) => mapJ.set(m, coJ.cs[k])); const d = dDeg * Math.PI / 180; let t = coI.c0 * coJ.c0;
  coI.ms.forEach((m, k) => { const cj = mapJ.get(m); if (cj !== undefined) t += 0.5 * coI.cs[k] * cj * Math.cos(m * d); }); return t;
}
function exchangeStep(S) {   // S = {Th, n, v, A, phi（幅 φ_i）, delta（δ_i）, co, k}。純関数。返り値 {n, v, d, ov, F, alpha, E, undef, dp, Pbefore, Pafter, Kbefore, Kafter}
  const N = S.A.length, k = S.k === undefined ? kapOn() : S.k, m = S.A.map(A => A * A);                      // 質量相当 m_i = A_i²
  const pm = S.v.map((v, i) => v3mul(S.n[i], m[i] * v));                                                   // 運動量 p_i = m_i v_i n_i
  const Z = () => S.A.map(() => new Array(N).fill(0));
  const d = Z(), ov = Z(), F = Z(), alpha = Z(), E = S.A.map(() => new Array(N).fill(null)), undef = Z(), dp = S.A.map(() => [0, 0, 0]);
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
    const dij = Math.abs(S.Th[i] - S.Th[j]); d[i][j] = d[j][i] = dij;                                       // 位相の差 d_ij = |Θ_i − Θ_j|（度、巻数込み）
    const o = dij < 0.5 * (S.phi[i] + S.phi[j]) ? 1 : 0; ov[i][j] = ov[j][i] = o;                            // 重なり：台が交わる d < (φ_i + φ_j)/2
    if (S.co) { const cii = xcorr(S.co[i], S.co[i], 0), cjj = xcorr(S.co[j], S.co[j], 0); F[i][j] = F[j][i] = (cii > 0 && cjj > 0) ? 0.5 * (1 + dot(S.n[i], S.n[j])) * xcorr(S.co[i], S.co[j], dij) / Math.sqrt(cii * cjj) : 0; }   // 表示用
    if (!o || !k) continue;
    const a = k * (S.delta[i] - S.delta[j]) * Math.PI / 180; alpha[i][j] = a; alpha[j][i] = -a;              // 一歩の角 α_ij = δ_i − δ_j（内部位相差の 1 tick 分）。入替えで反転
    if (Math.abs(a) < 1e-300) continue;
    const cx = crs(S.n[i], S.n[j]), cl = v3len(cx);
    if (cl < 1e-12) { undef[i][j] = undef[j][i] = 1; continue; }                                             // 軸 n_i × n_j が零（向きが平行）なら回せない。記録して回さない
    const e = v3mul(cx, 1 / cl); E[i][j] = e; E[j][i] = v3mul(e, -1);                                        // 軸 e_ij = (n_i × n_j)/|n_i × n_j|。入替えで反転（角も反転するので写像は同じ）
    const M = m[i] + m[j], Ptot = v3add(pm[i], pm[j]), q = v3mul(v3sub(v3mul(pm[i], m[j]), v3mul(pm[j], m[i])), 1 / M);   // 相対運動量 q = (m_j p_i − m_i p_j)/M
    if (v3len(q) < 1e-12) continue;
    const q2 = rotAxis(q, e, a);                                                                             // 重心系で q を e まわりに α 回す（P は不変、|q| 不変）
    const pi2 = v3add(v3mul(Ptot, m[i] / M), q2), pj2 = v3sub(v3mul(Ptot, m[j] / M), q2);                     // p_i′ = (m_i/M) P + q′、p_j′ = (m_j/M) P − q′
    dp[i] = v3add(dp[i], v3sub(pi2, pm[i])); dp[j] = v3add(dp[j], v3sub(pj2, pm[j]));                         // 全対の変化を足す（順序なし。対ごとに Σp 保存）
  }
  const n2 = [], v2 = [];
  for (let i = 0; i < N; i++) { const pn = v3add(pm[i], dp[i]), un = v3mul(pn, 1 / m[i]), sp = v3len(un);     // u_i′ = p_i′/m_i
    if (sp < 1e-12) { n2.push(S.n[i]); v2.push(S.v[i]); continue; }
    const sg = S.v[i] < 0 ? -1 : 1; n2.push(v3mul(un, sg / sp)); v2.push(sg * sp); }                          // v の符号の規約は保つ：n_i′ = sign(v_i) u′/|u′|、v_i′ = sign(v_i)|u′|。反転なら n → −n
  const K = vv => vv.reduce((t, v, i) => t + 0.5 * m[i] * v * v, 0), Psum = pp => pp.reduce((t, p) => v3add(t, p), [0, 0, 0]);
  return {n: n2, v: v2, d, ov, F, alpha, E, undef, dp, Pbefore: Psum(pm), Pafter: Psum(n2.map((n, i) => v3mul(n, m[i] * v2[i]))), Kbefore: K(S.v), Kafter: K(v2)};
}
const xdesc = () => { if (!lastX) return `相互作用 ${kapOn() ? 'あり' : 'なし'}（まだ交換なし）`;
  const P = lastX.Pafter; return `相互作用 ${kapOn() ? 'あり' : 'なし'}　d: ` + PAIRS.map(([i, j], k) => `${PAIRL[k]} ${lastX.d[i][j].toFixed(1)}°${lastX.ov[i][j] ? '◯' : ''}`).join(', ') + `　α: ` + PAIRS.map(([i, j], k) => `${PAIRL[k]} ${(lastX.alpha[i][j] * 180 / Math.PI).toFixed(2)}°${lastX.undef[i][j] ? '(軸なし)' : ''}`).join(', ') + `　核 F（表示）: ` + PAIRS.map(([i, j], k) => `${PAIRL[k]} ${lastX.F[i][j].toFixed(3)}`).join(', ') + `　Σ m v n = (${P.map(x => x.toFixed(4)).join(', ')})　Σ ½ m v² = ${lastX.Kafter.toFixed(6)}`; };
function selfTest() {   // ?test=1：K1〜K7（算出根拠 v0.6 §9）。結果は <pre id="ktest">
  const out = {}, x = [1, 0, 0], ex = [1, 0, 0.05], y = [0, 1, 0], ey = [0, 1, 0.05], un = v => v3mul(v, 1 / v3len(v));
  const mk = (Th, n, v, A, delta, phi) => ({Th, n: n.map(un), v, A, delta, phi: phi || Th.map(() => 20), k: 1});
  const angle = (a, b) => Math.acos(Math.max(-1, Math.min(1, dot(a, b) / (v3len(a) * v3len(b))))) * 180 / Math.PI;
  // K1：4 波・2 対が同時に重なる 1 tick。Σp と Σ ½ m v² の保存
  let S = mk([110, 111, 280, 281], [x, ex, y, ey], [1, -1, 1, -1], [0.5, 1, 1.5, 2], [0, 2, 0, 2]), X = exchangeStep(S);
  out.K1_P_conserved = v3len(v3sub(X.Pafter, X.Pbefore)); out.K1_K_rel = Math.abs(X.Kafter - X.Kbefore) / X.Kbefore;
  // K3：一歩の角と門。重なり（d = 1 < 20）で α₁₂ = δ₁ − δ₂ = −2°、重なりの外（d = 25）で 0
  out.K3_alpha12_deg_overlap = X.alpha[0][1] * 180 / Math.PI; out.K3_ov = [X.ov[0][1], X.ov[2][3], X.ov[0][2]];
  S = mk([100, 125, 280, 281], [x, ex, y, ey], [1, -1, 1, -1], [0.5, 1, 1.5, 2], [0, 2, 0, 2]); X = exchangeStep(S); out.K3_alpha12_deg_apart = X.alpha[0][1] * 180 / Math.PI; out.K3_n1_apart_unchanged = v3len(v3sub(X.n[0], S.n[0]));
  // K4：δ が全部同じ → 何も変わらない
  S = mk([110, 111, 280, 281], [x, ex, y, ey], [1, -1, 1, -1], [0.5, 1, 1.5, 2], [2, 2, 2, 2]); X = exchangeStep(S);
  out.K4_identity = Math.max(...X.n.map((n, i) => v3len(v3sub(n, S.n[i]))), ...X.v.map((v, i) => Math.abs(v - S.v[i])));
  // K5：通過と閉包。2 波、Θ 61/180、v = ±1 → d = 19, 17, …, 1, 1, …, 19 の 20 tick が重なる（n = 20）。nδ = 360（δ₂ = 18）で U^n = I、180（δ₂ = 9）で反転か
  const passage = (dl2, A, sp, nt, Th0, nn) => { let T = {Th: Th0 ? Th0.slice() : [61, 180], n: (nn || [x, ex]).map(un), v: [sp, -sp], A, delta: [0, dl2], phi: [20, 20], k: 1};
    const u0 = T.n.map((n, i) => v3mul(n, T.v[i])); const tr = []; let cum = 0, nov = 0, flips = 0, lastE = null;
    for (let t = 0; t < nt; t++) { T.Th = [T.Th[0] + T.v[0], T.Th[1] + T.v[1]]; const Y = exchangeStep(T);
      if (Y.ov[0][1]) { nov++; cum += Y.alpha[0][1] * 180 / Math.PI; if (Y.E[0][1]) { if (lastE && dot(lastE, Y.E[0][1]) < 0) flips++; lastE = Y.E[0][1]; } tr.push([t, +Math.abs(T.Th[0] - T.Th[1]).toFixed(1), +(Y.alpha[0][1] * 180 / Math.PI).toFixed(2), Y.E[0][1] ? +Y.E[0][1][1].toFixed(3) : null, Y.undef[0][1], +Y.n[0][0].toFixed(3), +Y.n[0][2].toFixed(3), +Y.v[0].toFixed(3), +Y.n[1][0].toFixed(3), +Y.n[1][2].toFixed(3), +Y.v[1].toFixed(3)]); }
      T = {Th: T.Th, n: Y.n, v: Y.v, A: T.A, delta: T.delta, phi: T.phi, k: 1}; }
    const u1 = T.n.map((n, i) => v3mul(n, T.v[i]));
    return {n_overlap: nov, sum_alpha_deg: +cum.toFixed(2), axis_flips: flips, turn1_deg: +angle(u1[0], u0[0]).toFixed(2), turn2_deg: +angle(u1[1], u0[1]).toFixed(2), closure_err: +Math.max(v3len(v3sub(u1[0], u0[0])), v3len(v3sub(u1[1], u0[1]))).toFixed(4), reversal_err: +Math.max(v3len(v3add(u1[0], u0[0])), v3len(v3add(u1[1], u0[1]))).toFixed(4), v_final: T.v.map(v => +v.toFixed(4)), n1_final: T.n[0].map(c => +c.toFixed(4)), n2_final: T.n[1].map(c => +c.toFixed(4)), trace: tr}; };
  out.K5a_closure_dl18_equal = passage(18, [1, 1], 1, 130); out.K5b_half_dl9_equal = passage(9, [1, 1], 1, 130); out.K5c_dl2_equal = passage(2, [1, 1], 1, 130);
  out.K5d_closure_dl18_mass1_4 = passage(18, [1, 2], 1, 130); out.K5e_half_dl9_mass1_4 = passage(9, [1, 2], 1, 130);
  // K6：添字の置換（4 波、2 対重なり）
  S = mk([110, 111, 280, 281], [x, ex, y, ey], [1, -1, 1, -1], [0.5, 1, 1.5, 2], [0, 2, 0, 2]); X = exchangeStep(S);
  const perm = [2, 0, 3, 1], Sp = mk(perm.map(i => S.Th[i]), perm.map(i => S.n[i]), perm.map(i => S.v[i]), perm.map(i => S.A[i]), perm.map(i => S.delta[i])), Xp = exchangeStep(Sp);
  out.K6_perm_max_diff = Math.max(...perm.map((i, k) => Math.max(v3len(v3sub(Xp.n[k], X.n[i])), Math.abs(Xp.v[k] - X.v[i]))));
  // K7：δ₂ の走査（0〜20、0.5 刻み）。同質量と質量 1:4、相対速度 2（n = 20）と 1（Θ 60.5/180 → n = 40）。終状態の向きの回転と閉包・反転の誤差、軸の反転回数
  const scan = (A, sp, Th0) => { const rows = []; for (let dl = 0; dl <= 20; dl += 0.5) { const r = passage(dl, A, sp, Math.round(130 / sp), Th0); rows.push([dl, r.n_overlap, r.sum_alpha_deg, r.axis_flips, r.turn1_deg, r.turn2_deg, r.closure_err, r.reversal_err]); } return rows; };
  out.K7_scan_cols = ["dl2", "n_overlap", "sum_alpha", "axis_flips", "turn1", "turn2", "closure_err", "reversal_err"];
  out.K7_scan_equal_v2 = scan([1, 1], 1); out.K7_scan_mass1_4_v2 = scan([1, 2], 1); out.K7_scan_equal_v1 = scan([1, 1], 0.5, [60.5, 180]);
  $('ktest').textContent = 'KTEST ' + JSON.stringify(out);
}

function params() {
'''
rep("\nfunction params() {\n", FUNCS)
rep("  for (const u of P) {\n    let t = Number(u.th.value) + Number(u.vv.value); t = ((t % 360) + 360) % 360; u.th.value = String(t);   // 位置 θ ← θ + v（負も折り返す）\n",
    "  P.forEach((u, i) => {\n    const vi = vDyn ? vDyn[i] : Number(u.vv.value);\n    let t = Number(u.th.value) + vi; t = ((t % 360) + 360) % 360; u.th.value = String(t);   // 位置 θ ← θ + v（負も折り返す。v は相互作用で変わる）\n    if (rDyn) rDyn[i] = rDyn[i].map((x, q) => x + vi * nDyn[i][q]);                           // 位置 r ← r + v n（向きの積分）\n")
rep("    u.phase0 = ((u.phase0 + Number(u.dl.value)) % 360 + 360) % 360;                                    // 内部位相 φ₀ ← φ₀ + δ（位置とは独立）\n  }\n",
    "    u.phase0 = ((u.phase0 + Number(u.dl.value)) % 360 + 360) % 360;                                    // 内部位相 φ₀ ← φ₀ + δ（位置とは独立）\n  });\n")
rep("    hist[k].push([last[0] + Number(P[i].vv.value) - Number(P[j].vv.value), last[1] + Number(P[i].dl.value) - Number(P[j].dl.value)]);",
    "    hist[k].push([last[0] + vDyn[i] - vDyn[j], last[1] + Number(P[i].dl.value) - Number(P[j].dl.value)]);")
rep("  P.forEach((u, i) => { const last = cum[i][cum[i].length - 1]; cum[i].push([last[0] + Number(u.vv.value), last[1] + Number(u.dl.value)]); if (cum[i].length > 4000) cum[i].shift(); });\n  draw();\n}\nfunction setAuto(on) {",
    "  P.forEach((u, i) => { const last = cum[i][cum[i].length - 1]; cum[i].push([last[0] + vDyn[i], last[1] + Number(u.dl.value)]); if (cum[i].length > 4000) cum[i].shift(); });\n"
    "  if (kapOn()) {   // 相互作用：自由運動のあとに、重なっている対を一歩だけ回す。変わるのは n_i と v_i だけ（A・φ₀・位相 Θ は変えない）\n"
    "    const Th = P.map((u, i) => base[i][0] + cum[i][cum[i].length - 1][0]);\n"
    "    const co = P.map(u => coeffs(Number(u.ph.value), Number(u.mm.value), u.odd.checked, u.even.checked));   // 核 F の表示用（包絡は原本の式そのまま）\n"
    "    const S = {Th, n: nDyn.map(x => x.slice()), v: vDyn.slice(), A: P.map(u => Number(u.aa.value)), phi: P.map(u => Number(u.ph.value)), delta: P.map(u => Number(u.dl.value)), co};\n"
    "    const X = exchangeStep(S); lastX = X;\n"
    "    P.forEach((u, i) => { nDyn[i] = X.n[i]; vDyn[i] = X.v[i]; });\n"
    "    ihist.push({t: ticks, d: X.d, ov: X.ov, alpha: X.alpha, n: X.n.map(v => v.slice()), v: X.v.slice(), P: X.Pafter, K: X.Kafter}); if (ihist.length > 4000) ihist.shift();\n"
    "  }\n  draw();\n}\nfunction setAuto(on) {")
rep("const resetAll = () => { ticks = 0; orb = [[], [], [], []]; orbTick = -1; hist = PAIRS.map(() => [[0, 0]]); cum = [0, 1, 2, 3].map(() => [[0, 0]]); win0 = null; for (const w of P) w.phase0 = 0; base = P.map(u => [Number(u.th.value), 0]); draw(); };",
    "const initDyn = () => { nDyn = P.map(u => nvec(u)); vDyn = P.map(u => Number(u.vv.value)); rDyn = P.map((u, i) => nDyn[i].map(x => x * Number(u.th.value))); ihist = []; lastX = null; };   // 動的状態の初期化：n・v はスライダー、r = θ n\n"
    "const resetAll = () => { ticks = 0; orb = [[], [], [], []]; orbTick = -1; hist = PAIRS.map(() => [[0, 0]]); cum = [0, 1, 2, 3].map(() => [[0, 0]]); win0 = null; for (const w of P) w.phase0 = 0; base = P.map(u => [Number(u.th.value), 0]); initDyn(); draw(); };\n"
    "kap.addEventListener('change', resetAll);   // 相互作用の有無を変えたら最初から")
rep("  if (q.has('s')) ss.value = q.get('s');   // s = 球の倍率（撮影用）\n",
    "  if (q.has('s')) ss.value = q.get('s');   // s = 球の倍率（撮影用）\n"
    "  if (q.has('kap')) kap.checked = q.get('kap') !== '0';   // kap=0 で相互作用なし（原本と同一）\n")
rep("  base = P.map(u => [Number(u.th.value), u.phase0]);\n  draw();\n  setAuto(q.get('auto') !== '0');\n",
    "  base = P.map(u => [Number(u.th.value), u.phase0]); initDyn();\n  draw();\n  if (q.get('test') === '1') selfTest();\n  setAuto(q.get('auto') !== '0');\n")
assert "aDyn" not in s and "expmAnti" not in s and "name=xm" not in s and "wrap180(S.phi0" not in s
io.open(p, "w", encoding="utf-8").write(s); print("ok")
