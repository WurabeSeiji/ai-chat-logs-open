# -*- coding: utf-8 -*-
# index_v11q4.html → index_v12q4i.html：相互作用版 v0.5（2026-11-10 未明）。
# 相互作用はゲージ変換＝虚数単位 n（位相の傾斜の向き）を回す。位相 φ₀ と振幅 A は触らない（木原氏 11/9〜10：位相をいくら回しても向きは変わらない、A は損失がないので動かす理由がない、質量相当は A²）。
# 自由運動は従来どおり（Θ は v で、φ₀ は δ で進む）。位置 r は向きの積分 r ← r + v n（n が変われば今いる点から新しい向きへ）。
# 回す角 α_ij = κ F_ij Δ_F w（7/10 の二チャネル写像の角を向きの回転角と読む）。散乱面の法線 e = unit(q × b)（q 相対運動量、b 衝突径数）。
# 重心系で q を e まわりに α 回し、P = p_i + p_j を保つ。p_i = m_i v_i n_i、m_i = A_i²。全対の運動量変化を足して一度に適用（順序なし、Σ p 保存）。
import io
p = "index_v12q4i.html"; s = io.open(p, encoding="utf-8").read()
def rep(old, new, n=1):
    global s
    c = s.count(old); assert c == n, (old[:90], c, n); s = s.replace(old, new)
rep("<title>閉じた弦の四元数パルス 4 波（index_v8q3.html の 4 パルス版）</title>", "<title>閉じた弦の四元数パルス 4 波 相互作用版（index_v11q4.html に重なりの核とゲージ回転）</title>")
rep('  <div class="sl"><span>s</span><input id="ss" type="range" min="1" max="100" step="1" value="10"><span class="v" id="ssv">10</span></div>\n',
    '  <div class="sl"><span>s</span><input id="ss" type="range" min="1" max="100" step="1" value="10"><span class="v" id="ssv">10</span></div>\n'
    '  <div class="sl"><span>κ 相互作用</span><input id="kap" type="range" min="0" max="1" step="0.01" value="1"><span class="v" id="kapv">1.00</span></div>\n'
    '  <span class="grp">交換 <label class="ck"><input type="radio" name="xm" value="rate" checked>率（通過で合計 Δ_F）</label><label class="ck"><input type="radio" name="xm" value="full">全量（毎 tick F Δ_F）</label></span>\n')
rep("<script>\nconst $ = id => document.getElementById(id);", '<pre id="ktest" style="font-size:12px;white-space:pre-wrap;margin:4px 8px"></pre>\n<script>\nconst $ = id => document.getElementById(id);')
rep("const rr = $('rr'), ss = $('ss'), btn = $('auto');   // ss = 球の倍率 s（1〜100、初期 10）\n",
    "const rr = $('rr'), ss = $('ss'), btn = $('auto'), kap = $('kap');   // ss = 球の倍率 s（1〜100、初期 10）、kap = 相互作用の強さ κ（0〜1、0 で相互作用なし版と同一）\n"
    "let nDyn = null, vDyn = null, rDyn = null, ihist = [], lastX = null;   // 相互作用で変わる向き n_i と速さ v_i（スライダーは初期値）、位置 r_i（向きの積分、図の点）、記録\n"
    "const xmOf = () => document.querySelector('input[name=xm]:checked').value;   // 交換の方式：rate（率）| full（全量）\n"
    "const fmtv = v => Number.isInteger(v) ? String(v) : v.toFixed(2);   // 速さの表示：整数はそのまま（κ = 0 で原本と同じ文字列）、変われば小数 2 桁\n")
rep("  const pulses = P.map(u => {\n", "  const pulses = P.map((u, i) => {\n")
rep("    const theta = Number(u.th.value), phi = Number(u.ph.value), v = Number(u.vv.value), M = Number(u.mm.value), delta = Number(u.dl.value), phase0 = u.phase0;",
    "    const theta = Number(u.th.value), phi = Number(u.ph.value), v = vDyn ? vDyn[i] : Number(u.vv.value), M = Number(u.mm.value), delta = Number(u.dl.value), phase0 = u.phase0;   // v は相互作用で変わる（vDyn）")
rep("    return {theta, phi, v, M, delta, phase0, co, A, nv: nvec(u), R: Math.sqrt(R2), R2,", "    return {theta, phi, v, M, delta, phase0, co, A, nv: nDyn ? nDyn[i] : nvec(u), R: Math.sqrt(R2), R2,")
rep("  $('rrv').textContent = p.k; $('ssv').textContent = ss.value;\n", "  $('rrv').textContent = p.k; $('ssv').textContent = ss.value; $('kapv').textContent = Number(kap.value).toFixed(2);\n")
rep("u.vvv.textContent = q.v + '°/t';", "u.vvv.textContent = fmtv(q.v) + '°/t';")
rep("δ = ${u.delta}°/t, M = ${u.M}（奇", "δ = ${u.delta}°/t, φ₀ = ${u.phase0.toFixed(1)}°, M = ${u.M}（奇")
rep("v = ${u.v}°/t, δ", "v = ${fmtv(u.v)}°/t, δ")
rep("  $('h1').textContent = `${desc}　　和：R = ", "  $('h1').textContent = `${desc}　　${xdesc()}　　和：R = ")
# 右パネルの点：Θ_i n_i → 位置 r_i（向きの積分）
rep("  const orbPts = p.pulses.map((u, i) => { const nv = u.nv; return [nv[0] * Th[i], nv[1] * Th[i], nv[2] * Th[i]]; });   // 波 i の点 = Θ_i n_i = (i_i Θ_i, j_i Θ_i, k_i Θ_i)",
    "  const orbPts = p.pulses.map((u, i) => rDyn ? rDyn[i].slice() : [u.nv[0] * Th[i], u.nv[1] * Th[i], u.nv[2] * Th[i]]);   // 波 i の点 = 位置 r_i（向きの積分 r ← r + v n。n が変わらなければ Θ_i n_i）")
FUNCS = r'''
// ---------- 相互作用（算出根拠_相互作用版 v0.5 §2〜§4）。ゲージ回転：向き n_i（と速さ v_i）だけを変える。A と φ₀ は触らない ----------
const v3add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]], v3sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], v3mul = (a, k) => [a[0] * k, a[1] * k, a[2] * k], v3len = a => Math.hypot(a[0], a[1], a[2]);
const rotAxis = (x, e, a) => { const c = Math.cos(a), s = Math.sin(a), k = dot(e, x), cx = crs(e, x); return [x[0] * c + cx[0] * s + e[0] * k * (1 - c), x[1] * c + cx[1] * s + e[1] * k * (1 - c), x[2] * c + cx[2] * s + e[2] * k * (1 - c)]; };   // Rodrigues：x を単位軸 e まわりに a 回す
function xcorr(coI, coJ, dDeg) {   // 原本の包絡どうしの重なり ⟨S_i S_j⟩(d) = c₀i c₀j + ½ Σ_m c_i,m c_j,m cos(m d)（共有する m、直流を含む）
  const mapJ = new Map(); coJ.ms.forEach((m, k) => mapJ.set(m, coJ.cs[k])); const d = dDeg * Math.PI / 180; let t = coI.c0 * coJ.c0;
  coI.ms.forEach((m, k) => { const cj = mapJ.get(m); if (cj !== undefined) t += 0.5 * coI.cs[k] * cj * Math.cos(m * d); }); return t;
}
function exchangeStep(S) {   // S = {Th, n, v, A, phi0（度）, r, co, mode}。純関数。返り値 {n, v, F, dF, alpha, dp, Pbefore, Pafter, Kbefore, Kafter}
  const N = S.A.length, k = Number(kap.value), m = S.A.map(A => A * A);                                     // 質量相当 m_i = A_i²（パワー。木原氏 11/10）
  const pm = S.v.map((v, i) => v3mul(S.n[i], m[i] * v));                                                   // 運動量 p_i = m_i v_i n_i
  const F = S.A.map(() => new Array(N).fill(0)), dF = S.A.map(() => new Array(N).fill(0)), alpha = S.A.map(() => new Array(N).fill(0)), dp = S.A.map(() => [0, 0, 0]), E = S.A.map(() => new Array(N).fill(null));
  for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
    const d = Math.abs(S.Th[i] - S.Th[j]);                                                                  // 位相の差 d_ij = |Θ_i − Θ_j|（度、巻数込み）。位置は位相、向きは別
    const c = dot(S.n[i], S.n[j]), cii = xcorr(S.co[i], S.co[i], 0), cjj = xcorr(S.co[j], S.co[j], 0);
    const Ft = (cii > 0 && cjj > 0) ? xcorr(S.co[i], S.co[j], d) / Math.sqrt(cii * cjj) : 0;
    const Fij = 0.5 * (1 + c) * Ft; F[i][j] = F[j][i] = Fij;                                                 // 核 F_ij = (1 + n_i·n_j)/2 · ⟨S_i S_j⟩(d)/√(⟨S_i²⟩⟨S_j²⟩)
    const us = v3add(S.n[i], S.n[j]), L = v3len(us);
    if (L < 1e-9 || Math.abs(Fij) < 1e-300) continue;
    const u = v3mul(us, 1 / L);                                                                              // 共有軸 û_ij ∝ n_i + n_j（向きの規約にだけ使う）
    const rel = v3sub(v3mul(S.n[i], S.v[i]), v3mul(S.n[j], S.v[j])), sg = Math.sign(dot(rel, u));             // 向きの規約 s_ij = sign((v_i n_i − v_j n_j)·û)
    const dFij = sg * wrap180(S.phi0[i] - S.phi0[j]) * Math.PI / 180; dF[i][j] = dF[j][i] = dFij;            // Δ_F = s (φ₀i − φ₀j)
    let w = 1;                                                                                              // 全量：毎 tick κ F Δ_F
    if (S.mode === 'rate') {                                                                                 // 率：κ F Δ_F |Δd|/W、Δd = |v_i − v_j|、W = ∫ F⁺ dd（1° 刻み）。通過で合計 Δ_F
      const dd = Math.abs(S.v[i] - S.v[j]), nrm = 0.5 * (1 + c) / Math.sqrt(cii * cjj);
      let W = 0; for (let g = -180; g <= 180; g++) W += Math.max(0, nrm * xcorr(S.co[i], S.co[j], g));
      w = W > 0 ? dd / W : 0; }
    const a = k * Fij * dFij * w; alpha[i][j] = alpha[j][i] = a;                                             // 回す角 α_ij = κ F Δ_F w（二チャネル写像の角を向きの回転角と読む）
    if (Math.abs(a) < 1e-300) continue;
    const M = m[i] + m[j], Ptot = v3add(pm[i], pm[j]), q = v3mul(v3sub(v3mul(pm[i], m[j]), v3mul(pm[j], m[i])), 1 / M);   // 相対運動量 q = (m_j p_i − m_i p_j)/M
    const ql = v3len(q); if (ql < 1e-12) continue;
    const qh = v3mul(q, 1 / ql), rr = v3sub(S.r[i], S.r[j]), b = v3sub(rr, v3mul(qh, dot(rr, qh)));              // 衝突径数 b = 相対位置の q に垂直な成分
    let e = null; if (v3len(b) > 1e-9) e = v3mul(crs(qh, v3mul(b, 1 / v3len(b))), 1); else { const cx = crs(S.n[i], S.n[j]); if (v3len(cx) > 1e-9) e = v3mul(cx, 1 / v3len(cx)); }   // 散乱面の法線 e = q̂ × b̂（なければ n_i × n_j）
    if (!e) continue;                                                                                        // 散乱面が決まらない（完全に同一直線）なら回さない
    E[i][j] = E[j][i] = e;
    const q2 = rotAxis(q, e, a);                                                                             // 重心系で q を e まわりに α 回す（P は不変）
    const pi2 = v3add(v3mul(Ptot, m[i] / M), q2), pj2 = v3sub(v3mul(Ptot, m[j] / M), q2);                     // p_i′ = (m_i/M) P + q′、p_j′ = (m_j/M) P − q′
    dp[i] = v3add(dp[i], v3sub(pi2, pm[i])); dp[j] = v3add(dp[j], v3sub(pj2, pm[j]));                         // 全対の変化を足す（順序なし。対ごとに Σp 保存なので全体も保存）
  }
  const n2 = [], v2 = [];
  for (let i = 0; i < N; i++) { const pn = v3add(pm[i], dp[i]), un = v3mul(pn, 1 / m[i]), sp = v3len(un);     // 新しい速度ベクトル u_i′ = p_i′/m_i
    if (sp < 1e-12) { n2.push(S.n[i]); v2.push(S.v[i]); continue; }
    const sg = S.v[i] < 0 ? -1 : 1; n2.push(v3mul(un, sg / sp)); v2.push(sg * sp); }                          // v の符号の規約は保つ：n_i′ = sign(v_i) u′/|u′|、v_i′ = sign(v_i)|u′|。反転なら n → −n
  const K = (vv, nn) => vv.reduce((t, v, i) => t + 0.5 * m[i] * v * v, 0), Psum = pp => pp.reduce((t, p) => v3add(t, p), [0, 0, 0]);
  return {n: n2, v: v2, F, dF, alpha, E, dp, Pbefore: Psum(pm), Pafter: Psum(n2.map((n, i) => v3mul(n, m[i] * v2[i]))), Kbefore: K(S.v), Kafter: K(v2)};
}
const xdesc = () => { const kv = Number(kap.value).toFixed(2); if (!lastX) return `κ = ${kv}、交換 ${xmOf() === 'rate' ? '率' : '全量'}（まだ交換なし）`;
  const P = lastX.Pafter; return `κ = ${kv}、交換 ${xmOf() === 'rate' ? '率' : '全量'}　核 F: ` + PAIRS.map(([i, j], k) => `${PAIRL[k]} ${lastX.F[i][j].toFixed(4)}`).join(', ') + `　Δ_F: ` + PAIRS.map(([i, j], k) => `${PAIRL[k]} ${(lastX.dF[i][j] * 180 / Math.PI).toFixed(1)}°`).join(', ') + `　α: ` + PAIRS.map(([i, j], k) => `${PAIRL[k]} ${(lastX.alpha[i][j] * 180 / Math.PI).toFixed(2)}°`).join(', ') + `　Σ m v n = (${P.map(x => x.toFixed(4)).join(', ')})　Σ ½ m v² = ${lastX.Kafter.toFixed(6)}`; };
function selfTest() {   // ?test=1：K1〜K7（算出根拠 v0.5 §9）。結果は <pre id="ktest">。包絡は原本の式（φ, M, 偶奇はスライダーの値）
  const out = {}, x = [1, 0, 0], y = [0, 1, 0], mx = [-1, 0, 0], ex = [1, 0, 0.05], ey = [0, 1, 0.05];
  const un = v => v3mul(v, 1 / v3len(v));
  const mk = (Th, n, v, A, phi0, r, mode) => ({Th, n: n.map(un), v, A, phi0, r, co: P.map(u => coeffs(Number(u.ph.value), Number(u.mm.value), u.odd.checked, u.even.checked)), mode: mode || 'full'});
  const pos = (Th, n) => Th.map((T, i) => v3mul(un(n[i]), T));
  // K3：同形同大（A = 1, 1）、X 軸の正面、わずかな z オフセット、重なり最大、Δ_F = 180。対 3–4 は軸 −x で見えない
  let n0 = [x, ex, mx, mx], Th0 = [110, 110, 220, 190]; let S = mk(Th0, n0, [1, -1, 2, -1], [1, 1, 1, 1], [0, 180, 30, 30], pos(Th0, n0)), X = exchangeStep(S);
  out.K1_P_conserved = v3len(v3sub(X.Pafter, X.Pbefore)); out.K1_K_conserved_rel = Math.abs(X.Kafter - X.Kbefore) / X.Kbefore;
  out.K3_F12 = X.F[0][1]; out.K3_alpha12_deg = X.alpha[0][1] * 180 / Math.PI; out.K3_n1_after = X.n[0]; out.K3_n2_after = X.n[1]; out.K3_v_after = X.v.slice(0, 2);
  const u0 = S.n.map((n, i) => v3mul(n, S.v[i])), u1 = X.n.map((n, i) => v3mul(n, X.v[i]));
  out.K3_reversed = dot(u1[0], u0[0]) < 0 && dot(u1[1], u0[1]) < 0; out.K3_swap_err = Math.max(v3len(v3sub(u1[0], u0[1])), v3len(v3sub(u1[1], u0[0])));   // 同質量の正面：速度ベクトルの入替え
  // K4：全対 Δ_F = 0 → 何も変わらない
  S = mk(Th0, n0, [1, -1, 2, -1], [1, 1, 1, 1], [40, 40, 40, 40], pos(Th0, n0)); X = exchangeStep(S);
  out.K4_identity = Math.max(...X.n.map((n, i) => v3len(v3sub(n, S.n[i]))), ...X.v.map((v, i) => Math.abs(v - S.v[i])));
  // K5：質量 m = A² が違う正面（A = 1 と 2、v = ±1、Δ_F = 180）：古典の弾性衝突 v₁′ = ((m₁−m₂)v₁ + 2m₂v₂)/(m₁+m₂)
  S = mk(Th0, n0, [1, -1, 2, -1], [1, 2, 1, 1], [0, 180, 30, 30], pos(Th0, n0)); X = exchangeStep(S);
  const m1 = 1, m2 = 4, v1 = 1, v2 = -1, v1c = ((m1 - m2) * v1 + 2 * m2 * v2) / (m1 + m2), v2c = ((m2 - m1) * v2 + 2 * m1 * v1) / (m1 + m2);   // 速度ベクトルの x 成分（ほぼ x 軸）
  const ux1 = X.v[0] * X.n[0][0], ux2 = X.v[1] * X.n[1][0];
  out.K5_classical_vx = [v1c, v2c]; out.K5_model_vx = [ux1, ux2]; out.K5_err = Math.max(Math.abs(ux1 - v1c), Math.abs(ux2 - v2c)); out.K5_P = [X.Pbefore, X.Pafter]; out.K5_K_rel = Math.abs(X.Kafter - X.Kbefore) / X.Kbefore;
  // K6：添字の置換
  const perm = [2, 0, 3, 1], Sp = mk(perm.map(i => S.Th[i]), perm.map(i => S.n[i]), perm.map(i => S.v[i]), perm.map(i => S.A[i]), perm.map(i => S.phi0[i]), perm.map(i => S.r[i])); Sp.co = perm.map(i => S.co[i]);
  const Xp = exchangeStep(Sp);
  out.K6_perm_max_diff = Math.max(...perm.map((i, k) => Math.max(v3len(v3sub(Xp.n[k], X.n[i])), Math.abs(Xp.v[k] - X.v[i]))));
  // K7：率の通過（2 波、X 軸、z オフセット、Δ_F = 150、相対速度 2 と 1）：終状態の向きが速度に依らないか。全量との対比
  const passage = (mode, sp, nt) => { let T = {Th: [50, 170], n: [x, un(ex)], v: [sp, -sp], A: [1, 1], phi0: [0, 150], co: [S.co[0], S.co[1]], mode}; T.r = pos(T.Th, T.n);
    const tr = []; let cum = 0;
    for (let t = 0; t < nt; t++) { T.Th = [T.Th[0] + T.v[0], T.Th[1] + T.v[1]]; T.r = T.r.map((r, i) => v3add(r, v3mul(T.n[i], T.v[i]))); const Y = exchangeStep(T); T = {Th: T.Th, n: Y.n, v: Y.v, A: T.A, phi0: T.phi0, r: T.r, co: T.co, mode};
      const a = Y.alpha[0][1] * 180 / Math.PI; if (Math.abs(a) > 1e-3) { cum += a; tr.push([t, Math.abs(T.Th[0] - T.Th[1]), +a.toFixed(2), Y.E[0][1] ? +Y.E[0][1][1].toFixed(3) : null, +cum.toFixed(1)]); } }   // [tick, d, α°, e_y, Σα°]
    return {n1: T.n[0], n2: T.n[1], v: T.v, turn1_deg: Math.acos(Math.max(-1, Math.min(1, dot(T.n[0], x)))) * 180 / Math.PI, sum_alpha_deg: cum, trace: tr}; };
  const r1 = passage('rate', 1, 120), r2 = passage('rate', 0.5, 240), f1 = passage('full', 1, 120), f2 = passage('full', 0.5, 240);
  out.K7_rate_v1 = r1; out.K7_rate_v05 = r2; out.K7_rate_speed_dep_deg = Math.abs(r1.turn1_deg - r2.turn1_deg);
  out.K7_full_v1 = f1; out.K7_full_v05 = f2; out.K7_full_speed_dep_deg = Math.abs(f1.turn1_deg - f2.turn1_deg);
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
    "  if (Number(kap.value) > 0) {   // 相互作用：自由運動のあとにゲージ回転。変わるのは n_i と v_i だけ（A・φ₀・位相 Θ は変えない）\n"
    "    const Th = P.map((u, i) => base[i][0] + cum[i][cum[i].length - 1][0]);\n"
    "    const co = P.map(u => coeffs(Number(u.ph.value), Number(u.mm.value), u.odd.checked, u.even.checked));   // 包絡は原本の式そのまま\n"
    "    const S = {Th, n: nDyn.map(x => x.slice()), v: vDyn.slice(), A: P.map(u => Number(u.aa.value)), phi0: P.map(u => u.phase0), r: rDyn.map(x => x.slice()), co, mode: xmOf()};\n"
    "    const X = exchangeStep(S); lastX = X;\n"
    "    P.forEach((u, i) => { nDyn[i] = X.n[i]; vDyn[i] = X.v[i]; });\n"
    "    ihist.push({t: ticks, F: X.F, dF: X.dF, alpha: X.alpha, n: X.n.map(v => v.slice()), v: X.v.slice(), P: X.Pafter, K: X.Kafter}); if (ihist.length > 4000) ihist.shift();\n"
    "  }\n  draw();\n}\nfunction setAuto(on) {")
rep("const resetAll = () => { ticks = 0; orb = [[], [], [], []]; orbTick = -1; hist = PAIRS.map(() => [[0, 0]]); cum = [0, 1, 2, 3].map(() => [[0, 0]]); win0 = null; for (const w of P) w.phase0 = 0; base = P.map(u => [Number(u.th.value), 0]); draw(); };",
    "const initDyn = () => { nDyn = P.map(u => nvec(u)); vDyn = P.map(u => Number(u.vv.value)); rDyn = P.map((u, i) => nDyn[i].map(x => x * Number(u.th.value))); ihist = []; lastX = null; };   // 動的状態の初期化：n・v はスライダー、r = θ n\n"
    "const resetAll = () => { ticks = 0; orb = [[], [], [], []]; orbTick = -1; hist = PAIRS.map(() => [[0, 0]]); cum = [0, 1, 2, 3].map(() => [[0, 0]]); win0 = null; for (const w of P) w.phase0 = 0; base = P.map(u => [Number(u.th.value), 0]); initDyn(); draw(); };\n"
    "kap.addEventListener('input', resetAll);   // κ を変えたら最初から（resetAll の定義のあと）\n"
    "for (const r of document.querySelectorAll('input[name=xm]')) r.addEventListener('change', resetAll);")
rep("  if (q.has('s')) ss.value = q.get('s');   // s = 球の倍率（撮影用）\n",
    "  if (q.has('s')) ss.value = q.get('s');   // s = 球の倍率（撮影用）\n"
    "  if (q.has('kap')) kap.value = q.get('kap');   // kap = 相互作用の強さ（0 で相互作用なし版と同一）\n"
    "  if (q.has('xm')) document.querySelector(`input[name=xm][value=\"${q.get('xm')}\"]`).checked = true;   // xm = rate | full\n")
rep("  base = P.map(u => [Number(u.th.value), u.phase0]);\n  draw();\n  setAuto(q.get('auto') !== '0');\n",
    "  base = P.map(u => [Number(u.th.value), u.phase0]); initDyn();\n  draw();\n  if (q.get('test') === '1') selfTest();\n  setAuto(q.get('auto') !== '0');\n")
assert "coeffsCos2" not in s and 'name="env"' not in s and "aDyn" not in s
io.open(p, "w", encoding="utf-8").write(s); print("ok")
