# -*- coding: utf-8 -*-
"""対照 1：7/13 の原本プログラム（原本/ に無変更コピー）を import し、その関数だけで片側高次倍音 (N_A, N_B) = (1, 63)、R = 0.55 の再帰散乱を
128 回走らせ、論文 §6 の表（衝突 0・16・32・64・128 の A/B の L と N_eff）を再現する。原本の recursive_transfer_rows と同じ手順（r a + t b、t a + r b、チャネルごとに normalize）。"""
import importlib.util, json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "原本", "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py")
spec = importlib.util.spec_from_file_location("orig", SRC); orig = importlib.util.module_from_spec(spec); sys.modules["orig"] = orig; spec.loader.exec_module(orig)   # dataclass（from __future__ annotations）のため sys.modules に登録してから実行
params = orig.Params()
delta_f = 2.0 * math.asin(math.sqrt(0.55))
t, r, T, R = orig.scattering_coefficients(delta_f)
a = orig.make_state(params, 1, params.q_A, params.m_A, True, params.A_A)
b = orig.make_state(params, 63, params.q_B, params.m_B, True, params.A_B)
# 論文 §6 の表（R = 0.55, T = 0.45）
paper = {0: (0.000183105, 1, 0.00520897, 32), 16: (0.00167148, 16.9941, 0.00151127, 16.0059), 32: (0.00519937, 31.9685, 0.000183722, 1.0315),
         64: (0.000185626, 1.12587, 0.00517067, 31.8741), 128: (0.000194012, 1.50145, 0.00505728, 31.4985)}
rows = {}; bad = 0
for k in range(params.recursive_collision_count + 1):
    if k in paper:
        La = orig.localization(a); Na = orig.effective_n(orig.harmonic_distribution(params, a))[0]
        Lb = orig.localization(b); Nb = orig.effective_n(orig.harmonic_distribution(params, b))[0]
        pa = paper[k]; ok = abs(La - pa[0]) < 1e-6 and abs(Na - pa[1]) < 1e-3 and abs(Lb - pa[2]) < 1e-6 and abs(Nb - pa[3]) < 1e-3
        bad += (not ok); rows[k] = {"A_L": La, "A_N_eff": Na, "B_L": Lb, "B_N_eff": Nb, "paper": pa, "match": ok}
        print(f"collision {k:3d}: A L={La:.6g} N_eff={Na:.5g} | B L={Lb:.6g} N_eff={Nb:.5g} | 論文 {pa} | {'一致' if ok else '不一致'}")
    if k >= params.recursive_collision_count: break
    a_next = orig.normalize(r * a + t * b); b_next = orig.normalize(t * a + r * b); a, b = a_next, b_next
print("R =", R, "T =", T, "| 対照 1", "一致" if bad == 0 else f"不一致 {bad}")
json.dump({"R": R, "T": T, "rows": rows, "ok": bad == 0}, open(os.path.join(HERE, "results", "control_original_v1.json"), "w"), ensure_ascii=False, indent=2)
sys.exit(0 if bad == 0 else 2)
