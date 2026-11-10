# -*- coding: utf-8 -*-
"""図 1・図 2：4 パルス実験（index_v12q4i.html）の既存の静止画から左右のパネルを切り抜く。
元の静止画は実験フォルダの shoot.sh が撮ったもの（既定の初期条件、tick 130、1800×1400）。
画素の加工はしない。切り抜きだけ。
"""
import os
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "有界な定常波の再検討_閉じた弦のパルス_フーリエ級数_20261028",
                   "四元数化の検討_4パルス_相互作用版_20261109", "figs", "v12q4i_default_tick130.png")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
im = Image.open(SRC)
assert im.size == (1800, 1400), im.size
im.crop((60, 900, 800, 1400)).save(os.path.join(OUT, "fig1_experiment_left_panel_power_on_circle.png"))
im.crop((960, 900, 1700, 1400)).save(os.path.join(OUT, "fig2_experiment_right_panel_internal_phase_space.png"))
print("saved fig1, fig2")
