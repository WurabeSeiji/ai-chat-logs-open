# -*- coding: utf-8 -*-
"""note 記事用の図。論文の図 1〜3 の PNG を切り分けて並べるだけで、画素の加工はしない。
使い方：論文フォルダで python3 note記事_観測者の分解能はゲージである_20261110/make_note_figs.py
"""
import os
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# 図 A：数値実験の左右のパネルを横に並べる
l = Image.open(os.path.join(ROOT, "図表/figures/fig1_experiment_left_panel_power_on_circle.png")).convert("RGB")
r = Image.open(os.path.join(ROOT, "図表/figures/fig2_experiment_right_panel_internal_phase_space.png")).convert("RGB")
s = Image.new("RGB", (l.width + r.width, max(l.height, r.height)), "white")
s.paste(l, (0, 0)); s.paste(r, (l.width, 0))
s.save(os.path.join(HERE, "note_fig1_experiment_two_panels.png"))
# 図 B〜F：概念図（論文の図 3）をパネルごとに切る。切る行は、パネルの間の白い帯の中
c = Image.open(os.path.join(ROOT, "概念図_次数を逆に読む/figs/fig_reverse_reading_v0.1.png")).convert("RGB")
assert c.size == (1356, 1798), c.size
for name, y0, y1 in (("note_fig3_ladder_two_readings", 0, 390), ("note_fig2_first_reading_harmonics", 390, 670),
                     ("note_fig4_rigid_crystal_lowpass", 670, 965), ("note_fig5_coarse_wave_centre_phase", 965, 1260),
                     ("note_fig6_power_tooth_jump", 1260, 1798)):
    c.crop((0, y0, c.width, y1)).save(os.path.join(HERE, name + ".png"))
print("saved 6")
