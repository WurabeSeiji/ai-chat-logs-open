"""note 記事用の図（PNG）を作る。日本語版は figures/ の SVG から、英語版は同じ SVG の文字列を
make_figures_en.py の規則で英語に置き換えてから、rsvg-convert で PNG にする。図形と数値は変えない。
使い方：python3 make_note_figs.py
"""
import pathlib
import subprocess
import tempfile

import make_figures_en as en

HERE = pathlib.Path(__file__).parent
FIGS = [
    ("fig_cup_coin_schematic_covered", "schematic_covered"),
    ("fig_cup_coin_schematic_beside", "schematic_beside"),
    ("fig_cup_coin_s1_phase90", "peak_cup_over_coin"),
    ("cross_term/fig_cup_coin_s1_cross_phase090", "peak_cross_overlap"),
    ("cross_term/fig_cup_coin_s1_cross_phase060", "peak_cross_apart"),
    ("well_cross/fig_cup_well_cross_phase090", "well_hidden"),
    ("well_cross/fig_cup_well_cross_phase075", "well_partly_visible"),
    ("well_cross/fig_cup_well_cross_phase060", "well_visible"),
]
out = HERE / "note_figs"
out.mkdir(exist_ok=True)
for i, (rel, name) in enumerate(FIGS, 1):
    svg = HERE / "figures" / (rel + ".svg")
    subprocess.run(["rsvg-convert", "-w", "1400", "-b", "white", str(svg), "-o", str(out / ("note_fig%d_%s.png" % (i, name)))], check=True)
    text = en.fit(en.translate(svg.read_text(encoding="utf-8")))
    assert not en.JP.search(text), rel
    with tempfile.NamedTemporaryFile("w", suffix=".svg", encoding="utf-8", delete=False) as f:
        f.write(text)
    subprocess.run(["rsvg-convert", "-w", "1400", "-b", "white", f.name, "-o", str(out / ("note_fig%d_%s_en.png" % (i, name)))], check=True)
    pathlib.Path(f.name).unlink()
    print("ok", i, name)
