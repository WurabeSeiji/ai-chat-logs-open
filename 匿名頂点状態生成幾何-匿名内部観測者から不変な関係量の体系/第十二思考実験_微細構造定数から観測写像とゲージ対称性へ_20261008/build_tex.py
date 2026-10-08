import sys, pathlib, subprocess
src = pathlib.Path(sys.argv[1])
def prep(name, title, subtitle, author, date, out):
    md = (src / name).read_text(encoding="utf-8").splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    body = "\n".join(md[2:]).replace("](figures/", "](figs/").replace("──", "――")
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (title, subtitle, author, date)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
prep("thought_experiment_12_observation_mapping_gauge_symmetry_ja_v1.0.md",
     "第十二思考実験：微細構造定数から観測写像とゲージ対称性へ",
     "― 微細構造定数から始めて、光円錐の方程式と観測可能な軸の選択から、時空・質量・電荷・色荷・ゲージ対称性を読み直す ―",
     "木原 範昭（WF System Co., Ltd.）　ORCID 0009-0004-6753-4020", "2026年10月8日　Version DOI 10.5281/zenodo.23226259", "ja.md")
prep("thought_experiment_12_observation_mapping_gauge_symmetry_en_v1.0.md",
     "The Twelfth Thought Experiment: From the Fine-Structure Constant to Observation Maps and Gauge Symmetry",
     "― Starting from the Fine-Structure Constant, Rereading Spacetime, Mass, Charge, Color Charge and Gauge Symmetry from the Light-Cone Equation and the Choice of Observable Axes ―",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "October 8, 2026  Version DOI 10.5281/zenodo.23226259", "en.md")
common = ["-s", "-f", "markdown", "--top-level-division=section", "-V", "papersize=a4", "-V", "geometry:margin=25mm",
          "-V", "mainfont=DejaVu Serif", "-V", "sansfont=DejaVu Sans", "-V", "monofont=DejaVu Sans Mono", "-V", "colorlinks=true"]
subprocess.run(["pandoc", "ja.md", "-o", "ja.tex", "-V", "documentclass=ltjsarticle"] + common, check=True)
subprocess.run(["pandoc", "en.md", "-o", "en.tex", "-V", "documentclass=article", "-V", "header-includes=\\usepackage{luatexja}"] + common, check=True)
for n in ("ja.tex", "en.tex"):
    p = pathlib.Path(n); t = p.read_text(encoding="utf-8")
    p.write_text(t.replace("\\begin{document}", "\\ltjsetparameter{jacharrange={-2,-3}}\n\\begin{document}", 1), encoding="utf-8")
print("tex ready")
