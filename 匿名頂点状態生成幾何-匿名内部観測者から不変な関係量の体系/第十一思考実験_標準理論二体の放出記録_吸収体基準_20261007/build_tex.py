import sys, pathlib, subprocess
src = pathlib.Path(sys.argv[1])
def prep(name, title, subtitle, author, date, out):
    md = (src / name).read_text(encoding="utf-8").splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    body = "\n".join(md[2:]).replace("](figures/", "](figs/").replace("──", "――")
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (title, subtitle, author, date)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
prep("thought_experiment_11_emission_record_standard_two_body_ja_v1.0.md",
     "第十一思考実験：背景時空を置かない二体系の放出記録",
     "― 標準理論（Dirac–Coulomb 束縛状態・1PN 二体重力・Einstein の A 係数）を無限遠の吸収体の静止系だけを基準に数値計算し、水素原子・電子–電子・陽子–陽子・水素分子の 4 系で放出と吸収の関係量を読み出す ―",
     "木原 範昭（WF System Co., Ltd.）　ORCID 0009-0004-6753-4020", "2026年10月7日　Version DOI 10.5281/zenodo.23199300", "ja.md")
prep("thought_experiment_11_emission_record_standard_two_body_en_v1.0.md",
     "The Eleventh Thought Experiment: Emission Records of Two-Body Systems without a Background Spacetime",
     "- Computing Standard Theory (Dirac–Coulomb Bound States, 1PN Two-Body Gravity, Einstein A Coefficients) Referenced Only to the Rest Frame of an Absorber at Infinity, and Reading Out Emission–Absorption Relations for the Hydrogen Atom, Electron–Electron, Proton–Proton and the Hydrogen Molecule -",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "October 7, 2026  Version DOI 10.5281/zenodo.23199300", "en.md")
common = ["-s", "-f", "markdown", "--top-level-division=section", "-V", "papersize=a4", "-V", "geometry:margin=25mm",
          "-V", "mainfont=DejaVu Serif", "-V", "sansfont=DejaVu Sans", "-V", "monofont=DejaVu Sans Mono", "-V", "colorlinks=true"]
subprocess.run(["pandoc", "ja.md", "-o", "ja.tex", "-V", "documentclass=ltjsarticle"] + common, check=True)
subprocess.run(["pandoc", "en.md", "-o", "en.tex", "-V", "documentclass=article", "-V", "header-includes=\\usepackage{luatexja}"] + common, check=True)
for n in ("ja.tex", "en.tex"):
    p = pathlib.Path(n); t = p.read_text(encoding="utf-8")
    p.write_text(t.replace("\\begin{document}", "\\ltjsetparameter{jacharrange={-2,-3}}\n\\begin{document}", 1), encoding="utf-8")
print("tex ready")
