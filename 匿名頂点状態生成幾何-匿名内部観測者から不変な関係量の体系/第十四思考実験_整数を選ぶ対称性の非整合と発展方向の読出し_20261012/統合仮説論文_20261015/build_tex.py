import sys, pathlib, subprocess, re
src = pathlib.Path(sys.argv[1])
def prep(name, title, subtitle, author, date, out):
    md = (src / name).read_text(encoding="utf-8").splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    body = "\n".join(md[2:]).replace("](figures/", "](figs/").replace("──", "――")
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (title, subtitle, author, date)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
prep("unified_hypothesis_paper_ja_v1.1.md",
     "背景時空を置かない複素波の関係と干渉から、時間・空間・質量・電荷・並進速度・スピンを読み出す",
     "――閉じた円周の住民から四次元球殻の住民へ、そして α・四次元格子・E₈ に至る、思考実験の順序に沿った統合仮説――",
     "木原 範昭（WF System Co., Ltd.）　ORCID 0009-0004-6753-4020", "2026年10月15日　Version DOI 10.5281/zenodo.23386447", "ja.md")
prep("unified_hypothesis_paper_en_v1.1.md",
     "Reading Out Time, Space, Mass, Charge, Translational Velocity and Spin from the Relations and Interference of Complex Waves without a Background Spacetime",
     "― An Integrated Hypothesis in the Order of the Thought Experiment, from Residents of a Closed Circle to Residents of a Four-Dimensional Spherical Shell, and on to α, the Four-Dimensional Lattice and E₈ ―",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "October 15, 2026  Version DOI 10.5281/zenodo.23386447", "en.md")
common = ["-s", "-f", "markdown", "--top-level-division=section", "-V", "papersize=a4", "-V", "geometry:margin=25mm",
          "-V", "mainfont=DejaVu Serif", "-V", "sansfont=DejaVu Sans", "-V", "monofont=DejaVu Sans Mono", "-V", "colorlinks=true",
          "-V", "header-includes=\\usepackage{xurl}"]
subprocess.run(["pandoc", "ja.md", "-o", "ja.tex", "-V", "documentclass=ltjsarticle"] + common, check=True)
subprocess.run(["pandoc", "en.md", "-o", "en.tex", "-V", "documentclass=article", "-V", "header-includes=\\usepackage{luatexja}"] + common, check=True)
for n in ("ja.tex", "en.tex"):
    p = pathlib.Path(n); t = p.read_text(encoding="utf-8")
    t = re.sub(r"\\texttt\{([^{}]*)\}", lambda m: "\\texttt{" + m.group(1).replace("\\_", "\\_\\allowbreak{}").replace("/", "/\\allowbreak{}") + "}", t)
    t = re.sub(r"(?<![\w{/])(https?://[^\s{}]+)", r"\\url{\1}", t)
    p.write_text(t.replace("\\begin{document}", "\\ltjsetparameter{jacharrange={-2,-3}}\n\\begin{document}", 1), encoding="utf-8")
print("tex ready")
