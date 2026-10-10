import sys, pathlib, subprocess, re
src = pathlib.Path(sys.argv[1])
def prep(name, title, subtitle, author, date, out):
    md = (src / name).read_text(encoding="utf-8").splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    body = "\n".join(md[2:]).replace("](figures/", "](figs/").replace("──", "――")
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (title, subtitle, author, date)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
prep("thought_experiment_13_coherence_harmonic_exchange_mass_lifetime_ja_v2.0.md",
     "第十三思考実験：質量を獲得した粒子は、なぜ散逸して消滅しないのか？",
     "――コヒーレンスの関係性、局在波の倍音交換、内部時計、粒子寿命を一つの相互作用から考える――",
     "木原 範昭（WF System Co., Ltd.）　ORCID 0009-0004-6753-4020", "2026年10月10日　Version DOI 10.5281/zenodo.23280070", "ja.md")
prep("thought_experiment_13_coherence_harmonic_exchange_mass_lifetime_en_v2.0.md",
     "The Thirteenth Thought Experiment: Why Does a Particle That Has Acquired Mass Not Dissipate and Vanish?",
     "― Coherence as a Relation, Harmonic Exchange of Localized Waves, Internal Clocks and Particle Lifetimes Considered from a Single Interaction ―",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "October 10, 2026  Version DOI 10.5281/zenodo.23280070", "en.md")
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
