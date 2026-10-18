import sys, pathlib, subprocess, re
src = pathlib.Path(sys.argv[1])
def prep(name, title, subtitle, author, date, out):
    md = (src / name).read_text(encoding="utf-8").splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    body = "\n".join(md[2:]).replace("](図表/figures/", "](figs/").replace("──", "――").replace("ℓ", "$\\ell$").replace("≅", "$\\cong$")
    body = re.sub(r"!\[(?:図 \d+　|Figure \d+\. )", "![", body)  # 図の番号は LaTeX が付けるので、説明文の先頭の番号を外す
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (title, subtitle, author, date)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
prep("quaternion_spacetime_map_ja_v1.md",
     "複素数の集合から時空 (1＋3) への写像",
     "― 四元数による読み出しと部分投射 ―",
     "木原 範昭（WF System Co., Ltd.）　ORCID 0009-0004-6753-4020", "2026年10月18日　Version DOI 10.5281/zenodo.23452319", "ja.md")
prep("quaternion_spacetime_map_en_v1.md",
     "A Map from a Set of Complex Numbers to Spacetime (1+3)",
     "― Readout by Quaternions and Partial Projection ―",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "October 18, 2026  Version DOI 10.5281/zenodo.23452319", "en.md")
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
