import sys, pathlib, subprocess, re
# 使い方：コンパイル用フォルダ（/tmp 下など）で python3 build_tex.py <論文フォルダ> → ja.tex, en.tex（その後 lualatex を 2 回ずつ）
# 「QCD を八元数で表現したら」(2026-11-19) の build_tex.py をもとにした。変えたのは、入力ファイル名・題名・日付だけである。
src = pathlib.Path(sys.argv[1])
BB = {"ℂ": "$\\mathbb{C}$", "ℍ": "$\\mathbb{H}$", "𝕆": "$\\mathbb{O}$"}
def bb(s):
    for k, v in BB.items():
        s = s.replace(k, v)
    return s
def prep(name, title, subtitle, author, date_plain, date_fmt, out):
    text = (src / name).read_text(encoding="utf-8")
    m = re.search(r"Version DOI[：:]\*\*\s*(10\.5281/zenodo\.\d+)", text)
    date = date_fmt % m.group(1) if m else date_plain
    md = text.splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    head = md[2:md.index("---")]  # 冒頭の著者・日付・DOI の行は行ごとに改行する
    head = [l + "\\" if l.startswith("**") and i < len(head) - 1 and head[i + 1].startswith("**") else l for i, l in enumerate(head)]
    body = "\n".join(head + md[md.index("---"):]).replace("──", "――").replace("ℓ", "$\\ell$").replace("≅", "$\\cong$")
    body = re.sub(r"\n(?=\[(?:S|X)\d+\] )", "\n\n", body)  # 参考文献の各項目を別の段落にする（md では 1 行ずつ）
    body = bb(body)
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (bb(title), bb(subtitle), author, date)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
prep("八元数の積からヤンミルズとディラックを読み出す.md",
     "八元数の積から、ヤン＝ミルズとディラックを読み出す",
     "― 方向を逆にした思考実験：作用・距離・結合係数が一つの数であるという仮説 ―",
     "木原 範昭（株式会社WFシステム）　ORCID 0009-0004-6753-4020", "2026年11月20日", "2026年11月20日　Version DOI %s", "ja.md")
prep("reading_out_yang_mills_dirac_from_octonion_product_en.md",
     "Reading Out Yang–Mills and Dirac from the Product of Octonions",
     "― A Thought Experiment in the Reverse Direction: The Hypothesis That Action, Distance, and Coupling Are One Number ―",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "November 20, 2026", "November 20, 2026  Version DOI %s", "en.md")
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
