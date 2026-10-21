import sys, pathlib, subprocess, re, shutil
# 使い方：コンパイル用フォルダで python build_tex.py <論文フォルダ> → ja.tex, en.tex（その後 lualatex を 2 回）
src = pathlib.Path(sys.argv[1])
def prep(name, figdir, title, subtitle, author, date, out):
    md = (src / name).read_text(encoding="utf-8").splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    head = md[2:md.index("---")]  # 冒頭の著者・日付・DOI の行は行ごとに改行する
    head = [l + "\\" if l.startswith("**") and i < len(head) - 1 and head[i + 1].startswith("**") else l for i, l in enumerate(head)]
    body = "\n".join(head + md[md.index("---"):]).replace("──", "――").replace("ℓ", "$\\ell$").replace("≅", "$\\cong$")
    # 図の直後の説明文の段落（**図 N** …／**Figure N** …）を図の説明文にまとめる。番号は LaTeX が付ける
    body = re.sub(r"!\[[^\]]*\]\((?:図表/en/)?([^)/]+)\.svg\)\n\n\*\*(?:図 \d+|Figure \d+)\*\*\s*(.+)", r"![\2](figs/\1.png)", body)
    assert "](fig" not in body.replace("](figs/", "")
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (title, subtitle, author, date)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
    pathlib.Path("figs").mkdir(exist_ok=True)
    for p in (src / figdir).glob("*.png"):
        shutil.copy(p, "figs/" + p.name)
prep("quaternion_crystallization_inflation_ja_v1.1.md", "図表/ja",
     "四元数の関係読み出しと、結晶化としてのインフレーション",
     "― 無名性を保った自律的構造化と部分投射 ―",
     "木原 範昭（WF System Co., Ltd.）　ORCID 0009-0004-6753-4020", "2026年10月21日　Version DOI 10.5281/zenodo.23523238", "ja.md")
prep("quaternion_crystallization_inflation_en_v1.1.md", "図表/en",
     "Relational Readout of Quaternions and Crystallization as Inflation",
     "― Anonymity-Preserving Autonomous Structuring and Partial Projection ―",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "October 21, 2026  Version DOI 10.5281/zenodo.23523238", "en.md")
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
