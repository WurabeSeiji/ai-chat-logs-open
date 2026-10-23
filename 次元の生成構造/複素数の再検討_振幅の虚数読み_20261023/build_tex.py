import sys, pathlib, subprocess, re, shutil
# 使い方：コンパイル用フォルダ（/tmp 下）で python3 build_tex.py <論文フォルダ> → ja.tex, en.tex（その後 lualatex を 2 回ずつ）
# 前稿（四元数の関係読み出しと、結晶化としてのインフレーション）の build_tex.py を元に、
# (a) 図の参照先を 図表/figures/*.svg（英語版は *_en.svg）に、(b) Version DOI を md のヘッダーから読む、
# (c) 黒板太字 ℂ ℍ 𝕆 を \mathbb に置く（本文フォントに字形がないため）、の 3 点を変えた。
src = pathlib.Path(sys.argv[1])
BB = {"ℂ": "$\\mathbb{C}$", "ℍ": "$\\mathbb{H}$", "𝕆": "$\\mathbb{O}$"}
def bb(s):
    for k, v in BB.items():
        s = s.replace(k, v)
    return s
def prep(name, title, subtitle, author, date_fmt, out):
    text = (src / name).read_text(encoding="utf-8")
    doi = re.search(r"Version DOI[：:]\*\*\s*(10\.5281/zenodo\.\d+)", text).group(1)
    md = text.splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    head = md[2:md.index("---")]  # 冒頭の著者・日付・DOI の行は行ごとに改行する
    head = [l + "\\" if l.startswith("**") and i < len(head) - 1 and head[i + 1].startswith("**") else l for i, l in enumerate(head)]
    body = "\n".join(head + md[md.index("---"):]).replace("──", "――").replace("ℓ", "$\\ell$").replace("≅", "$\\cong$")
    # 図の直後の説明文の段落（**図 N** …／**Figure N** …）を図の説明文にまとめる。番号は LaTeX が付ける
    body = re.sub(r"!\[[^\]]*\]\((?:図表/figures/)?([^)/]+)\.svg\)\n\n\*\*(?:図 \d+|Figure \d+)\*\*\s*(.+)", r"![\2](figs/\1.png)", body)
    assert "](fig" not in body.replace("](figs/", "")
    body = re.sub(r"\n(?=\[(?:S|X)\d+\] )", "\n\n", body)  # 参考文献の各項目を別の段落にする（md では 1 行ずつ）
    body = bb(body)
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (bb(title), bb(subtitle), author, date_fmt % doi)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
    pathlib.Path("figs").mkdir(exist_ok=True)
    for p in (src / "図表" / "figures").glob("*.png"):
        shutil.copy(p, "figs/" + p.name)
prep("complex_revisited_amplitude_imaginary_ja_v1.md",
     "複素数を正しく理解していたか",
     "― 振幅はどこに住むか、虚数読み、ℂ → ℍ → 𝕆 の帰属 ―",
     "木原 範昭（WF System Co., Ltd.）　ORCID 0009-0004-6753-4020", "2026年10月23日　Version DOI %s", "ja.md")
prep("complex_revisited_amplitude_imaginary_en_v1.md",
     "Did We Understand Complex Numbers Correctly?",
     "― Where the Amplitude Lives, the Imaginary Reading, and Attribution across ℂ → ℍ → 𝕆 ―",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "October 23, 2026  Version DOI %s", "en.md")
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
