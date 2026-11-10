import sys, pathlib, subprocess, re, shutil
# 使い方：コンパイル用フォルダ（/tmp 下）で python3 build_tex.py <論文フォルダ> → ja.tex, en.tex（その後 lualatex を 2 回ずつ）
# 前稿（観測は放射と吸収の同時成立である、2026-11-01）の build_tex.py の忠実コピー。変えたのは、入力ファイル名・題名・日付の引数と、
# 英題の ' の重ね書き、図の扱いの 2 か所だけ（図が PNG で二つのフォルダにあるので、正規表現のパスと拡張子、figs/ へ写すフォルダの一覧を替えた）。
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
    body = re.sub(r"!\[[^\]]*\]\((?:[^)]*/)?([^)/]+)\.png\)\n\n\*\*(?:図 \d+|Figure \d+|Fig\. \d+)\*\*\s*(.+)", r"![\2](figs/\1.png)", body)
    assert "](fig" not in body.replace("](figs/", "")
    body = re.sub(r"\n(?=\[(?:S|X)\d+\] )", "\n\n", body)  # 参考文献の各項目を別の段落にする（md では 1 行ずつ）
    body = bb(body)
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (bb(title), bb(subtitle), author, date_fmt % doi)
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")
    pathlib.Path("figs").mkdir(exist_ok=True)
    for d in (src / "図表" / "figures", src / "概念図_次数を逆に読む" / "figs"):
        for p in d.glob("*.png"):
            shutil.copy(p, "figs/" + p.name)
prep("paper_resolution_gauge_reverse_reading_ja_v1.0.md",
     "観測者の分解能はゲージである",
     "― 次数を逆に読む思考実験と、既知の対称性から組み上がる理論の予言 ―",
     "木原 範昭（WF System Co., Ltd.）　ORCID 0009-0004-6753-4020", "2026年11月10日　Version DOI %s", "ja.md")
prep("paper_resolution_gauge_reverse_reading_en_v1.0.md",
     "The Observer''s Resolution Is a Gauge",  # YAML の単引用符の中なので ' を二つ重ねる
     "― A Thought Experiment That Reads the Harmonic Order in Reverse, and the Prediction of a Theory Assembled from Known Symmetries ―",
     "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020", "November 10, 2026  Version DOI %s", "en.md")
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
