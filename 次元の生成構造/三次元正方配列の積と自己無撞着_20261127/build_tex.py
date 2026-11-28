import sys, pathlib, subprocess, re, shutil
# Usage: run in a temp build directory (e.g. /tmp/build_tri_loop_tex):
#   /usr/bin/python3 build_tex.py <paper_dir>
# Generates ja.tex, en.tex and copies figures into figs_ja/, figs_en/.
# Then run: lualatex ja.tex (x2) and lualatex en.tex (x2).

src = pathlib.Path(sys.argv[1]).resolve()
BB = {"ℂ": "$\\mathbb{C}$", "ℍ": "$\\mathbb{H}$", "𝕆": "$\\mathbb{O}$", "𝕊": "$\\mathbb{S}$"}


def bb(s):
    for k, v in BB.items():
        s = s.replace(k, v)
    return s


def prep(name, title, subtitle, author, date_plain, date_fmt, out, figdir):
    text = (src / name).read_text(encoding="utf-8")
    m = re.search(r"Version DOI[：:]\*\*\s*(10\.5281/zenodo\.\d+)", text)
    date = date_fmt % m.group(1) if m else date_plain
    md = text.splitlines()
    assert md[0].startswith("# ") and md[1].startswith("## ")
    sep_idx = md.index("---")
    head = md[2:sep_idx]
    head = [
        l.rstrip() + "\\"
        if l.startswith("**") and i < len(head) - 1 and head[i + 1].startswith("**")
        else l
        for i, l in enumerate(head)
    ]
    body = "\n".join(head + md[sep_idx:]).replace("──", "――").replace("ℓ", "$\\ell$").replace("≅", "$\\cong$")
    pathlib.Path(figdir).mkdir(exist_ok=True)

    def fig(mm):
        cap = mm.group(1).strip()
        rel = mm.group(2).strip()
        png_src = src / rel
        assert png_src.exists(), f"Missing figure: {png_src}"
        png_dst = pathlib.Path(figdir) / png_src.name
        shutil.copy2(png_src, png_dst)
        return "![%s](%s){width=95%%}" % (cap, png_dst)

    body, n = re.subn(r"!\[([^\]]*)\]\(([^)]+\.png)\)", fig, body)
    assert n == 5, f"Expected 5 figures, found {n}"
    body = bb(body)
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (
        bb(title),
        bb(subtitle),
        author,
        date,
    )
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")


prep(
    "三次元正方配列の三角形ループ積による結晶選択と融解_速報.md",
    "三次元正方配列の三角形ループ積による四元数・八元数結晶の選択、カスケード融解、および視野角スペクトルの二極分裂（速報）",
    "― 背景時空・外部演算子・乱数を置かない自己更新系からの粒子凝縮と空間希薄化の数値実証 ―",
    "木原 範昭（株式会社WFシステム）　ORCID 0009-0004-6753-4020",
    "2026年11月28日",
    "2026年11月28日　Version DOI %s",
    "ja.md",
    "figs_ja",
)
prep(
    "triangle_loop_product_crystal_selection_melting_rapid_en.md",
    "Spontaneous Selection of Quaternionic and Octonionic Crystals, Cascade Melting, and Bimodal View-Angle Splitting under the Triangle-Loop Product on Cubic Arrays (Rapid Communication)",
    "― Numerical Demonstration of Particle Condensation and Spatial Rarefaction in a Self-Updating System without Background Spacetime, External Operators, or Randomness ―",
    "Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020",
    "November 28, 2026",
    "November 28, 2026  Version DOI %s",
    "en.md",
    "figs_en",
)

common = [
    "-s",
    "-f",
    "markdown",
    "--top-level-division=section",
    "-V",
    "papersize=a4",
    "-V",
    "geometry:margin=22mm",
    "-V",
    "mainfont=DejaVu Serif",
    "-V",
    "sansfont=DejaVu Sans",
    "-V",
    "monofont=DejaVu Sans Mono",
    "-V",
    "colorlinks=true",
    "-V",
    "header-includes=\\usepackage{xurl}",
]
subprocess.run(
    ["/opt/homebrew/bin/pandoc", "ja.md", "-o", "ja.tex", "-V", "documentclass=ltjsarticle"] + common,
    check=True,
)
subprocess.run(
    [
        "/opt/homebrew/bin/pandoc",
        "en.md",
        "-o",
        "en.tex",
        "-V",
        "documentclass=article",
        "-V",
        "header-includes=\\usepackage{luatexja}",
    ]
    + common,
    check=True,
)
for n in ("ja.tex", "en.tex"):
    p = pathlib.Path(n)
    t = p.read_text(encoding="utf-8")
    t = re.sub(
        r"\\texttt\{([^{}]*)\}",
        lambda m: "\\texttt{"
        + m.group(1).replace("\\_", "\\_\\allowbreak{}").replace("/", "/\\allowbreak{}")
        + "}",
        t,
    )
    t = re.sub(r"(?<![\w{/])(https?://[^\s{}]+)", r"\\url{\1}", t)
    p.write_text(
        t.replace("\\begin{document}", "\\ltjsetparameter{jacharrange={-2,-3}}\n\\begin{document}", 1),
        encoding="utf-8",
    )
print("tex ready")
