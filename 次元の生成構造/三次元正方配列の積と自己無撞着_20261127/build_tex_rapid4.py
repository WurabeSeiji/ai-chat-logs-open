#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 4 速報の日英 TeX 生成スクリプト。
"""
import sys, pathlib, subprocess, re, shutil

src = pathlib.Path(sys.argv[1]).resolve()
BB = {"ℂ": "$\\mathbb{C}$", "ℍ": "$\\mathbb{H}$", "𝕆": "$\\mathbb{O}$", "𝕊": "$\\mathbb{S}$", "ℤ": "$\\mathbb{Z}$", "ℝ": "$\\mathbb{R}$"}


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
    body = (
        "\n".join(head + md[sep_idx:])
        .replace("──", "――")
        .replace("ℓ", "$\\ell$")
        .replace("≅", "$\\cong$")
        .replace("→", "$\\to$")
    )
    pathlib.Path(figdir).mkdir(exist_ok=True)

    def fig(mm):
        cap = mm.group(1).strip()
        rel = mm.group(2).strip()
        png_src = src / rel
        assert png_src.exists(), f"Missing figure: {png_src}"
        png_dst = pathlib.Path(figdir) / png_src.name
        shutil.copy2(png_src, png_dst)
        return "![%s](%s){width=96%%}" % (cap, png_dst)

    body, n = re.subn(r"!\[([^\]]*)\]\(([^)]+\.png)\)", fig, body)
    assert n == 2, f"Expected 2 figures, found {n}"
    body = bb(body)
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (
        bb(title),
        bb(subtitle),
        author,
        date,
    )
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")


prep(
    "三次元正方配列の三角形ループ積からの量子構造の創発_第4速報.md",
    "三次元正方配列の三角形ループ積からの量子構造の創発：厳密離散化、単一複素構造の自発選択、Born 則、および離散不確定性飽和（第 4 速報）",
    "― 量子力学の公理を仮定しない実数 3 階テンソル反復からの複素ヒルベルト構造・確率二乗則・最小分解能の解析導出 ―",
    "木原 範昭（株式会社WFシステム）　ORCID 0009-0004-6753-4020",
    "2026年11月28日",
    "2026年11月28日　Version DOI %s",
    "ja.md",
    "figs_ja",
)
prep(
    "triangle_loop_emergent_quantum_structure_rapid4_en.md",
    "Emergent Quantum Structure from the 3D Triangle-Loop Product: Exact Discretization, Spontaneous Single Complex Structure, Born Rule, and Discrete Uncertainty Saturation (Fourth Rapid Communication)",
    "― Analytical Derivation of Complex Hilbert Structure, Probability-Squared Law, and Minimum Resolution from a Real 3-Tensor Without Quantum Postulates ―",
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
    "geometry:margin=20mm",
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
for fn in ["ja.tex", "en.tex"]:
    s = pathlib.Path(fn).read_text(encoding="utf-8")
    s = re.sub(r"\\begin\{figure\}(\[[^\]]*\])?", r"\\begin{figure}[H]", s)
    if "\\usepackage{float}" not in s:
        s = s.replace("\\usepackage{xurl}", "\\usepackage{xurl}\n\\usepackage{float}")
    pathlib.Path(fn).write_text(s, encoding="utf-8")
print("Generated ja.tex and en.tex successfully.")
