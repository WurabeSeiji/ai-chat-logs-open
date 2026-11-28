#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
統合完成論文（Full Paper）の日英 TeX 生成スクリプト（全 10 図・目次付き）。
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
        return "![%s](%s){width=95%%}" % (cap, png_dst)

    body, n = re.subn(r"!\[([^\n]+)\]\(([^)]+\.png)\)", fig, body)
    assert n == 10, f"Expected 10 figures, found {n}"
    body = bb(body)
    yaml = "---\ntitle: '%s'\nsubtitle: '%s'\nauthor: '%s'\ndate: '%s'\n---\n\n" % (
        bb(title),
        bb(subtitle),
        author,
        date,
    )
    pathlib.Path(out).write_text(yaml + body, encoding="utf-8")


prep(
    "三次元正方配列の三角形ループ積からの量子・ゲージ・重力・フェルミオン構造の自発創発_統合完成論文.md",
    "三次元正方配列の三角形ループ積からの量子・ゲージ・重力・フェルミオン構造の自発創発",
    "― 背景不要な実数 3 階テンソル反復の統合解析（第 1〜第 7 速報の体系的総括） ―",
    "木原 範昭（株式会社WFシステム）　ORCID 0009-0004-6753-4020",
    "2026年11月28日",
    "2026年11月28日　Version DOI %s",
    "ja.md",
    "figs_ja",
)
prep(
    "triangle_loop_emergent_quantum_gauge_gravity_fermion_unified_full_paper_en.md",
    "Spontaneous Emergence of Quantum, Gauge, Geometric, and Fermionic Structures from the 3D Triangle-Loop Product",
    "― A Unified Analysis of a Background-Free Real 3-Tensor Iteration (Synthesis of Rapid Communications 1–7) ―",
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
    "--toc",
    "--toc-depth=2",
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
