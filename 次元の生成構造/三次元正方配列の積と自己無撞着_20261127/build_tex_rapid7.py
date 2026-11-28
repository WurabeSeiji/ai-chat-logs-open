#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 7 速報の日英 TeX 生成スクリプト。
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
    "三次元正方配列の三角形ループ積における八元数Trialityスピン1_2、CARフェルミオン電荷、および時間結晶カイラル分裂_第7速報.md",
    "三次元正方配列の三角形ループ積における八元数 Triality スピン $1/2$、CAR フェルミオン電荷スペクトル（$0, 1/3, 2/3, 1$）、および時間結晶カイラル分裂（第 7 速報）",
    "― 28 次元 $Spin(8)$ Triality と $2\\pi$ 回転符号反転、パウリ排他則・$SU(3)$ カシミールの 1:3:3:1 構造、および内側 $M_8$ によるカイラル非対称性 ―",
    "木原 範昭（株式会社WFシステム）　ORCID 0009-0004-6753-4020",
    "2026年11月28日",
    "2026年11月28日　Version DOI %s",
    "ja.md",
    "figs_ja",
)
prep(
    "triangle_loop_triality_spin_half_fermion_charge_chiral_rapid7_en.md",
    "Octonionic Spin(8) Triality Spin-1/2, CAR Fermion Charge Spectrum (0, 1/3, 2/3, 1), and Time-Crystal Chiral Splitting in the 3D Triangle-Loop Product (Seventh Rapid Communication)",
    "― Exact 28D Triality Sign Flip Under $2\\pi$ Rotation, Pauli Exclusion and $SU(3)$ Casimirs in a 1:3:3:1 Hierarchy, and Inner-$M_8$ Chiral Asymmetry ―",
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
