#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 3 速報の日英 TeX 生成スクリプト。
"""
import sys, pathlib, subprocess, re, shutil

src = pathlib.Path(sys.argv[1]).resolve()
BB = {"ℂ": "$\\mathbb{C}$", "ℍ": "$\\mathbb{H}$", "𝕆": "$\\mathbb{O}$", "𝕊": "$\\mathbb{S}$", "ℤ": "$\\mathbb{Z}$"}


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
    "三次元正方配列の三角形ループ積におけるスピン選別と倍音局在_第3速報.md",
    "三次元正方配列の三角形ループ積における SO(3) スピン選別、G2/SU(3) ゲージ対称性、空間倍音による局在化、および 2 粒子の完全共存（第 3 速報）",
    "― 空間回転表現の無質量・有質量分離 (1 vs 1/3)、離散 $1/\\sqrt{s+1}$ 則、3:3:1 倍音二項比、および勝者総取りを超えた複数粒子共存の厳密実証 ―",
    "木原 範昭（株式会社WFシステム）　ORCID 0009-0004-6753-4020",
    "2026年11月28日",
    "2026年11月28日　Version DOI %s",
    "ja.md",
    "figs_ja",
)
prep(
    "triangle_loop_spin_selection_harmonic_localization_rapid3_en.md",
    "SO(3) Spin Selection, G2/SU(3) Gauge Symmetries, Spatial Harmonic Localization, and Two-Particle Coexistence in the 3D Triangle-Loop Product (Third Rapid Communication)",
    "― Exact Massless/Gapped Split (1 vs 1/3), Discrete $1/\\sqrt{s+1}$ Law, 3:3:1 Spatial Harmonic Binomial Ratio, and Coexistence Beyond Winner-Take-All ―",
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
    t = t.replace("⊗", "$\\otimes$").replace("√", "$\\sqrt{}$")
    p.write_text(
        t.replace("\\begin{document}", "\\ltjsetparameter{jacharrange={-2,-3}}\n\\begin{document}", 1),
        encoding="utf-8",
    )
print("rapid3 tex ready")
