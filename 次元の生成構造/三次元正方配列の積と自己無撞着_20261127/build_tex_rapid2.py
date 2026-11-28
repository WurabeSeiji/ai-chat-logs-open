#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 2 速報の日英 TeX 生成スクリプト（build_tex.py と同一フォーマット）。
Usage:
  mkdir -p /tmp/build_tri_loop_rapid2 && cd /tmp/build_tri_loop_rapid2
  /usr/bin/python3 <paper_dir>/build_tex_rapid2.py <paper_dir>
  lualatex ja.tex (x2) && lualatex en.tex (x2)
"""
import sys, pathlib, subprocess, re, shutil

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
    "三次元正方配列の三角形ループ積における4段階相転移とプローブ選別_第2速報.md",
    "三次元正方配列の三角形ループ積における過冷却時間結晶背景、プローブ選別、および 4 段階・倍々エンタングルメント相転移（第 2 速報）",
    "― N = 64 (O8 ⊗ O8) 結合大結晶の τ = 10,000 厳密追跡から現れた自己修復・衣の形成・3 軸凝縮 ―",
    "木原 範昭（株式会社WFシステム）　ORCID 0009-0004-6753-4020",
    "2026年11月28日",
    "2026年11月28日　Version DOI %s",
    "ja.md",
    "figs_ja",
)
prep(
    "triangle_loop_four_stage_transitions_probe_selection_rapid2_en.md",
    "Supercooled Time-Crystal Background, Probe Selection, and Four-Stage Doubling Entanglement Phase Transitions in the 3D Triangle-Loop Product (Second Rapid Communication)",
    "― Self-Repair, Dressed-Cloud Formation, and 3-Axis Condensation Revealed by Exact τ = 10,000 Tracking of N = 64 (O8 ⊗ O8) ―",
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
    # ⊗ や τ がテキスト中に裸で残った場合に安全に数式化
    t = t.replace("⊗", "$\\otimes$")
    p.write_text(
        t.replace("\\begin{document}", "\\ltjsetparameter{jacharrange={-2,-3}}\n\\begin{document}", 1),
        encoding="utf-8",
    )
print("rapid2 tex ready")
