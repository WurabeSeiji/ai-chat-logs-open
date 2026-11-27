#!/usr/bin/env python3
"""arXiv の LaTeXML HTML を Markdown に変換する。
使い方: python3 -I html2md.py paper.html out.md figs_dirname
数式は LaTeXML の TeX 注釈をそのまま $...$ / $$...$$ に入れる。
"""
import re
import sys
from bs4 import BeautifulSoup, NavigableString, Tag

src, dst, figdir = sys.argv[1], sys.argv[2], sys.argv[3]
soup = BeautifulSoup(open(src, encoding="utf-8").read(), "lxml")
art = soup.find("article")

footnotes = []          # (番号, 本文)
out = []                # ブロックの列
stats = dict(inline_math=0, display_eq=0, tagged=0, figures=0, footnotes=0, bib=0, paras=0, heads=0)


def cls(t):
    return t.get("class", []) if isinstance(t, Tag) else []


def tex_of(m):
    ann = m.find("annotation", attrs={"encoding": "application/x-tex"})
    t = ann.get_text() if ann else m.get("alttext", "")
    t = t.replace("\n", " ")
    t = t.replace("\\displaystyle}{", "")     # LaTeXML が複数行の式の継ぎ目に残す断片
    t = re.sub(r"%\s*", "", t) if "%" in t and "\\%" not in t else t
    t = re.sub(r"\s+", " ", t).strip()
    return re.sub(r"(?<!\\)\\$", "", t).strip()      # 末尾の「\ 」が削れて残る \ を除く


def esc(s):
    s = s.replace("\u2061", "").replace("\u2062", "").replace("\u2063", "").replace("\u200b", "")
    s = s.replace("\xa0", " ")
    s = re.sub(r"\s+", " ", s)
    return re.sub(r"([*_`<$\\])", r"\\\1", s)


def eq_rows(table):
    """数式の表から (tex, tag) の列を作る。"""
    rows = []
    row_tags = ("tr", "span")
    for tr in table.find_all(lambda t: t.name in row_tags and "ltx_eqn_row" in cls(t)):
        maths = [tex_of(m) for m in tr.find_all("math")]
        maths = [m for m in maths if m]
        tagel = tr.find(lambda t: "ltx_tag_equation" in cls(t) or "ltx_tag_equationgroup" in cls(t))
        tag = tagel.get_text().strip().strip("()") if tagel else None
        if maths or tag:
            rows.append((" ".join(maths), tag))
    return rows


def inline(node, in_note=False):
    """行内要素を Markdown の文字列にする。"""
    if isinstance(node, NavigableString):
        return esc(str(node))
    c = cls(node)
    n = node.name
    if n == "math":
        t = tex_of(node)
        if "pgfpicture" in t or "lxSVG" in t:      # 図の TikZ が数式として入っているもの。図は別に画像で付く
            return ""
        stats["inline_math"] += 1
        if node.get("display") == "block":
            return f" $\\displaystyle {t}$ "
        return f"${t}$"
    if n in ("svg", "script", "style", "annotation"):
        return ""
    if "ltx_note" in c and "ltx_role_footnote" in c:
        content = node.find(class_="ltx_note_content")
        mark = node.find(class_="ltx_note_mark")
        num = mark.get_text().strip() if mark else str(len(footnotes) + 1)
        if content is not None:
            for x in content.find_all(class_=["ltx_note_mark", "ltx_tag_note"]):
                x.extract()
            body = "".join(inline(ch, True) for ch in content.children).strip()
            footnotes.append((num, body))
            stats["footnotes"] += 1
        return f"[^{num}]"
    if "ltx_note" in c:      # thanks など
        return ""
    if "ltx_equation" in c or "ltx_equationgroup" in c:   # 脚注の中の別行立て数式
        parts = []
        for tex, tag in eq_rows(node):
            stats["display_eq"] += 1
            if tag:
                stats["tagged"] += 1
            parts.append(f" $\\displaystyle {tex}$" + (f" ({tag})" if tag else ""))
        return " ".join(parts) + " "
    if n == "a":
        href = node.get("href", "")
        txt = "".join(inline(ch, in_note) for ch in node.children).strip()
        if href.startswith("http"):
            return f"[{txt}]({href})" if txt else href
        return txt
    if n == "br":
        return " "
    inner = "".join(inline(ch, in_note) for ch in node.children)
    if n in ("em", "i") or "ltx_font_italic" in c or "ltx_emph" in c:
        s = inner.strip()
        return f"*{s}*" if s and "$" not in s[:1] else inner
    if n in ("b", "strong") or "ltx_font_bold" in c:
        s = inner.strip()
        return f"**{s}**" if s else inner
    if n == "sup":
        return f"^{inner}"
    return inner


def clean(s):
    s = re.sub(r"[ \t]+", " ", s).strip()
    s = re.sub(r"\s+([,.;:)])", r"\1", s)
    return s


def emit(s):
    s = s.strip()
    if s:
        out.append(s)


def intertext_of(tr):
    """数式の行に地の文(\\intertext)が入っていれば、その文を返す。"""
    import copy
    for td in tr.find_all(["td", "span"], recursive=False) or tr.find_all("td"):
        if "ltx_eqn_eqno" in cls(td):
            continue
        t2 = copy.copy(td)
        for m in t2.find_all("math"):
            m.decompose()
        if re.sub(r"[\s\u2061-\u2063]", "", t2.get_text()):
            return clean("".join(inline(ch) for ch in td.children))
    return None


def emit_equations(table):
    rows = table.find_all(lambda t: t.name in ("tr", "span") and "ltx_eqn_row" in cls(t))
    for tr in rows:
        txt = intertext_of(tr)
        if txt is not None:
            stats["intertext"] = stats.get("intertext", 0) + 1
            emit(txt)
            continue
        for tex, tag in eq_rows_one(tr):
            stats["display_eq"] += 1
            if tag:
                stats["tagged"] += 1
                emit(f"$$\n{tex} \\tag{{{tag}}}\n$$")
            else:
                emit(f"$$\n{tex}\n$$")


def eq_rows_one(tr):
    maths = [tex_of(m) for m in tr.find_all("math")]
    maths = [m for m in maths if m]
    tagel = tr.find(lambda t: "ltx_tag_equation" in cls(t) or "ltx_tag_equationgroup" in cls(t))
    tag = tagel.get_text().strip().strip("()") if tagel else None
    return [(" ".join(maths), tag)] if (maths or tag) else []


def emit_equations_old(table):
    for tex, tag in eq_rows(table):
        stats["display_eq"] += 1
        if tag:
            stats["tagged"] += 1
            emit(f"$$\n{tex} \\tag{{{tag}}}\n$$")
        else:
            emit(f"$$\n{tex}\n$$")


def emit_figure(fig):
    tagel = fig.find(class_="ltx_tag_figure")
    m = re.search(r"(\d+)", tagel.get_text()) if tagel else None
    if not m:
        for ch in fig.children:
            block(ch)
        return
    num = int(m.group(1))
    stats["figures"] += 1
    caps = []
    for cap in fig.find_all("figcaption"):
        t = clean("".join(inline(ch) for ch in cap.children))
        if t:
            caps.append(t)
    main = [c for c in caps if c.startswith("Figure")]
    subs = [c for c in caps if not c.startswith("Figure")]
    emit(f"![Figure {num}]({figdir}/fig{num:02d}.png)")
    cap = " ".join(main)
    cap = re.sub(r"^Figure (\d+):", r"**Figure \1:**", cap)
    if subs:
        cap += " (" + " / ".join(subs) + ")"
    emit(cap)


def emit_table(tab):
    rows = []
    for tr in tab.find_all("tr"):
        cells = [clean("".join(inline(ch) for ch in td.children)) for td in tr.find_all(["td", "th"])]
        rows.append(cells)
    if not rows:
        return
    w = max(len(r) for r in rows)
    lines = []
    for i, r in enumerate(rows):
        r = r + [""] * (w - len(r))
        lines.append("| " + " | ".join(x.replace("|", "\\|") for x in r) + " |")
        if i == 0:
            lines.append("|" + "---|" * w)
    emit("\n".join(lines))


def block(node):
    if isinstance(node, NavigableString):
        if str(node).strip():
            emit(esc(str(node)))
        return
    c = cls(node)
    n = node.name
    if n == "nav" or "ltx_TOC" in c or "ltx_rdf" in c or "ltx_pagination" in c:
        return
    if n in ("h1", "h2", "h3", "h4", "h5", "h6"):
        if "ltx_title_document" in c:
            emit("# " + clean("".join(inline(ch) for ch in node.children)))
        elif "ltx_title_abstract" in c:
            emit("## Abstract")
        else:
            level = {"h2": 2, "h3": 3, "h4": 4, "h5": 5, "h6": 6}.get(n, 2)
            tag = node.find(class_="ltx_tag")
            tagtxt = ""
            if tag is not None:
                tagtxt = clean(tag.get_text())
                tag.extract()
            title = clean("".join(inline(ch) for ch in node.children))
            emit("#" * level + " " + (tagtxt + " " if tagtxt else "") + title)
        stats["heads"] += 1
        return
    if "ltx_authors" in c:
        name = node.find(class_="ltx_personname")
        aff = node.find(class_="ltx_role_affiliation")
        if name:
            emit("**" + clean(name.get_text()) + "**")
        if aff:
            emit("*" + re.sub(r"^Affiliation:\s*", "", clean(aff.get_text(" "))) + "*")
        return
    if "ltx_dates" in c:
        emit(clean(node.get_text(" ")).strip("()"))
        return
    if n == "figure":
        emit_figure(node)
        return
    if n == "table":
        if "ltx_eqn_table" in c:
            emit_equations(node)
        else:
            emit_table(node)
        return
    if n == "p":
        t = clean("".join(inline(ch) for ch in node.children))
        if t:
            stats["paras"] += 1
            emit(t)
        return
    if n in ("ol", "ul"):
        if "ltx_biblist" in c:
            for li in node.find_all("li", class_="ltx_bibitem"):
                key = li.find(class_="ltx_tag_bibitem")
                k = clean(key.get_text()) if key else ""
                if key:
                    key.extract()
                parts = [clean("".join(inline(ch) for ch in b.children)) for b in li.find_all(class_="ltx_bibblock")]
                parts = [p for p in parts if p]
                emit(f"- **{k}** " + " ".join(parts))
                stats["bib"] += 1
            return
        items = []
        for li in node.find_all("li", recursive=False):
            tag = li.find(class_="ltx_tag_item")
            lab = clean(tag.get_text()) if tag else "-"
            if tag:
                tag.extract()
            if not re.match(r"^\d+\.$", lab):
                lab = "-"
            sub_start = len(out)
            for ch in li.children:
                block(ch)
            body = out[sub_start:]
            del out[sub_start:]
            items.append(lab + " " + "\n\n   ".join(b.replace("\n", "\n   ") for b in body))
        emit("\n".join(items))
        return
    # 入れ物はそのまま下る
    block_children = any(isinstance(ch, Tag) and ch.name in ("p", "div", "section", "table", "figure", "ol", "ul",
                                                            "h1", "h2", "h3", "h4", "h5", "h6", "nav") for ch in node.children)
    if block_children:
        for ch in node.children:
            block(ch)
    else:
        t = clean(inline(node))
        if t:
            emit(t)


for ch in art.children:
    block(ch)

header_note = (
    "> Source: arXiv:2611.34924v1 [hep-th], submitted 25 Nov 2026. "
    "https://arxiv.org/abs/2611.34924  \n"
    "> This Markdown was converted mechanically from the arXiv HTML version (LaTeXML); "
    "the equations are the TeX of that version, and the figures are cropped from the arXiv PDF. "
    "The PDF is the authoritative text. Reference numbers [n] follow the HTML version (alphabetical) "
    "and differ from those printed in the PDF; they are consistent with the reference list below.  \n"
    "> License of the original: arXiv.org perpetual, non-exclusive license "
    "(http://arxiv.org/licenses/nonexclusive-distrib/1.0/)."
)
out.insert(1, header_note)

if footnotes:
    out.append("---")
    out.append("## Footnotes")
    for num, body in footnotes:
        out.append(f"[^{num}]: {body}")

open(dst, "w", encoding="utf-8").write("\n\n".join(out) + "\n")
print(stats)
