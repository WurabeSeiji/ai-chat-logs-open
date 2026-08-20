#!/usr/bin/env python3
"""論文目次ジェネレータ
Zenodo API（ORCID 0009-0004-6753-4020）から最新版レコードを取得し、
  PAPERS_ja.md / PAPERS_en.md      … 人間向け目次（リポジトリルート）
  docs/index.html, docs/p/<id>.html … Google Scholar 向け（citation_* メタタグ付き、GitHub Pages）
を生成する。リポジトリ内の md を DOI で照合し、日本語題・フォルダ・RELEASE_NOTES を紐付ける。

使い方:  python3 tools/gen_paper_index.py        （ZENODO_TOKEN があれば 100 件/頁、無ければ 25 件/頁）
公開シーケンス末尾で実行 → git add PAPERS_*.md docs → commit/push
"""
import json, os, re, sys, html, subprocess, collections, urllib.request, urllib.parse, unicodedata
NFC = lambda x: unicodedata.normalize('NFC', x)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
ORCID = "0009-0004-6753-4020"
GH = "https://github.com/WurabeSeiji/ai-chat-logs-open"
GH_BLOB = GH + "/blob/main/"
GH_TREE = GH + "/tree/main/"
PAGES = "https://wurabeseiji.github.io/ai-chat-logs-open"
TOKEN = os.environ.get("ZENODO_TOKEN", "")

SERIES_EN = {
    "グノモン正写像による4次元時空の幾何学的定式化": "Central Projection Framework (Gnomonic Spacetime Geometry)",
    "中心投影による宇宙の3層モデル": "Three-Layer Projection Model of the Universe",
    "波動方程式": "Wave Equation / Phase Equation & Delay-Circuit Model",
    "BH熱力学プログラム": "Black-Hole Thermodynamics Programme (incl. α-Identity)",
    "波長空間と周波数空間の双対幾何": "Dual Geometry of Wavelength and Frequency Space",
    "平方数を基本量とした場合の検討": "Square-Quantity Readout Series",
    "矩形ソリトンの面積保存と交換子": "Rectangular Solitons: Area Conservation and Commutators",
    "波の情報読出し": "Information Readout of Waves (Closed Two-Channel Systems)",
    "次元の生成構造": "Generative Structure of Dimensions",
    "時間軸Q軸とフェルミオンの生成構造": "Time Axis, Q Axis, and Fermion Generation",
    "新版量子論の基礎": "Reading Notes on Foundations of Quantum Theory (Shimizu)",
    "ヒッグス波についての考察": "Considerations on the Higgs Wave",
    "-": "Other / Not stored in this repository",
}
SERIES_ORDER = list(SERIES_EN.keys())
EXCLUDE_TITLE = re.compile(r"^(Withdrawal Notice|Repository Snapshot)", re.I)
CJK = re.compile(r"[぀-ヿ一-鿿]")


def fetch_all():
    size = 100 if TOKEN else 25
    hits, page = [], 1
    q = urllib.parse.quote(f'metadata.creators.person_or_org.identifiers.identifier:"{ORCID}"')
    while True:
        url = f"https://zenodo.org/api/records?q={q}&size={size}&page={page}&sort=newest"
        req = urllib.request.Request(url, headers={"Authorization": f"Bearer {TOKEN}"} if TOKEN else {})
        d = json.load(urllib.request.urlopen(req))
        h = d["hits"]["hits"]
        hits += h
        if len(h) < size:
            break
        page += 1
    return hits


def head(f, n=40):
    try:
        with open(f, encoding="utf-8", errors="ignore") as fp:
            return [next(fp) for _ in range(n)]
    except (StopIteration, OSError):
        try:
            return open(f, encoding="utf-8", errors="ignore").read().split("\n")
        except OSError:
            return []


def build_local_index():
    out = subprocess.run("grep -rIl --include='*.md' 'zenodo\\.' . | grep -v '^./.git/'",
                         shell=True, capture_output=True, text=True).stdout.split("\n")
    idx = {}
    for f in filter(None, out):
        try:
            t = open(f, encoding="utf-8", errors="ignore").read()
        except OSError:
            continue
        for m in set(re.findall(r"zenodo\.(\d{7,8})", t)):
            idx.setdefault(m, []).append(f)
    return idx


def strip_tags(s):
    return re.sub(r"<[^>]+>", "", html.unescape(s or "")).strip()


def build_name_index():
    out = subprocess.run("find . -name '*.md' -not -path './.git/*' -not -path './articles/*'",
                         shell=True, capture_output=True, text=True).stdout.split("\n")
    idx = {}
    for f in filter(None, out):
        idx.setdefault(NFC(os.path.basename(f)), []).append(f)
    return idx


def h1_of(lines):
    for l in lines:
        if l.startswith("# "):
            return NFC(l[2:].strip())
    return ""


def fetch_text(url):
    try:
        return urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "ignore").split("\n")[:40]
    except Exception:
        return []


def resolve(h, idx, nidx):
    m = h["metadata"]
    rid = str(h["id"])
    files = h.get("files", [])
    keys = [f["key"] for f in files]
    mds = [k for k in keys if k.lower().endswith(".md")]
    is_ja = lambda k: bool(CJK.search(k)) or bool(re.search(r"(_|\b)ja(\b|_|\.)", k))
    # --- local files by Zenodo filename (authoritative) ---
    local = []
    for k in mds:
        local += nidx.get(NFC(k), [])
    local = [f for f in local if f.count("/") >= 2]          # exclude root-level files
    # fallback: DOI mention in first 40 lines, excluding root/portal files
    if not local:
        hits = [f for f in idx.get(rid, []) if f.count("/") >= 2 and not f.startswith("./articles/")
                and "note_portal" not in f and "handout" not in f and "作業一覧" not in f and "概要" not in f]
        local = [f for f in hits if any(rid in l for l in head(f))]
    folder = os.path.dirname(collections.Counter(local).most_common(1)[0][0])[2:] if local else ""
    series = NFC(folder.split("/")[0]) if folder else "-"
    if series not in SERIES_EN:
        SERIES_EN[series] = series; SERIES_ORDER.insert(-1, series)
    # --- Japanese title: H1 of the Japanese md (local first, else Zenodo) ---
    ja_title = ""
    ja_mds = [k for k in mds if is_ja(k)] or [k for k in mds if not re.search(r"_en", k)]
    for k in ja_mds:
        for f in nidx.get(NFC(k), []):
            ja_title = h1_of(head(f))
            if ja_title: break
        if ja_title: break
    if not ja_title:
        for k in ja_mds:
            ja_title = h1_of(fetch_text(f"https://zenodo.org/api/records/{rid}/files/{urllib.parse.quote(k)}/content"))
            if ja_title: break
    if ja_title and not CJK.search(ja_title):
        ja_title = ""
    # --- nearest RELEASE_NOTES mentioning this record ---
    rn = ""
    d = folder
    while d:
        rns = sorted(x for x in os.listdir(d) if x.startswith("RELEASE_NOTES") and x.endswith(".md"))
        rns_hit = [x for x in rns if any(rid in l for l in open(os.path.join(d, x), encoding="utf-8", errors="ignore"))]
        if rns_hit:
            rn = os.path.join(d, rns_hit[0]); break
        if len(rns) == 1 and d == series:
            rn = os.path.join(d, rns[0]); break
        d = os.path.dirname(d)
    pdfs = [k for k in keys if k.lower().endswith(".pdf")]
    ja_pdf = next((p for p in pdfs if is_ja(p)), None)
    en_pdf = next((p for p in pdfs if p != ja_pdf), None)
    pdf_url = lambda k: f"https://zenodo.org/records/{rid}/files/{urllib.parse.quote(k)}" if k else ""
    return dict(
        id=rid, concept=str(h.get("conceptrecid", "")), doi=m["doi"],
        concept_doi=m.get("concept_doi") or (f"10.5281/zenodo.{h['conceptrecid']}" if h.get("conceptrecid") else ""),
        title=NFC(m["title"]), ja_title=ja_title or (NFC(m["title"]) if CJK.search(m["title"]) else ""),
        date=m.get("publication_date", ""), version=str(m.get("version") or ""),
        abstract=strip_tags(m.get("description", "")), keywords=m.get("keywords", []),
        series=series, folder=folder, release_notes=rn,
        en_pdf=pdf_url(en_pdf), ja_pdf=pdf_url(ja_pdf),
    )


def gh_link(path):
    return urllib.parse.quote(path, safe="/")


def md_index(recs, lang):
    ja = lang == "ja"
    L = []
    if ja:
        L += ["# 論文目次（全公開論文・最新版）", "",
              f"木原 範昭（Noriaki Kihara）/ WF System Co., Ltd. / ORCID [{ORCID}](https://orcid.org/{ORCID})  ",
              "全論文は Zenodo で CC BY 4.0 公開。各行：題名（日本語 / 英語）・Version DOI・Concept DOI・公開日・版・PDF・リポジトリ内フォルダ・リリースノート。",
              f"本ファイルは `tools/gen_paper_index.py` により Zenodo API から自動生成。English: [PAPERS_en.md](PAPERS_en.md) / Scholar 用ページ: {PAGES}/", ""]
    else:
        L += ["# Index of Papers (all public papers, latest versions)", "",
              f"Noriaki Kihara / WF System Co., Ltd. / ORCID [{ORCID}](https://orcid.org/{ORCID})  ",
              "All papers are published on Zenodo under CC BY 4.0. Columns: title (English / Japanese), Version DOI, Concept DOI, date, version, PDFs, folder in this repository, release notes.",
              f"Auto-generated from the Zenodo API by `tools/gen_paper_index.py`. 日本語: [PAPERS_ja.md](PAPERS_ja.md) / Scholar landing pages: {PAGES}/", ""]
    L += [f"**{'論文数' if ja else 'Papers'}: {len(recs)}**", ""]
    by = collections.defaultdict(list)
    for r in recs:
        by[r["series"]].append(r)
    for s in SERIES_ORDER:
        rs = sorted(by.get(s, []), key=lambda r: (r["date"], r["id"]))
        if not rs:
            continue
        if s == '-':
            L += ["## " + ("その他（本リポジトリ未収録） / Other (not stored in this repository)" if ja else "Other (not stored in this repository) / その他（本リポジトリ未収録）"), ""]
        else:
            L += [f"## {s} / {SERIES_EN[s]}" if ja else f"## {SERIES_EN[s]} / {s}", ""]
        for i, r in enumerate(rs, 1):
            t1 = (r["ja_title"] or r["title"]) if ja else r["title"]
            t2 = (r["title"] if r["ja_title"] and r["ja_title"] != r["title"] else "") if ja else (r["ja_title"] if r["ja_title"] != r["title"] else "")
            line = f"{i}. **{t1}**"
            if t2:
                line += f"  \n   {t2}"
            line += f"  \n   DOI [{r['doi']}](https://doi.org/{r['doi']})"
            if r["concept_doi"]:
                line += f" · Concept [{r['concept_doi']}](https://doi.org/{r['concept_doi']})"
            line += f" · {r['date']}" + (f" · v{r['version'].lstrip('v')}" if r["version"] else "")
            pdfs = []
            if r["en_pdf"]:
                pdfs.append(f"[PDF en]({r['en_pdf']})")
            if r["ja_pdf"]:
                pdfs.append(f"[PDF ja]({r['ja_pdf']})")
            if pdfs:
                line += " · " + " ".join(pdfs)
            if r["folder"]:
                line += f"  \n   {'フォルダ' if ja else 'folder'}: [{r['folder']}]({GH_TREE}{gh_link(r['folder'])})"
            if r["release_notes"]:
                line += f" · [{'リリースノート' if ja else 'release notes'}]({GH_BLOB}{gh_link(r['release_notes'])})"
            L.append(line)
        L.append("")
    return "\n".join(L)


def page_html(r):
    e = html.escape
    metas = [("citation_title", r["title"]), ("citation_author", "Kihara, Noriaki"),
             ("citation_publication_date", r["date"].replace("-", "/")), ("citation_doi", r["doi"]),
             ("citation_publisher", "Zenodo"), ("citation_technical_report_institution", "WF System Co., Ltd."),
             ("citation_language", "en"), ("citation_abstract_html_url", f"https://doi.org/{r['doi']}")]
    if r["en_pdf"]:
        metas.append(("citation_pdf_url", r["en_pdf"]))
    for k in r["keywords"]:
        metas.append(("citation_keywords", k))
    mt = "\n".join(f'<meta name="{e(k)}" content="{e(v)}">' for k, v in metas)
    links = [f'<a href="https://doi.org/{e(r["doi"])}">Zenodo record (DOI {e(r["doi"])})</a>']
    if r["concept_doi"]:
        links.append(f'<a href="https://doi.org/{e(r["concept_doi"])}">Concept DOI (all versions)</a>')
    if r["en_pdf"]:
        links.append(f'<a href="{e(r["en_pdf"])}">PDF (English)</a>')
    if r["ja_pdf"]:
        links.append(f'<a href="{e(r["ja_pdf"])}">PDF (日本語)</a>')
    if r["folder"]:
        links.append(f'<a href="{GH_TREE}{gh_link(r["folder"])}">Source folder on GitHub</a>')
    if r["release_notes"]:
        links.append(f'<a href="{GH_BLOB}{gh_link(r["release_notes"])}">Release notes</a>')
    ja = f"<p class='ja'>{e(r['ja_title'])}</p>" if r["ja_title"] and r["ja_title"] != r["title"] else ""
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{e(r['title'])}</title>
{mt}
<link rel="canonical" href="{PAGES}/p/{r['id']}.html">
<style>body{{font-family:Georgia,serif;max-width:46em;margin:2em auto;padding:0 1em;line-height:1.55}}h1{{font-size:1.4em}}.ja{{color:#555}}.meta{{color:#333}}ul{{padding-left:1.2em}}</style>
</head><body>
<p><a href="../index.html">← Index of papers</a></p>
<h1>{e(r['title'])}</h1>{ja}
<p class="meta">Noriaki Kihara (WF System Co., Ltd.) · ORCID <a href="https://orcid.org/{ORCID}">{ORCID}</a> · {e(r['date'])}{(' · v' + e(r['version'].lstrip('v'))) if r['version'] else ''} · Series: {e(SERIES_EN[r['series']])}</p>
<h2>Abstract</h2>
<p>{e(r['abstract'])}</p>
<h2>Links</h2>
<ul>{''.join(f'<li>{l}</li>' for l in links)}</ul>
<p class="meta">License: CC BY 4.0. Published on Zenodo; this page is a landing page for indexing.</p>
</body></html>
"""


def index_html(recs):
    e = html.escape
    by = collections.defaultdict(list)
    for r in recs:
        by[r["series"]].append(r)
    secs = []
    for s in SERIES_ORDER:
        rs = sorted(by.get(s, []), key=lambda r: (r["date"], r["id"]))
        if not rs:
            continue
        items = "".join(
            f'<li><a href="p/{r["id"]}.html">{e(r["title"])}</a>'
            + (f'<br><span class="ja">{e(r["ja_title"])}</span>' if r["ja_title"] and r["ja_title"] != r["title"] else "")
            + f'<br><span class="meta">{e(r["date"])} · <a href="https://doi.org/{e(r["doi"])}">{e(r["doi"])}</a></span></li>'
            for r in rs)
        secs.append(f"<h2>{e(SERIES_EN[s])}" + (f" <span class='ja'>/ {e(s)}</span>" if s != '-' else "") + f"</h2><ol>{items}</ol>")
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Noriaki Kihara — Papers</title>
<meta name="description" content="Index of {len(recs)} papers by Noriaki Kihara (WF System Co., Ltd.) published on Zenodo: central projection geometry, closed wave systems, generative structure of dimensions, α-identity.">
<style>body{{font-family:Georgia,serif;max-width:52em;margin:2em auto;padding:0 1em;line-height:1.5}}.ja{{color:#555}}.meta{{color:#666;font-size:.9em}}li{{margin:.5em 0}}</style>
</head><body>
<h1>Noriaki Kihara <span class="ja">木原 範昭</span></h1>
<p><b>Independent researcher · WF System Co., Ltd. (Nara, Japan)</b><br>
B.Eng. Osaka University (Mechanical Engineering, 1983). Forty years as a computer-graphics programmer, software entrepreneur and management consultant; since 2026 working full-time on the foundations of physics.</p>
<p><b>Research.</b> Reconstruction of the axioms of quantum theory and of spacetime geometry from a smaller set of wave/signal-theoretic principles: <i>central (gnomonic) projection</i> as the generator of 4-dimensional geometry, <i>closed two-channel wave systems</i> as a minimal model in which Born-type weights, spin/statistics, an inverse-square law and a three-direction space emerge, and a self-consistent geometric identity for the fine-structure constant (α⁻¹ = 137 + (π²/2)α). All results are published as numbered, reproducible preprints with code.</p>
<p class="meta">ORCID <a href="https://orcid.org/{ORCID}">{ORCID}</a> · <a href="https://scholar.google.com/citations?user=XUpD3aYAAAAJ">Google Scholar</a> · <a href="https://zenodo.org/search?q=metadata.creators.person_or_org.identifiers.identifier%3A%22{ORCID}%22&l=list&p=1&s=10&sort=newest">Zenodo</a> · <a href="{GH}">GitHub</a> · <a href="https://zenn.dev/noriaki_kihara">Zenn (articles, ja)</a> · <a href="https://note.com/kiharanoriaki">note (ja/en)</a> · <a href="{GH_BLOB}PAPERS_en.md">PAPERS_en.md</a> / <a href="{GH_BLOB}PAPERS_ja.md">PAPERS_ja.md</a></p>
<h2 id="papers">Papers</h2>
<p>{len(recs)} papers, latest versions, all on Zenodo under CC BY 4.0. Each entry links to a landing page with abstract, DOI, PDFs, source folder and release notes.</p>
{''.join(secs)}
</body></html>
"""


def main():
    hits = fetch_all()
    idx = build_local_index()
    nidx = build_name_index()
    recs = [resolve(h, idx, nidx) for h in hits if not EXCLUDE_TITLE.search(h["metadata"]["title"])]
    print(f"records: {len(hits)}  included: {len(recs)}", file=sys.stderr)
    open("PAPERS_ja.md", "w", encoding="utf-8").write(md_index(recs, "ja"))
    open("PAPERS_en.md", "w", encoding="utf-8").write(md_index(recs, "en"))
    os.makedirs("docs/p", exist_ok=True)
    open("docs/.nojekyll", "w").close()
    open("docs/index.html", "w", encoding="utf-8").write(index_html(recs))
    for r in recs:
        open(f"docs/p/{r['id']}.html", "w", encoding="utf-8").write(page_html(r))
    json.dump(recs, open("docs/papers.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    for r in recs:
        if not r["ja_title"] or r["series"] == "-":
            print(f"  check: {r['id']} series={r['series']} ja_title={'ok' if r['ja_title'] else 'MISSING'} {r['title'][:60]}", file=sys.stderr)


if __name__ == "__main__":
    main()
