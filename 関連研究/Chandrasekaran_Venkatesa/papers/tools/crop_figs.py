#!/usr/bin/env python3
"""PDF から図を切り出す。使い方: python3 -I crop_figs.py paper.pdf outdir workdir
pdftotext -bbox-layout で「Figure N:」の説明文の位置を取り、
その上にある幅の広い本文ブロックの下端から説明文の上端までを図とみなす。
"""
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

pdf, outdir, work = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
outdir.mkdir(parents=True, exist_ok=True)
work.mkdir(parents=True, exist_ok=True)
DPI = 200
xml = work / "bbox.html"
subprocess.run(["pdftotext", "-bbox-layout", pdf, str(xml)], check=True)
txt = xml.read_text(encoding="utf-8")
txt = re.sub(r"<!DOCTYPE[^>]*>", "", txt)
txt = txt.replace('xmlns="http://www.w3.org/1999/xhtml"', "")
root = ET.fromstring(txt)

found = {}
for pno, page in enumerate(root.iter("page"), start=1):
    W, H = float(page.get("width")), float(page.get("height"))
    blocks = []
    for b in page.iter("block"):
        wl = list(b.iter("word"))
        words = [w.text or "" for w in wl]
        capy, capn = None, None
        for i, w in enumerate(wl[:-1]):
            if w.text == "Figure" and re.fullmatch(r"\d+:", wl[i + 1].text or ""):
                capy, capn = float(w.get("yMin")), int(wl[i + 1].text[:-1])
                break
        blocks.append(dict(x0=float(b.get("xMin")), y0=float(b.get("yMin")), x1=float(b.get("xMax")),
                           y1=float(b.get("yMax")), text=" ".join(words), capy=capy, capn=capn))
    textw = W - 2 * 72
    caps = [b for b in blocks if b["capn"] is not None]
    for cap in sorted(caps, key=lambda b: b["y0"]):
        num = cap["capn"]
        top = 60.0
        for b in blocks:
            if b is cap or b["y1"] > cap["capy"] + 1:
                continue
            wide = (b["x1"] - b["x0"]) > 0.82 * textw
            is_cap = b["capn"] is not None
            body = wide and (len(b["text"]) > 150 or is_cap)
            # 番号つきの式や短い本文行も本文として扱う
            eqno = re.search(r"\(\d+\.\d+[a-z]?\)\s*$|\([A-E]\.\d+[a-z]?\)\s*$", b["text"]) is not None and b["x1"] > W - 80
            if body or eqno:
                top = max(top, b["y1"])
        found[num] = (pno, top + 3, cap["capy"] - 3, W, H)

for num, (pno, y0, y1, W, H) in sorted(found.items()):
    s = DPI / 72.0
    x, y, w, h = int(50 * s), int(y0 * s), int((W - 100) * s), int((y1 - y0) * s)
    base = outdir / f"fig{num:02d}"
    subprocess.run(["pdftoppm", "-r", str(DPI), "-f", str(pno), "-l", str(pno), "-x", str(x), "-y", str(y),
                    "-W", str(w), "-H", str(h), "-png", "-singlefile", pdf, str(base)], check=True)
    print(num, "page", pno, "height_pt", round(y1 - y0, 1))
print("figures:", len(found))
