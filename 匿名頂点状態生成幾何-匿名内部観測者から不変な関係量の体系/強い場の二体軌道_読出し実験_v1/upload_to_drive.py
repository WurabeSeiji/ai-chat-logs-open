#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
upload_to_drive.py  —  読出し実験の図と JSON を Google Drive に保存する（標準ライブラリのみ）

保存するもの（このフォルダの中のファイルをそのまま送る）
  figures/readout_summary_ja|en.png|svg, figures/readout_details_ja|en.png|svg  → Drive の figures/（無ければ作る）
  results/readout_results.json                                                    → Drive の results/（既存）
保存先: Drive「強い場の二体軌道_読出し実験_v1」 folderId 1N6RnwSGlbjgRDXdo3FvhhCOd7PELmMBH

使い方
  python upload_to_drive.py --dry-run                          # 送るファイルの一覧だけ表示
  GOOGLE_ACCESS_TOKEN=<アクセストークン> python upload_to_drive.py
  （--token <アクセストークン> でも可。同じ名前のファイルが Drive にあれば送らない。--overwrite で上書き）

アクセストークンには Drive のスコープ https://www.googleapis.com/auth/drive が必要（有効期限はふつう1時間）。
送った後、Drive 上のサイズとこちらのサイズが一致するかを表示する。
"""
import argparse
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT_ID = "1N6RnwSGlbjgRDXdo3FvhhCOd7PELmMBH"      # 強い場の二体軌道_読出し実験_v1
RESULTS_ID = "1wbogr7d80Mp_sVdUNJls7BXQPNxbr3nO"   # その中の results/
API = "https://www.googleapis.com/drive/v3/files"
UPLOAD = "https://www.googleapis.com/upload/drive/v3/files"
FOLDER = "application/vnd.google-apps.folder"
MIME = {".png": "image/png", ".svg": "image/svg+xml", ".json": "application/json"}


def call(method, url, token, data=None, headers=None):
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", "Bearer " + token)
    for k, v in (headers or {}).items():
        req.add_header(k, v)
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")[:600]
        reason = e.headers.get("x-deny-reason")
        if reason:
            hint = "（実行環境のネットワーク設定で www.googleapis.com が許可されていない可能性）"
        elif e.code == 401:
            hint = "（アクセストークンが無効か期限切れ。新しいトークンが必要）"
        elif e.code in (403, 404):
            hint = "（トークンのスコープ不足か、保存先フォルダへの権限がない可能性）"
        else:
            hint = ""
        sys.exit(f"HTTP {e.code} {method} {url.split('?')[0]} {hint}\n  x-deny-reason: {reason}\n  {body}")
    except urllib.error.URLError as e:
        sys.exit(f"接続できません: {e.reason}（ネットワーク設定を確認）")


def find(token, parent, name, mime=None):
    q = f"'{parent}' in parents and name = '{name}' and trashed = false"
    if mime:
        q += f" and mimeType = '{mime}'"
    url = API + "?" + urllib.parse.urlencode({"q": q, "fields": "files(id,name,size)", "pageSize": "10"})
    return call("GET", url, token).get("files", [])


def ensure_folder(token, parent, name):
    hit = find(token, parent, name, FOLDER)
    if hit:
        return hit[0]["id"], "既存"
    body = json.dumps({"name": name, "mimeType": FOLDER, "parents": [parent]}).encode("utf-8")
    res = call("POST", API + "?fields=id", token, body, {"Content-Type": "application/json; charset=UTF-8"})
    return res["id"], "新規作成"


def upload(token, parent, path, overwrite):
    existing = find(token, parent, path.name)
    if existing and not overwrite:
        return existing[0], "既存のため送らない"
    boundary = "=====" + uuid.uuid4().hex
    meta = {"name": path.name} if existing else {"name": path.name, "parents": [parent]}
    mime = MIME.get(path.suffix.lower(), "application/octet-stream")
    body = (f"--{boundary}\r\nContent-Type: application/json; charset=UTF-8\r\n\r\n"
            f"{json.dumps(meta)}\r\n--{boundary}\r\nContent-Type: {mime}\r\n\r\n").encode("utf-8")
    body += path.read_bytes() + f"\r\n--{boundary}--\r\n".encode("utf-8")
    headers = {"Content-Type": f"multipart/related; boundary={boundary}"}
    if existing:
        return call("PATCH", f"{UPLOAD}/{existing[0]['id']}?uploadType=multipart&fields=id,name,size",
                    token, body, headers), "上書き"
    return call("POST", UPLOAD + "?uploadType=multipart&fields=id,name,size", token, body, headers), "新規"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--token", default=os.environ.get("GOOGLE_ACCESS_TOKEN"))
    ap.add_argument("--dry-run", action="store_true", help="送るファイルの一覧だけ表示する")
    ap.add_argument("--overwrite", action="store_true", help="同じ名前のファイルがあれば上書きする")
    args = ap.parse_args()

    plan = [("figures", HERE / "figures" / f"readout_{kind}_{lang}.{ext}")
            for kind in ("summary", "details") for lang in ("ja", "en") for ext in ("png", "svg")]
    plan.append(("results", HERE / "results" / "readout_results.json"))
    missing = [str(p.relative_to(HERE)) for _, p in plan if not p.exists()]
    if missing:
        sys.exit("見つからないファイル: " + ", ".join(missing))
    print(f"送るファイル（{len(plan)} 本）:")
    for dest, p in plan:
        print(f"  {dest}/{p.name}  {p.stat().st_size} bytes")
    if args.dry_run:
        return
    if not args.token:
        sys.exit("アクセストークンがありません（GOOGLE_ACCESS_TOKEN または --token）")
    fig_id, how = ensure_folder(args.token, ROOT_ID, "figures")
    print(f"Drive の figures/: {how}（folderId {fig_id}）")
    all_ok = True
    for dest, p in plan:
        info, how = upload(args.token, fig_id if dest == "figures" else RESULTS_ID, p, args.overwrite)
        size = int(info.get("size", -1))
        ok = size == p.stat().st_size
        all_ok &= ok
        print(f"  {how}: {dest}/{p.name}  id {info['id']}  Drive {size} / ローカル {p.stat().st_size} bytes"
              f"  {'一致' if ok else '不一致'}")
    print("すべてのファイルのサイズが一致しました。" if all_ok else "サイズが一致しないファイルがあります。")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
