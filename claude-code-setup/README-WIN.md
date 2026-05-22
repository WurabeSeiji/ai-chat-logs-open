# Claude Code memory → Google Drive 同期セットアップ（Windows 版）

このフォルダの `setup-memory-sync.bat` を使って、Windows マシン（VPS 含む）の Claude Code memory を Google Drive 経由で同期するためのガイド。

macOS 版と同じ設計思想で、**プロジェクト識別子による自然な isolation**（マシンごとに別フォルダに書き込み）と **memory のみの軽量同期** を実現する。

## 設計目的

3 経路からの memory アクセス：

```
GoogleDrive/.../OneDrive/ClaudeCode/memory/<project-id>/
  ├─ Claude Code（Windows VPS）       自動ロード + 書き込み（symlink 経由）
  ├─ Claude Code（macOS 等の別マシン） 別 namespace で自動ロード + 書き込み
  └─ Claude.ai（ブラウザ/モバイル）    Google Drive コネクタ経由で読み取り
```

別マシンは **プロジェクト識別子（作業ディレクトリの絶対パス由来）が必ず異なる** ため、書き込み競合は構造的に発生しない。

---

## ⚠ 重要：実テスト（2026-05-22, c:\xampp\htdocs\osj VPS）で判明した 4 件の罠

以下の制約は **配布版 bat が本来カバーすべき内容**として現在は反映済み。**この bat の編集・再生成を行う場合は、必ず以下を維持する**こと。

### 罠 1：Stream モードでは `mklink /J`（junction）が使えない

- Google Drive Stream モード（既定の `G:\` 仮想ドライブ）は **junction 作成不可**（「対応していないネットワークパス」エラー）
- → `mklink /D`（directory symlink）一択
- → `mklink /D` には **管理者権限または Developer Mode が必須**
- 本 bat は冒頭で `net session` による管理者権限チェックを行い、非管理者なら即座に exit

### 罠 2：`chcp 65001`（UTF-8 codepage）を使ってはいけない

- bat ファイル自体は **CP932（Shift-JIS）** で保存されている前提
- `chcp 65001` を実行すると：
  - 日本語 echo 出力が文字化け（CP932 バイトを UTF-8 解釈してしまう）
  - `if exist "G:\マイドライブ\..."` のパス解決が **失敗** する（実在するのに「見つからない」エラー）
- 本 bat は **`chcp` を実行せず Japanese Windows の既定 CP932 のまま**動作する

### 罠 3：ファイルエンコーディングは CP932 + CRLF を厳守

- VSCode 等で Edit すると UTF-8 + LF になりがち。cmd は CP932 で読み始めるため誤読する
- 症状例：`'o' は、内部コマンドまたは...`、`'ho.' は...`、`'_PARENT' は...`、`'eq' は...`、`&& の使い方が誤っています` 等
- **編集後は必ず CP932 + CRLF に再変換**する：

```powershell
$path = "...\setup-memory-sync.bat"
$cp932 = [System.Text.Encoding]::GetEncoding(932)
$c = [System.IO.File]::ReadAllText($path, [System.Text.Encoding]::UTF8)
$c = ($c -replace "`r`n", "`n" -replace "`r", "`n") -replace "`n", "`r`n"
[System.IO.File]::WriteAllText($path, $c, $cp932)
```

- 検証：`[System.IO.File]::ReadAllBytes($path)[0..15]` の最初の改行が `0D 0A`（CRLF）であること

### 罠 4：クロスドライブ移動には `move` ではなく `robocopy /MOVE`

- cmd の `move` コマンドは **クロスドライブ（C: → G:）のディレクトリ移動を扱えず失敗**する（「アクセスが拒否されました」）
- 内部実装は `MoveFile` API ベースで、別ボリュームのディレクトリ rename が構造的に不可
- → `robocopy "!SRC!" "!DST!" /E /MOVE /NFL /NDL /NJH /NJS /NC /NS /NP` を使用
- **重要**：robocopy の errorlevel は他コマンドと意味が違う。**0-7 が成功（警告含む）、8 以上がエラー**。`if errorlevel 1` ではなく `if errorlevel 8` で判定

---

## 前提条件

### 1. Google Drive for desktop

- **Stream モード**でも動作可（実際に c:\xampp\htdocs\osj VPS で実証済み）
- Mirror モードでも動作する（macOS 版と同等の動作）
- 同期対象に `OneDrive/ClaudeCode/memory/` フォルダが含まれていること
- マウントポイントが正常（通常は `G:\` 仮想ドライブまたは `%USERPROFILE%\Google Drive\` のいずれか）

### 2. 管理者権限が必須

bat 冒頭の `net session` チェックで非管理者なら即 exit する。**必ず「管理者として実行」する**こと。

- エクスプローラで `setup-memory-sync.bat` を右クリック → 「管理者として実行」
- もしくは管理者権限の cmd / PowerShell で起動

Developer Mode を有効化していても、`mklink /D` 経由の symlink 作成は基本的に管理者権限を要求する Windows ビルドが多いため、本 bat は管理者権限を必須としている。

### 3. Claude Code（claude.ai/code）インストール済み

`%USERPROFILE%\.claude\projects\` が存在していること（Claude Code を一度でも起動していれば作成される）。

### 4. VSCode と全 Claude Code セッションが終了していること

実行中だとファイルハンドルの競合が起こる可能性があるため。`tasklist` で `node.exe` / `Code.exe` / `claude.exe` が動作中なら警告を出す。

## 配置

```
G:\マイドライブ\OneDrive\                                  ← Google Drive 上
├── GitHub\
│   └── ai-chat-logs-open\
│       └── claude-code-setup\
│           ├── setup-memory-sync.sh     (macOS 版)
│           ├── setup-memory-sync.bat    ← このスクリプト（CP932 + CRLF）
│           ├── README.md                 (macOS 版)
│           └── README-WIN.md            ← この文書
└── ClaudeCode\
    └── memory\                          ← 各プロジェクトの memory が集約される
        ├── c--xampp-htdocs-osj\        ← Windows VPS の例
        ├── c--xampp-htdocs-fm3\
        └── -Users-...-ai-chat-logs-open\  ← macOS の例
```

## 実行手順

### 1. スクリプトを取得

リポジトリを `git clone` するか、Google Drive 上の `OneDrive\GitHub\ai-chat-logs-open\claude-code-setup\setup-memory-sync.bat` を直接実行する。

### 2. 事前確認（PowerShell）

実行前に以下を確認：

```powershell
# Google Drive のマウント状態 (Stream モードかの判定)
Get-Process | Where-Object { $_.ProcessName -match "GoogleDrive" }
# → GoogleDriveFS.exe が動作中なら Stream モード（G:\ は仮想ドライブ）

# Developer Mode の確認
Get-ItemProperty -Path "HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\AppModelUnlock" -Name "AllowDevelopmentWithoutDevLicense" -ErrorAction SilentlyContinue

# Claude Code / VSCode が動作中なら終了させてから実行
Get-Process | Where-Object { $_.ProcessName -match "^(node|Code|claude)$" }
```

### 3. 設定を確認・調整

`setup-memory-sync.bat` の冒頭部分（管理者権限チェックの直後）を環境に合わせて編集：

```batch
REM Google Drive のローカルパス
set "GDRIVE_TARGET=G:\マイドライブ\OneDrive\ClaudeCode\memory"

REM Claude Code projects ディレクトリ
set "LOCAL_PROJECTS=%USERPROFILE%\.claude\projects"
```

#### Google Drive パスのパターン例

| Google Drive 設定 | パス例 |
|---|---|
| Stream モード（既定、仮想ドライブ） | `G:\マイドライブ\OneDrive\ClaudeCode\memory` |
| Mirror モード（実体マウント） | `%USERPROFILE%\Google Drive\マイドライブ\OneDrive\ClaudeCode\memory` |
| 英語 UI で My Drive 表示 | `G:\My Drive\OneDrive\ClaudeCode\memory` |

エクスプローラで Google Drive を開き、`OneDrive\ClaudeCode` のフルパスをアドレスバーから確認するのが確実。

**編集する場合の注意**：CP932 対応エディタ（メモ帳、サクラエディタ等）を使うか、編集後に上記 PowerShell コードで CP932 + CRLF に再変換する。VSCode で編集して保存すると UTF-8 + LF になり、bat が動作しなくなる。

### 4. 実行

エクスプローラで `setup-memory-sync.bat` を右クリック → **「管理者として実行」**

非管理者で起動した場合、冒頭の `net session` チェックで「[NG] 管理者権限で実行されていません」と表示され即 exit する。

### 5. ログ確認

スクリプトが pre-flight チェック → 移動 → symlink 作成 → 検証 を順に行い、最後に「[OK] セットアップ完了」が表示されれば成功。

### 6. VSCode 起動・動作確認

VSCode を起動し、何かしらのプロジェクトを開いて Claude Code が正常動作するか確認。

```
Claude Code に「メモリ確認」など指示 → 既存の memory が読まれれば成功
```

#### PowerShell での確認（推奨）

```powershell
# 全プロジェクトの symlink タイプと Target 確認
foreach ($p in (Get-ChildItem "$env:USERPROFILE\.claude\projects" -Directory).Name) {
    $link = "$env:USERPROFILE\.claude\projects\$p\memory"
    if (Test-Path $link) {
        $item = Get-Item $link -Force
        Write-Output "$p : $($item.LinkType) -> $($item.Target)"
    }
}

# MD5 比較で symlink 経由とGoogle Drive 直接で同一ファイルか確認
$local = "$env:USERPROFILE\.claude\projects\c--xampp-htdocs-osj\memory\MEMORY.md"
$drive = "G:\マイドライブ\OneDrive\ClaudeCode\memory\c--xampp-htdocs-osj\MEMORY.md"
(Get-FileHash $local).Hash -eq (Get-FileHash $drive).Hash
```

## 期待される出力

```
===============================================
  Claude Code memory => Google Drive 同期セットアップ (Windows)
===============================================

[1/6] 設定確認
  GDRIVE_TARGET : G:\マイドライブ\OneDrive\ClaudeCode\memory
  LOCAL_PROJECTS: C:\Users\USERNAME\.claude\projects

[2/6] Google Drive マウント確認
  [OK] G:\マイドライブ\OneDrive\ClaudeCode\

[3/6] %USERPROFILE%\.claude\projects 存在確認
  [OK] C:\Users\USERNAME\.claude\projects

[4/6] Claude Code 関連プロセスの確認
  [OK] Claude Code 関連プロセス見当たらず

[5/6] Google Drive 側ターゲットディレクトリ作成
  [OK] G:\マイドライブ\OneDrive\ClaudeCode\memory

[6/6] 各プロジェクトの memory\ を移動・symlink 化
  [MOVE] c--xampp-htdocs-osj
    [OK] symlink 作成
  [MOVE] c--xampp-htdocs-fm3
    [OK] symlink 作成

  集計: 移動 2 / 既 symlink 0 / memory 無 0 / エラー 0

===============================================
  [OK] セットアップ完了
===============================================
```

## トラブルシューティング

### 「[NG] 管理者権限で実行されていません」

- bat 冒頭の `net session` チェックで検出。**右クリック → 管理者として実行**してください

### `mklink /D` が「アクセスが拒否されました」エラー

- 管理者権限なしで実行している（本 bat では冒頭でブロックされるはずなので、通常は発生しない）
- 対処：管理者として bat を再実行

### 文字化けで日本語が `??`、`縺ｭ` 等で表示される

- bat ファイルが UTF-8 で保存されている（CP932 ではない）
- 対処：上記の PowerShell スクリプトで CP932 + CRLF に再変換

### `'o' は、内部コマンドまたは...`、`&& の使い方が誤っています` 等の構文エラー連発

- 同じく bat ファイルのエンコーディング問題
- 対処：CP932 + CRLF に再変換

### `[NG] Google Drive 親ディレクトリが見つかりません`（実在するのに）

- `chcp 65001` が実行されてしまっている可能性（パス解決失敗）
- 対処：`chcp 65001` の行を削除またはコメント化

### `move` が「アクセスが拒否されました」で失敗

- bat が古い版（`move` を使用）
- 対処：本 bat の最新版に差し替え（`robocopy /MOVE` を使用）

### Claude Code が memory を読まなくなった

- symlink が壊れている可能性
- 確認：`Get-Item "$env:USERPROFILE\.claude\projects\<プロジェクト識別子>\memory" -Force` で `LinkType` と `Target` を確認
- 復旧：下記ロールバック手順を実行してから再セットアップ

## Claude.ai 側からの参照

セットアップ後、ブラウザ Claude.ai（Google Drive コネクタ有効化済み）から：

```
「Google Drive の MEMORY.md を見て、〇〇について教えて」
「Drive で 'feedback_gdrive_memory_sync_setup' を検索して内容を要約」
```

同一 Google アカウントなら、絶対パス指定は不要。ファイル名やキーワードで検索される。

## プロジェクト識別子による分離の理解

複数マシン間で同一論理プロジェクトの memory が衝突する心配は不要。識別子は **作業ディレクトリの絶対パス由来** なので：

- macOS：`/Users/.../OneDrive-GitHub-...` → `-Users-...-OneDrive-GitHub-...`
- Windows：`C:\xampp\htdocs\osj` → `c--xampp-htdocs-osj`

物理的に別フォルダになるため、書き込み競合は構造的に発生しない。`*[Conflict]*.md` 警告は同一マシンで Claude Code を並列起動した場合などに限定される話で、通常運用では発生しない。

## ロールバック手順

何らかの問題があった場合、以下の PowerShell スクリプトで元に戻せる：

```powershell
$GDRIVE_TARGET = "G:\マイドライブ\OneDrive\ClaudeCode\memory"
$LOCAL_PROJECTS = "$env:USERPROFILE\.claude\projects"

foreach ($proj in Get-ChildItem $LOCAL_PROJECTS -Directory) {
    $src = Join-Path $proj.FullName "memory"
    $dst = Join-Path $GDRIVE_TARGET $proj.Name

    if (Test-Path $src) {
        $item = Get-Item $src -Force
        if ($item.LinkType) {
            # symlink を削除
            (Get-Item $src -Force).Delete()
            # Google Drive 側から戻す
            if (Test-Path $dst) {
                robocopy $dst $src /E /MOVE /NFL /NDL /NJH /NJS /NC /NS /NP | Out-Null
                Write-Output "Restored: $($proj.Name)"
            }
        }
    }
}

# Google Drive 側の空フォルダを削除
if ((Get-ChildItem $GDRIVE_TARGET -ErrorAction SilentlyContinue).Count -eq 0) {
    Remove-Item $GDRIVE_TARGET -Force
}
```

## 注意事項

- **同一プロジェクトを 2 マシン同時に編集しない**：プロジェクト識別子が違えば自然分離されるが、同一マシン上で同一プロジェクトを並列で開くのは避ける
- **conflicted copy の監視**：稀に `*[Conflict]*.md` ファイルが Google Drive に生成される可能性。週次で確認推奨
- **Google アカウント保護**：2 要素認証必須。memory には機微情報を含むため
- **Windows VPS 特有の留意点**：
  - VPS の Google Drive クライアントが常時動作していることを確認
  - VPS 再起動時に Google Drive クライアントが自動起動するよう設定
  - 同期が止まると memory への書き込みは成功するがクラウド未反映となる
- **bat 編集時の絶対ルール**：CP932 + CRLF を維持する。Edit ツールや VSCode で保存すると UTF-8 + LF になるため、再変換 PowerShell を必ず実行
