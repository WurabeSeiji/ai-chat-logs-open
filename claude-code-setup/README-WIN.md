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

## 前提条件

### 1. Google Drive for desktop（必須）

- **Mirror モード**で動作していること
  - Stream モード（仮想ファイルシステム）では Claude Code の頻繁な書き込みでエラーが起きやすい
  - 確認方法：Google Drive for desktop 設定 → 「マイドライブの設定」→ 「ファイルをミラーリングする」を選択
- Google アカウントが `kihara.noriaki@gmail.com`（または共有対象アカウント）でサインインされていること
- マウントポイントが正常（通常は `G:\` 仮想ドライブまたは `%USERPROFILE%\Google Drive\` のいずれか）

### 2. Developer Mode または 管理者権限（必須）

Windows の `mklink /D`（directory symlink 作成）には以下のいずれかが必要：

- **Developer Mode を有効化**（推奨、一度きりの設定で済む）
  - 設定 → プライバシーとセキュリティ → 開発者向け → 開発者モード を ON
- **管理者権限で bat を実行**（Developer Mode を有効化しない場合）
  - エクスプローラから `setup-memory-sync.bat` を右クリック → 「管理者として実行」

`mklink /D` が失敗した場合、スクリプトは自動的に `mklink /J`（junction）にフォールバックする。junction は管理者権限不要だが、対象が同一ボリュームでないと作成不可。Google Drive がローカルドライブとして Mirror モードでマウントされていれば動作する。

### 3. Claude Code（claude.ai/code）インストール済み

`%USERPROFILE%\.claude\projects\` が存在していること（Claude Code を一度でも起動していれば作成される）。

### 4. VSCode と全 Claude Code セッションが終了していること

実行中だとファイルハンドルの競合が起こる可能性があるため。

## 配置

```
G:\マイドライブ\OneDrive\                                  ← Google Drive 上
├── GitHub\
│   └── ai-chat-logs-open\
│       └── claude-code-setup\
│           ├── setup-memory-sync.sh     (macOS 版)
│           ├── setup-memory-sync.bat    ← このスクリプト
│           ├── README.md                 (macOS 版)
│           └── README-WIN.md            ← この文書
└── ClaudeCode\
    └── memory\                          ← 各プロジェクトの memory が集約される
        ├── -C--Users-USERNAME-projects-PROJECT_A\
        ├── -C--Users-USERNAME-projects-PROJECT_B\
        └── ...
```

## 実行手順

### 1. スクリプトを取得

リポジトリを `git clone` するか、Google Drive 上の `OneDrive\GitHub\ai-chat-logs-open\claude-code-setup\setup-memory-sync.bat` を直接実行する。

### 2. 設定を確認・調整

`setup-memory-sync.bat` の冒頭部分を環境に合わせて編集：

```batch
REM Google Drive のローカルパス (Mirror モード前提)
set "GDRIVE_TARGET=G:\マイドライブ\OneDrive\ClaudeCode\memory"

REM Claude Code projects ディレクトリ
set "LOCAL_PROJECTS=%USERPROFILE%\.claude\projects"
```

#### Google Drive パスのパターン例

| Google Drive 設定 | パス例 |
|---|---|
| Stream モード（仮想ドライブ G:） | `G:\マイドライブ\OneDrive\ClaudeCode\memory` |
| Mirror モード（マイドキュメント等） | `%USERPROFILE%\Google Drive\マイドライブ\OneDrive\ClaudeCode\memory` |
| 英語 UI で My Drive 表示 | `G:\My Drive\OneDrive\ClaudeCode\memory` |

エクスプローラで Google Drive を開き、`OneDrive\ClaudeCode` のフルパスをアドレスバーから確認するのが確実。

### 3. 実行

- **Developer Mode 有効の場合**：エクスプローラで `setup-memory-sync.bat` をダブルクリック
- **管理者権限が必要な場合**：右クリック → 「管理者として実行」

### 4. ログ確認

スクリプトが pre-flight チェック → 移動 → symlink 作成 → 検証 を順に行い、最後に「[OK] セットアップ完了」が表示されれば成功。

### 5. VSCode 起動・動作確認

VSCode を起動し、何かしらのプロジェクトを開いて Claude Code が正常動作するか確認。

```
Claude Code に「メモリ確認」など指示 → 既存の memory が読まれれば成功
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
  [MOVE] -C--Users-USERNAME-projects-PROJECT_A
    [OK] symlink 作成
  [MOVE] -C--Users-USERNAME-projects-PROJECT_B
    [OK] symlink 作成

  集計: 移動 2 / 既 symlink 0 / memory 無 0 / エラー 0

===============================================
  [OK] セットアップ完了
===============================================
```

## トラブルシューティング

### `mklink /D` が「アクセスが拒否されました」エラー

- Developer Mode が無効、かつ管理者権限なしで実行している
- 対処：Developer Mode を有効化するか、管理者として bat を再実行
- スクリプトは自動的に junction (`/J`) にフォールバックを試みるが、それも失敗する場合は要対処

### `mklink /J` が「対応していないネットワークパス」エラー

- Google Drive が Stream モード（仮想ファイルシステム）で動作している
- Junction はローカルボリューム同士でないと作成できない
- 対処：Google Drive を Mirror モードに切り替えてから再実行

### 「移動先に既にあり」エラー

- 過去に同様の操作を行ったため、Google Drive 側に既にデータがある
- 対処：手動で `%GDRIVE_TARGET%\<プロジェクト識別子>\` の中身を確認し、最新のものを残して旧データを削除してから再実行

### 文字化けで日本語パスが正しく扱えない

- コードページが UTF-8 (`chcp 65001`) になっていない
- 対処：スクリプト冒頭で `chcp 65001` を実行している。それでも問題があればコマンドプロンプトを「フォント → MS ゴシック等」に変更

### Claude Code が memory を読まなくなった

- symlink が壊れている可能性
- 確認：`dir /AL %USERPROFILE%\.claude\projects\<プロジェクト識別子>` で `<SYMLINKD>` または `<JUNCTION>` 表示があるか
- 復旧：下記ロールバック手順を実行してから再セットアップ

## Claude.ai 側からの参照

セットアップ後、ブラウザ Claude.ai（Google Drive コネクタ有効化済み）から：

```
「Google Drive の MEMORY.md を見て、〇〇について教えて」
「Drive で 'user_contacts_morita' を検索して内容を要約」
```

同一 Google アカウントなら、絶対パス指定は不要。ファイル名やキーワードで検索される。

## ロールバック手順

何らかの問題があった場合、以下の batch で元に戻せる：

```batch
@echo off
setlocal enabledelayedexpansion
set "GDRIVE_TARGET=G:\マイドライブ\OneDrive\ClaudeCode\memory"
set "LOCAL_PROJECTS=%USERPROFILE%\.claude\projects"

for /d %%P in ("%LOCAL_PROJECTS%\*") do (
    set "PROJ_NAME=%%~nxP"
    set "SRC=%%P\memory"
    set "DST=%GDRIVE_TARGET%\!PROJ_NAME!"

    REM symlink/junction を判定
    dir "%%P" /AL 2>nul | findstr /i "memory" >nul
    if not errorlevel 1 (
        rmdir "!SRC!"
        if exist "!DST!" (
            move "!DST!" "!SRC!" >nul
            echo Restored: !PROJ_NAME!
        )
    )
)

REM Google Drive 側の空フォルダを削除
rmdir "%GDRIVE_TARGET%" 2>nul
endlocal
```

## 注意事項

- **同一プロジェクトを 2 マシン同時に編集しない**：プロジェクト識別子が違えば自然分離されるが、同一マシン上で同一プロジェクトを並列で開くのは避ける
- **conflicted copy の監視**：稀に `*[Conflict]*.md` ファイルが Google Drive に生成される可能性。週次で確認推奨
- **Google アカウント保護**：2 要素認証必須。memory には機微情報を含むため
- **Windows VPS 特有の留意点**：
  - VPS の Google Drive クライアントが常時動作していることを確認
  - VPS 再起動時に Google Drive クライアントが自動起動するよう設定
  - 同期が止まると memory への書き込みは成功するがクラウド未反映となる
