@echo off
REM ============================================
REM Claude Code memory -> Google Drive 同期セットアップ (Windows 版)
REM
REM 目的:
REM   - 複数マシンの Claude Code 間で memory を共有
REM   - Claude.ai (ブラウザ/モバイル) からも参照可能にする
REM
REM 動作:
REM   %USERPROFILE%\.claude\projects\<id>\memory\  (実体)
REM   -> Google Drive 上に移動
REM   -> 元の場所には directory symlink を残す
REM
REM 既に symlink になっているプロジェクトはスキップ (再実行安全)
REM ============================================

setlocal enabledelayedexpansion
REM chcp は実行しない (CP932 console と CP932 ファイルで統一)

REM ============================================
REM 管理者権限チェック (mklink /D に必須)
REM Google Drive が Stream モードの場合 /J フォールバック不可のため /D 一択
REM ============================================
net session >nul 2>&1
if errorlevel 1 (
    echo.
    echo ===============================================
    echo   [NG] 管理者権限で実行されていません
    echo ===============================================
    echo.
    echo   このスクリプトは mklink /D ^(directory symlink^) を使用するため、
    echo   管理者権限が必要です。
    echo.
    echo   対処方法:
    echo     1. エクスプローラで setup-memory-sync.bat を右クリック
    echo     2. 「管理者として実行」を選択
    echo.
    echo   ※ Google Drive が Stream モードの場合、junction ^(/J^) への
    echo      フォールバックが「対応していないネットワークパス」エラーで
    echo      失敗するため、本スクリプトでは /D 一択としています。
    echo.
    pause
    exit /b 1
)

REM ============================================
REM 設定 - 環境に合わせて変更してください
REM ============================================

REM Google Drive のローカルパス (Mirror モード前提)
REM 例: G:\マイドライブ\OneDrive\ClaudeCode\memory
REM 例: %USERPROFILE%\Google Drive\マイドライブ\OneDrive\ClaudeCode\memory
set "GDRIVE_TARGET=G:\マイドライブ\OneDrive\ClaudeCode\memory"

REM Claude Code projects ディレクトリ (通常は %USERPROFILE%\.claude\projects)
set "LOCAL_PROJECTS=%USERPROFILE%\.claude\projects"

REM ============================================

echo.
echo ===============================================
echo   Claude Code memory =^> Google Drive 同期セットアップ (Windows)
echo ===============================================
echo.

echo [1/6] 設定確認
echo   GDRIVE_TARGET : %GDRIVE_TARGET%
echo   LOCAL_PROJECTS: %LOCAL_PROJECTS%
echo.

REM ===== Step 2: 親ディレクトリ存在確認 =====
echo [2/6] Google Drive マウント確認
for %%i in ("%GDRIVE_TARGET%") do set "GDRIVE_PARENT=%%~dpi"
if not exist "%GDRIVE_PARENT%" (
    echo   [NG] Google Drive 親ディレクトリが見つかりません:
    echo        %GDRIVE_PARENT%
    echo        Google Drive for desktop が起動し、Mirror モードでマウントされているか確認してください
    pause
    exit /b 1
)
echo   [OK] %GDRIVE_PARENT%
echo.

REM ===== Step 3: LOCAL_PROJECTS 存在確認 =====
echo [3/6] %%USERPROFILE%%\.claude\projects 存在確認
if not exist "%LOCAL_PROJECTS%" (
    echo   [NG] %LOCAL_PROJECTS% が存在しません
    echo        Claude Code が一度も起動していない可能性があります
    pause
    exit /b 1
)
echo   [OK] %LOCAL_PROJECTS%
echo.

REM ===== Step 4: Claude Code プロセス確認 =====
echo [4/6] Claude Code 関連プロセスの確認
set CLAUDE_RUNNING=0
tasklist /FI "IMAGENAME eq node.exe" 2>nul | findstr /i "node.exe" >nul && set CLAUDE_RUNNING=1
tasklist /FI "IMAGENAME eq Code.exe" 2>nul | findstr /i "Code.exe" >nul && set CLAUDE_RUNNING=1
tasklist /FI "IMAGENAME eq claude.exe" 2>nul | findstr /i "claude.exe" >nul && set CLAUDE_RUNNING=1

if %CLAUDE_RUNNING% equ 1 (
    echo   [!!] node.exe / Code.exe / claude.exe が動作中の可能性
    echo        VSCode と全 Claude Code セッションを終了してから続行することを推奨します
    set /p CONTINUE="   このまま続行しますか? (y/N): "
    if /i not "!CONTINUE!"=="y" (
        echo   中止しました
        pause
        exit /b 1
    )
) else (
    echo   [OK] Claude Code 関連プロセス見当たらず
)
echo.

REM ===== Step 5: ターゲットディレクトリ作成 =====
echo [5/6] Google Drive 側ターゲットディレクトリ作成
if not exist "%GDRIVE_TARGET%" (
    mkdir "%GDRIVE_TARGET%"
    if errorlevel 1 (
        echo   [NG] mkdir 失敗
        pause
        exit /b 1
    )
)
echo   [OK] %GDRIVE_TARGET%
echo.

REM ===== Step 6: 各プロジェクトの memory を移動して symlink 化 =====
echo [6/6] 各プロジェクトの memory\ を移動・symlink 化
set MOVED=0
set SKIPPED_LINK=0
set SKIPPED_NOMEM=0
set ERRORS=0

for /d %%P in ("%LOCAL_PROJECTS%\*") do (
    set "PROJ_PATH=%%P"
    set "PROJ_NAME=%%~nxP"
    set "SRC=%%P\memory"
    set "DST=%GDRIVE_TARGET%\!PROJ_NAME!"

    REM 既に symlink/junction か判定 (dir の出力に <SYMLINKD> or <JUNCTION> が含まれるか)
    set IS_LINK=0
    dir /AL "%%P" 2>nul | findstr /i "memory" >nul && set IS_LINK=1

    if !IS_LINK! equ 1 (
        echo   [SKIP] 既に symlink/junction: !PROJ_NAME!
        set /a SKIPPED_LINK+=1
    ) else if not exist "!SRC!" (
        echo   [SKIP] memory なし: !PROJ_NAME!
        set /a SKIPPED_NOMEM+=1
    ) else if exist "!DST!" (
        echo   [ERROR] 移動先に既にあり: !DST!
        echo           手動で確認してください。スキップします。
        set /a ERRORS+=1
    ) else (
        echo   [MOVE] !PROJ_NAME!
        REM cmd の move はクロスドライブのディレクトリ移動を扱えないため robocopy /MOVE を使用
        REM robocopy の errorlevel: 0-7 は成功 (警告含む)、8 以上が真のエラー
        robocopy "!SRC!" "!DST!" /E /MOVE /NFL /NDL /NJH /NJS /NC /NS /NP >nul
        if errorlevel 8 (
            echo     [NG] robocopy 失敗 ^(errorlevel=!errorlevel!^)
            set /a ERRORS+=1
        ) else (
            REM mklink /D (directory symlink) のみ使用
            REM Stream モードでは /J (junction) がフォールバックできないため /D 一択
            REM 冒頭で管理者権限チェック済みのため、通常はここで成功する
            mklink /D "!SRC!" "!DST!" >nul 2>&1
            if errorlevel 1 (
                echo     [NG] mklink /D 失敗。データは Google Drive 側に移動済み:
                echo          !DST!
                echo          手動復旧コマンド:
                echo            robocopy "!DST!" "!SRC!" /E /MOVE /NFL /NDL /NJH /NJS /NC /NS /NP
                set /a ERRORS+=1
            ) else (
                echo     [OK] symlink 作成
                set /a MOVED+=1
            )
        )
    )
)

echo.
echo   集計: 移動 !MOVED! / 既 symlink !SKIPPED_LINK! / memory 無 !SKIPPED_NOMEM! / エラー !ERRORS!
echo.

REM ===== 検証 =====
echo --- ローカル側 (symlink/junction 確認) ---
dir /AL "%LOCAL_PROJECTS%" >nul 2>&1
for /d %%P in ("%LOCAL_PROJECTS%\*") do (
    dir "%%P" /AL 2>nul | findstr /i "memory" 2>nul
)
echo.
echo --- Google Drive 側 ---
dir /B "%GDRIVE_TARGET%" 2>nul
echo.

echo ===============================================
if !ERRORS! equ 0 (
    echo   [OK] セットアップ完了
    echo.
    echo   次のステップ:
    echo     1. VSCode を再起動
    echo     2. このプロジェクトを開いて Claude Code が正常起動するか確認
    echo     3. Google Drive for desktop は「ミラーリング」モード推奨
) else (
    echo   [!!] エラーが !ERRORS! 件あります。上記ログを確認してください。
)
echo ===============================================
echo.
pause
endlocal
