@echo off
setlocal EnableExtensions
chcp 65001 >nul

set "BASE_URL=https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/jir"
set "TARGET_DIR=C:\xampp\htdocs\jir"
if not "%~1"=="" set "TARGET_DIR=%~1"

where curl >nul 2>&1
if errorlevel 1 (
    echo [ERROR] curl command not found.
    exit /b 1
)

if not exist "%TARGET_DIR%" mkdir "%TARGET_DIR%"

set "REMOTE_NAME=style.css"
call :fetch "style.css" || exit /b 1

set "REMOTE_NAME=%%E5%%9B%%BD%%E4%%BC%%9A%%E4%%BC%%9A%%E8%%AD%%B0%%E9%%8C%%B2%%E3%%81%%AB%%E3%%81%%BF%%E3%%82%%8B%%E8%%BB%%8D%%E4%%BA%%8B%%E3%%83%%BB%%E3%%82%%A4%%E3%%83%%B3%%E3%%83%%86%%E3%%83%%AA%%E3%%82%%B8%%E3%%82%%A7%%E3%%83%%B3%%E3%%82%%B9%%E3%%81%%AE%%E8%%B3%%AA%%E7%%96%%91%%EF%%BC%%881947%%E3%%80%%9C1984%%E5%%B9%%B4%%EF%%BC%%89.html"
call :fetch "国会会議録にみる軍事・インテリジェンスの質疑（1947〜1984年）.html" || exit /b 1

set "REMOTE_NAME=%%E5%%9B%%BD%%E4%%BC%%9A%%E4%%BC%%9A%%E8%%AD%%B0%%E9%%8C%%B2%%E3%%81%%AB%%E3%%81%%BF%%E3%%82%%8B%%E8%%BB%%8D%%E4%%BA%%8B%%E3%%83%%BB%%E3%%82%%A4%%E3%%83%%B3%%E3%%83%%86%%E3%%83%%AA%%E3%%82%%B8%%E3%%82%%A7%%E3%%83%%B3%%E3%%82%%B9%%E3%%81%%AE%%E8%%B3%%AA%%E7%%96%%91_%%E7%%AC%%AC1%%E9%%83%%A8_1947%%E3%%80%%9C1959%%E5%%B9%%B4.html"
call :fetch "国会会議録にみる軍事・インテリジェンスの質疑_第1部_1947〜1959年.html" || exit /b 1

set "REMOTE_NAME=%%E5%%9B%%BD%%E4%%BC%%9A%%E4%%BC%%9A%%E8%%AD%%B0%%E9%%8C%%B2%%E3%%81%%AB%%E3%%81%%BF%%E3%%82%%8B%%E8%%BB%%8D%%E4%%BA%%8B%%E3%%83%%BB%%E3%%82%%A4%%E3%%83%%B3%%E3%%83%%86%%E3%%83%%AA%%E3%%82%%B8%%E3%%82%%A7%%E3%%83%%B3%%E3%%82%%B9%%E3%%81%%AE%%E8%%B3%%AA%%E7%%96%%91_%%E7%%AC%%AC2%%E9%%83%%A8_1960%%E3%%80%%9C1977%%E5%%B9%%B4.html"
call :fetch "国会会議録にみる軍事・インテリジェンスの質疑_第2部_1960〜1977年.html" || exit /b 1

set "REMOTE_NAME=%%E5%%9B%%BD%%E4%%BC%%9A%%E4%%BC%%9A%%E8%%AD%%B0%%E9%%8C%%B2%%E3%%81%%AB%%E3%%81%%BF%%E3%%82%%8B%%E8%%BB%%8D%%E4%%BA%%8B%%E3%%83%%BB%%E3%%82%%A4%%E3%%83%%B3%%E3%%83%%86%%E3%%83%%AA%%E3%%82%%B8%%E3%%82%%A7%%E3%%83%%B3%%E3%%82%%B9%%E3%%81%%AE%%E8%%B3%%AA%%E7%%96%%91_%%E7%%AC%%AC3%%E9%%83%%A8_1978%%E3%%80%%9C1984%%E5%%B9%%B4%%E3%%81%%A8%%E8%%80%%83%%E5%%AF%%9F.html"
call :fetch "国会会議録にみる軍事・インテリジェンスの質疑_第3部_1978〜1984年と考察.html" || exit /b 1

echo [OK] Sync completed: %TARGET_DIR%
endlocal
exit /b 0

:fetch
set "LOCAL_NAME=%~1"
set "DST=%TARGET_DIR%\%LOCAL_NAME%"
set "TMP=%DST%.tmp"

echo Downloading %LOCAL_NAME%...
curl -fsSL -o "%TMP%" "%BASE_URL%/%REMOTE_NAME%"
if errorlevel 1 (
    echo [ERROR] Failed to download: %LOCAL_NAME%
    if exist "%TMP%" del /f /q "%TMP%"
    exit /b 1
)
move /y "%TMP%" "%DST%" >nul
exit /b 0
