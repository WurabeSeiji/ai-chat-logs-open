# グノモン正写像 計算機検証プログラム

## 概要

本フォルダには、グノモン正写像による誘導計量の正しさを**計算機代数システム（SymPy）**で独立検証するPythonプログラムが格納されています。

## 必要環境

- Python 3.8 以上
- SymPy (`pip install sympy`)

## ファイル一覧

| ファイル | 内容 |
|---------|------|
| `verify_euclidean.py` | ユークリッド版の検証（n=2,3,4,5次元） |
| `verify_lorentzian.py` | Lorentzian版の検証（4次元、署名 +,+,+,-） |

## 実行方法

```bash
# ユークリッド版（P-1〜P-5）
python3 verify_euclidean.py

# Lorentzian版（P-8）
python3 verify_lorentzian.py
```

## 検証項目

| 項目 | 内容 | スクリプト |
|-----|------|-----------|
| P-1 | g_μν · g^μν = I_n | 両方 |
| P-2 | Christoffel 記号 Γ^α_μν | 両方 |
| P-3 | R_μν = (n-1)/R² g_μν, R = n(n-1)/R² | 両方 |
| P-4 | G_μν + Λg_μν = 0 (Λ = (n-1)(n-2)/(2R²)) | 両方 |
| P-5 | n=2,3,4,5 での一般検証 | verify_euclidean.py |
| P-8 | Lorentzian署名での検証 | verify_lorentzian.py |

## 実行結果（2026-04-06）

全項目 ✅ ALL PASS
