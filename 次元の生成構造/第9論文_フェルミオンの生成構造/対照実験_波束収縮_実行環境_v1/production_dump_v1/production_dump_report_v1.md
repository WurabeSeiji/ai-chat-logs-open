# 本番ダンプ報告 v1（倍音別複素係数の系統取得）

実施日: 2026-07-28
ランナー: `20260715/run_system_A_localization_exchange_R_sweep_instrumented_v1.py`
駆動: `run_production_dump_v1.py` ＋ B2全スペクトル追補（下記）
総容量: 608MB（npz 185ファイル）

## 取得一覧

| ラン | ケース | R点数 | M（記録倍音範囲） | npz | 捕捉率最小 |
|---|---|---|---|---|---|
| 01_femtofocus_R137_B12 | custom A=1,B=1,2（固定正規化・原本と同一121点） | 121 | 6 | 121 | ≈1（15桁） |
| 02_B12_keyR | 同上 | 代表9点 | 6 | 9 | ≈1 |
| 03_oddN_B1_keyR | odd_kernel 1:1 | 9 | 3 | 9 | ≈1 |
| 03_oddN_B2_keyR | odd_kernel 1:2 | 9 | 4 | 9 | **0.9770**（下記参照） |
| 03_oddN_B2_keyR_fullM | odd_kernel 1:2（追補） | 9 | **255（全スペクトル）** | 9 | 0.99995 |
| 03_oddN_B3_keyR | odd_kernel 1:3 | 9 | 5 | 9 | ≈1 |
| 03_oddN_B5_keyR | odd_kernel 1:5 | 9 | 7 | 9 | ≈1 |
| 03_oddN_B15_keyR | odd_kernel 1:15 | 9 | 17 | 9 | ≈1 |
| 03_oddN_B63_keyR | odd_kernel 1:63 | 9 | 65 | 9 | ≈1 |

代表9点の R: 0.0, 0.5, 1.0（対照）／R128厳密値 0.686671465671125／R137厳密値
0.6971778791282474／femto窓中心 0.697177927／off-resonance 0.68, 0.70, 0.71。
全ラン stride=1（全257衝突を記録）。

## 検証

1. **内蔵parity**: 01（femtofocus再現）の既存4CSVは、係数ダンプを行いながら
   **原本とバイト単位で同一**。計測が物理へ与える影響ゼロを本番規模で再確認。
2. **捕捉率監査**: 各npzの coverage 配列に全記録時点の捕捉パワー率を保存。

## 発見事項: B=2（偶数nh）カーネルの広帯域床

`odd_harmonic_kernel(u, nh)` = sin((nh+1)u)/((nh+1)sin u) は nh が偶数のとき
sin(u)=0 の格子点で符号不連続（±1 分岐処理）となり、B=2 状態は
主線 {−3, +1, −1}（各32.56%、キャリア−1シフト込み）に加えて
**全周波数に広がる一様な広帯域床（各モード≈4.6e-5、合計≈2.3%）**を持つ。
奇数 nh（1,3,5,15,63）のカーネルは滑らかで床を持たない。

対応: B2 は M=255（512モード中511、欠けは Nyquist n=−256 の1モードのみ≈5e-5）で
全スペクトル追補を取得（03_oddN_B2_keyR_fullM）。M=4 版も対照として保持。
偶数倍音対照の解析では fullM 版を正とすること。

## npz の読み方

```python
z = np.load("harmonic_coeffs_..._v1.npz")
coeffs = z["coeffs"]        # complex128 (n_records, 2, n_harm, 16)  ch軸=[A,B]
harms  = z["harmonics"]     # 符号付き倍音番号
colls  = z["collisions"]    # 衝突番号
meta   = json.loads(str(z["meta"]))   # R, T, t, r, パケット仕様, FFT規約
# 交差相関行列（eta縮約）: C[m,n] = Σ_η coeffs[k,ch,m,η] · conj(coeffs[k,ch,n,η])
C = coeffs[k, ch] @ coeffs[k, ch].conj().T
```

これで調査指示の実行順序⑤（倍音別複素係数の抽出）が完了。
次段階は⑥交差相関行列の再構成と、⑦局在性との時系列比較。
