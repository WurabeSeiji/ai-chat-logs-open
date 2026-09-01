# ChatGPT自身への厳密対照実験指図書

作成日: 2026-09-01

## 目的
正本 `pass2_run.py` を忠実に複製し、N=5, L=124 のみを 500 step 実行して、既存のChatGPTテスト系 N=5, L=124, 500 step と比較する。

## 絶対条件
1. 正本 `pass2_run.py` と正本 `original_engine.py` をGoogle Driveから直接取得する。
2. 正本 `data/hm_N5/parent_v.npz` をそのまま使い、初期値を再生成しない。
3. 正本コピーに対する力学条件の変更は `STEPS=40000` を `STEPS=500` にする一箇所だけとする。N=5, L=124, `np.linalg.eigh`, `complex128/float64`, H生成式、更新式、metricsを変更しない。
4. 差分を `diff -u` で保存する。
5. 正本は `key_steps` に 750以上が含まれるため、500 step化すると走行終了後のkey_steps出力でIndexErrorになり得る。その場合、力学走行CSVとstatesが501状態生成済みであることを確認し、力学コードには追加修正を加えない。後処理は別スクリプトで行う。
6. 生成された0..500の生データを、既存ChatGPTテスト系のN=5,L=124データと全step比較する。
7. さらに、当時保存された正本40000-stepデータの先頭501行とも比較する。
8. プログラム、入力、生成データ、比較CSV、図、環境情報、SHA256、分析をすべてGoogle Driveへ保存する。
9. 分析には保存先Google Driveパスを明記する。

## 判定
- exact-copy vs ChatGPT test の H_total が全501点で一致するか。
- H_perp/H_total の差の最大値を計算する。
- 当時保存正本との最初の不一致stepを記録する。
