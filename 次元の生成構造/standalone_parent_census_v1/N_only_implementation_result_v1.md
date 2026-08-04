# Nだけを理論入力とするmake_parent 実装結果 v1

## 実装

- 生成器: `make_parent_n_only_v1.py`
- 独立一覧器: `parent_wave_table_n_only_v1.py`
- 契約: `N_only_parent_contract_v1.md`
- テスト: `test_make_parent_n_only_v1.py`

公開関数は次の一つである。

```python
make_parent(N, seed=None, max_retries=3)
```

seed明示時は1回だけ実行する。seed未指定時はOS seedによる初回と最大3回の
リトライを行う。未収束値を成功として返さない。

## N=5 採用seedによる生成結果

seedを1から順に1回ずつ実行し、全5波が最初に収束したseed=11を採用した。
採用後の本計算は `make_parent(5, seed=11)` の1回である。

| 項目 | 結果 |
|---|---:|
| N | 5 |
| M=N(N−1)/2 | 10 |
| 波行列shape | 10×5 |
| 安定した周波数波 | 5 |
| 保存した関係・波成分 | 50 |
| 最大自己無撞着残差 | 5.453e−13 |
| 最大波別零閉塞絶対残差 | 2.902e−16 |
| 全波合成零閉塞絶対残差 | 1.471e−16 |
| 同seed再実行 | 配列bit一致 |

生成データは `parent_N5_seed11_v1/`、日本語一覧は
`wave_table_N5_seed11_v1/wave_table.md` にある。seedの選択記録は
`seed_selection_N5_v1.md` にある。

seed=11からPCG64で白色雑音近似の疑似乱数を生成しており、位相・振幅を事後に
手作業で変更していない。同じNとseedから入力と結果をbit単位で再現できる。

## N=5 OS seed本実行

`seed=None, max_retries=3`で一度実行した。初回＋3リトライの4候補はすべて
自己無撞着許容値へ収束せず、仕様どおりアボートした。成功結果への差し替えや
追加リトライは行っていない。

全4 seed、各残差、失敗周波数は
`actual_parent_N5_os_v2/failure_manifest.json` に保存した。この失敗履歴は数値実装の
監査情報であり、粒子生成率その他の物理量として解釈しない。
