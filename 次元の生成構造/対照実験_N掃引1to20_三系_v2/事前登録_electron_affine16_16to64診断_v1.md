# 事前登録 — electron affine16→64 不安定性診断 v1

状態: **コード・smoke のみ許可。本走行は未実行**  
対象: `N=12`, electron seed `δ=0.1`, `T=42000`

## 1. 目的と主張境界

既存 electron 強シードは、後半に

- 偶数帯: `η={0,2,4,6}`
- 奇数帯: `η={1,3,5,7}`

という parity-compatible 64 セル支持を持つ。一方、pump `(2,0)` と electron seed
`(1,3)` を通る理想アフィン軌道は

$$
\eta_*(k)=-3(k-2)\pmod 8
$$

の16セルである。本診断は、16セル軌道外の微小成分が成長して64セル支持へ移るのか、
どの sector が不安定なのかを判別する。

本診断は F v1 の証明や新粒子同定ではない。特に projection 腕は保存力学そのものではなく、
各 full-step 後に外部射影を入れる**診断介入**である。

## 2. 変更しない条件

| 項目 | 固定値 |
|---|---|
| N / M | 12 / 66 |
| レジスタ | Nn=16 / Neta=8 |
| 時間 | T=42000、後半窓 `[21000,42000)` |
| 親 | seed=2 |
| D | cell=(2,0)、order=6 |
| pump | `(k,η)=(2,0)` |
| electron seed | `(k,η)=(1,3)`, δ=0.1 |
| 力学 | `unified_interaction_v1.py` の F v1（無改変 import） |
| D/G/S | `unified_dimension_v1.py` / `unified_readout_v3.py` / `selection_v1.py` |
| 保存 | sector power は毎step、raw C は事前固定snapshot |

既存 `run_nsweep_three_series_v2.py`、F/D/G/S、親生成器は編集しない。専用runnerが
初期条件または full-step 後の診断介入だけを外側から与える。

## 3. sector の機械定義

128セルを排他的に次の3 sectorへ分ける。

1. `ideal16`: `η=η_*(k)` の16セル。
2. `allowed_off`: `k mod 2 = η mod 2` かつ `ideal16` 外の48セル。
3. `forbidden`: `k mod 2 != η mod 2` の64セル。

従って `ideal16 + allowed_off + forbidden = total` を毎stepで検算する。

## 4. 腕と唯一の差

| 腕 | 条件 | 本数 |
|---|---|---:|
| A | 標準 electron δ=0.1。介入なし | 1 |
| B | 各 F v1 full-step 後、`ideal16` 外を0化し、step前の総normへ全状態を再正規化 | 1 |
| C | parity-allowed off-orbit `(1,1)` に初期振幅 ε を注入 | 3 |
| D | parity-forbidden `(1,0)` に同じ ε を注入 | 3 |

固定する注入振幅は

$$
\epsilon\in\{10^{-15},10^{-10},10^{-5}\}.
$$

C/D は注入後の**全初期状態**を標準Aの総normへ再正規化する。`PF`, `PB`, pump power,
primary seed power, injected-cell power, `Ptotal` を再正規化後の実測値として保存する。
注入方向はprimary seedと同じ正規化親Fourierモード
`Csec[:,1] / ||Csec[:,1]||` とし、複素位相も同じに固定する。
ε最大でも注入powerは `10^-10` であり、δ=0.1 seed powerに比べ十分小さい。

B は各stepについて、F v1 step直後・射影前のsector power、除去power、射影後のsector
power、再正規化倍率、step前normとの差を保存する。F v1 の内部stepは変更しない。
既存 `rec_m_*` は無改変 `RecordingEngine` がF step直後に記録するためBでは射影前、
`m_*` のD/G/S readoutは射影後となる。この時点差はJSONにも明記する。

## 5. raw state と閾値

raw `C(M,Nn,Neta)` は初期step 0、step `1,2,10,100,1000`、および500stepごとに保存する。
B は同じstepで射影前と射影後を分けて保存する。sector powerは全42000 stepを保存する。

sector fraction の初回到達を次の固定閾値で読む。

$$
10^{-30},10^{-24},10^{-20},10^{-16},10^{-12},10^{-10},10^{-8},10^{-6},10^{-4},10^{-2}.
$$

各閾値について `allowed_off`, `forbidden`, 両者の和を別々に記録する。数値的非零と
物理的支持を混同しないため、閾値を事後選択しない。

## 6. 事前予測と判別

- Aが16セルの厳密不変軌道なら `allowed_off=forbidden=0` を保つ。64セル化は起きない。
- Aで `allowed_off` が先に増え、Cの到達時刻が `log ε` に沿って前進するなら、
  parity-allowed 横方向の不安定性を支持する。
- DがCより抑制されるなら、64セル支持は任意roundoffではなく parity sector 選択を持つ。
- Bで大域的な0.7遷移や `4/7,8/15` 近傍への再配分が消えるなら、軌道外自由度がその遷移に
  必要である。ただしBは射影・再正規化を伴うため、代替の自然力学とは呼ばない。
- Bでも主要readoutがAと同じなら、64セル支持はその大域readoutに必要でない可能性がある。

## 7. 実行ゲート

1. smoke でmask分割 `16+48+64=128`、初期norm、注入sector、B射影後off=0を確認する。
2. 本走行Aを最初に実行し、既存
   `nsweep_electron_T42000_d0.1_N12_v2.npz` の全99配列について、key、dtype、shape、
   NaN bitを含む連続byte列が一致しなければ停止する。
3. A合格前にB/C/Dを実行しない。
4. 全腕を逐次実行し、source/reference SHA-256を走行前後で照合する。
5. 既存・部分成果物は上書きしない。

## 8. 停止判定

- standard byte比較、norm、sector partition、内部真空のどれかが不一致なら即停止。
- B射影後にoff powerが厳密0でない場合は即停止。
- C/Dで初期注入セルまたはsectorが宣言と違う場合は即停止。
- 上記機械検査に合格した後にのみ、閾値時刻とreadout差を解釈する。
