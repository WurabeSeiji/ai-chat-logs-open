# Nのみを使う白色零閉塞make_parentと独立粒子表

作成日: 2026-08-04<br>
現在の正本: 生成・読出し実装 v2、粒子表レイアウト v3

## 1. 正本プログラム

- `make_parent_white_harmonics_n_only_v2.py`
- `particle_table_white_harmonics_n_only_v2.py`
- `test_white_harmonics_n_only_v2.py`

旧 `make_parent_n_only_v1.py`、`make_parent_white_v1.py` と旧一覧器群は正本ではない。
旧成果ディレクトリの `WITHDRAWN.md` / `SUPERSEDED.md` を優先して読む。

## 2. 入力契約

```python
make_parent(N, seed=None, max_retries=3)
```

- `N`: 唯一の理論入力。体数と標本分解能を同じNに固定する。
- `seed`: 疑似乱数系列の再現指定であり、理論パラメータではない。
- `max_retries`: seedを指定しない場合の失敗処理回数であり、理論パラメータではない。

Nからだけ

$$
M=\frac{N(N-1)}{2},\qquad \lambda_0=\frac{2\pi}{N}
$$

を導く。独立な倍音段数H、別の時間分解能、倍音次数、倍音振幅、倍音位相を入力しない。

## 3. seed規約

- seedを指定した場合: そのseedで1回だけ計算し、未収束ならアボートする。
- seedを指定しない場合: OS乱数から128 bit seedを得る。初回と最大`max_retries`回の再試行後も未収束ならアボートする。
- 実験用の逐次探索: `--search-seed-from 1` により1,2,3,...の順で試し、最初に収束したseedを採用する。失敗を含む全seedと残差をmanifestへ保存する。

## 4. 白色雑音からの二重零閉塞構成

既存make_parentと同じ自己無撞着方程式を解き、親ベクトル

$$
v\in\mathbb C^M,\qquad \sum_{m=1}^{M}v_m^2=0
$$

を得る。次に各関係波についてN×2実Gaussian疑似白色雑音を生成し、QR分解から
等ノルム直交実ベクトル $q_{m,1},q_{m,2}$ を得て

$$
w_m=\frac{v_m}{\sqrt{2}}(q_{m,1}+i q_{m,2})
$$

とする。したがって全M本で

$$
\sum_{j=1}^{N}w_{m,j}^2=0
$$

が成立し、さらに

$$
\sum_{m=1}^{M}\sum_{j=1}^{N}w_{m,j}^2=0
$$

も成立する。生成器は倍音を置かない。

## 5. 倍音の事後読出し

独立一覧器は生成器をimportせず、保存された $W\in\mathbb C^{M\times N}$ だけを読む。
各関係波にN点DFTを行い、次を分離する。

- 直流: $k=0$
- 基本波: $|k|=1$
- 倍音: $|k|\ge2$
- 正周回枝: $k>0$
- 負周回枝: $k<0$

主表は1関係波につき1行である。N=5なら10行、N=40なら780行になる。
全DFT成分の次数・振幅・強度・位相・住所は `harmonics_ja.csv` に一成分一行で保存する。
位相の異なる複素波形を同一視しない。

## 6. 粒子構造の候補読出し量

表示を削らず、次の内部量を全関係波について計算する。

- B/F/E位相網: B=全DFT成分が同相、F=二成分が逆相、E=中間位相網。外部の位相分解能や144セルは使わない。
- 倍音偶奇のB/F型: 奇数倍音強度と偶数倍音強度、および純奇数・純偶数・混合の分類。
- 粒子／反粒子型: 正負周回強度の偏極符号。
- 電荷型: 支配住所 $m/n$、$\sin^2(\pi m/n)$、周回符号を付けた住所電荷量。
- 質量型: 正負倍音枝の $\det\Gamma=N_+N_- - |\langle c_+|c_-\rangle|^2$ と平方根。
- スピン型: 半周期で $+w$ に戻る1:1被覆、$-w$になり二次量が戻る2:1被覆、状態と二次量の回帰次数。
- 寿命型: 循環自己相関が最初に $1/e$ 未満になる標本遅れ。
- 零閉塞: 各関係波の絶対残差・相対残差と全体系残差。

これらは「物理量として一致するものだけを残す」ための表ではない。物理量へ対応する
可能性のある内部量を削らずに比較する探索表である。計算式と数値は保存し、対応の身分は
列名の「型」で区別する。

## 7. 本実行

```bash
python3 make_parent_white_harmonics_n_only_v2.py 5 \
  --search-seed-from 1 --output parent_white_harmonics_N5_v2

python3 make_parent_white_harmonics_n_only_v2.py 40 \
  --search-seed-from 1 --output parent_white_harmonics_N40_v2

python3 particle_table_white_harmonics_n_only_v2.py \
  --input parent_white_harmonics_N5_v2 \
  --input parent_white_harmonics_N40_v2 \
  --output particle_tables_white_harmonics_N5_N40_v3
```

採用seedはN=5が2（seed 1は未収束）、N=40が1。逐次探索記録は各manifestにある。

## 8. 出力

- `parent_white_harmonics_N5_v2/`
- `parent_white_harmonics_N40_v2/`
- `particle_tables_white_harmonics_N5_N40_v3/summary.md`
- `particle_tables_white_harmonics_N5_N40_v3/N5/particle_table.md`
- `particle_tables_white_harmonics_N5_N40_v3/N40/particle_table.md`
- 各Nの `particle_waves_ja.csv`, `harmonics_ja.csv`, `census.json`

`summary.md` と各Nの `particle_table.md` は、波ごとに B/F/E、電荷型量、
質量型量、スピン型量を同一行で比較できる粒子主表から始まる。
