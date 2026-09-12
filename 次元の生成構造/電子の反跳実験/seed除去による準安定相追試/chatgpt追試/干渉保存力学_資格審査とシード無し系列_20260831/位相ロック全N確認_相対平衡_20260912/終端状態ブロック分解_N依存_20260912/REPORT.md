# 追試B：ロック終端状態の最小ゼロ閉塞ブロック分解の N 依存（論文9 の提案検証）

日付: 2026-09-12。読出しのみ・物理無変更・新規走行なし（既存 10000step 終端状態、N=3..40）。
動機: 論文9（検討）が「終状態のブロック分解の N 依存＝種の目録が $S_N$ 不変量に立つか」を
既存データで検定可能（未実施）とした。既存 10000step 終端状態で実施し、立つ／立たないを確定させる。

## 方法

各 $N$ の終端状態 $z=S[-1]$（ロック済み）について、全 unordered pair の二波自己閉塞残差
$\varepsilon_{ij}=|z_i^2+z_j^2|/(|z_i|^2+|z_j|^2)$ を総当たり監査（exact $\varepsilon<10^{-9}$／near $<10^{-6}$／最小 $\varepsilon$）。
総閉塞 $\Sigma z^2$ と等モジュラー偏差も記録。対照として同 $N$ の床（step0＝make_parent 床）でも同監査。

## 結果（`results/terminal_block_summary.csv`）

- **終端状態は全 $N=3$–$40$ で厳密等モジュラー**（振幅偏差 $\sim10^{-13}$）、$\Sigma z^2\approx0$ 保存（$\sim10^{-13}$）。
- **終端状態に exact $\pm i$ 対はゼロ**（全 $N$）。最小 $\varepsilon$ は $N=3$ の $0.5$ から大 $N$ の $\sim2\times10^{-6}$ まで漸減するのみ（near は $N=7$ に2本のみ）。
- **make_parent 床（step0）も exact $\pm i$ 対ゼロ**（最小 $\varepsilon$ $0.33$→$\sim2\times10^{-7}$）——論文4・実験15（make_parent は strict 光子型 pair 0、near は高 $N$ のみ）と整合。

## 判定（探究型物理学者ロール）

- 分類: 数値確認（提案検定の実施）。判定: 論文9 の提案は**データ系列の終端では立たない**（否定的結論を確定）。
- **最小ブロック（$\pm i$ 対）分解は巡回Fourier／整数K 高対称床（$N=8k$、論文4）に固有の構造**であり、実験データ系列（make_parent 床→ロック終端）には exact には存在しない。終端状態は厳密等モジュラー化するが、その等モジュラー環は $\pm i$ 対や三波へ exact 分解されない（generic 位相）。
- したがって「終端ブロック分解の N 依存＝種の目録」は、**終端状態の性質としては成立せず**、ブロック構造（もしあれば）は高対称床側（論文4）に宿る。これは論文9 の本旨（情報の担い手はブロックで、無名性と両立する型は床の閉塞構造にある）と整合し、かつ「終端では消える」ことを明確化する。

## 未決 / 留保

- near-pair（$\varepsilon\sim10^{-6}$）が大 $N$ で漸近することの意味（連続極限で exact 化するか）は本追試の範囲外。
- 3波以上のブロックの exact 監査は pair が exact 0 のため省略（pair が無い状態で高次 exact ブロックは期待できない）。

## 成果物

- `analyze_terminal_block_decomposition_20260912.py`（正本 adjacency を SHA 照合の上 import・読出しのみ）
- `results/terminal_block_summary.csv`, `results/実行ログ_20260912.log`
- `README.md`, `SHA256SUMS.txt`, `run_all.sh`
- 入力: 既存 `../../N3_N40_long10000_20260905/results/hm_N{N}_den_{N}_states_10000.npz`、床対照 `../../make_parent型初期値_自己無撞着構造_20260907/full_N3_N40_sweep/states/`（読み取りのみ）。依存: numpy。
