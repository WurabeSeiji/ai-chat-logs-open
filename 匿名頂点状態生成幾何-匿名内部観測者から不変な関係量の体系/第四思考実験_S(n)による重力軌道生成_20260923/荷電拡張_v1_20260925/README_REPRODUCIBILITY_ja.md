# Paper 4 荷電拡張 v1 — 再現手順

## 目的

第四思考実験の状態依存軌道生成系へ、離散電荷 `(nA,nB)=(+1,-1)` と無次元結合 `Gamma` を最小追加し、

- クーロン優勢二体軌道
- GW 読み出し
- EM 双極子波読み出し
- GW/EM 放射反作用
- 5つの永続状態の加法（log側）表現
- 5つの永続状態の指数・乗法表現

を同一データから再現する。

## 基準ケース

- mass ratio: `mA/mB = 1`
- charge integers: `(nA,nB)=(+1,-1)`
- `|F_C/F_G| = 100`
- `Gamma = 25`
- `alpha = 101`
- orbit dilation: `10`
- angular mode: `newtonian`
- GW radiation reaction: ON
- EM radiation reaction: ON
- radial orbits: `3`
- steps/orbit: `4000`

この条件は `results_rho100_d10_v2/charged_v1_qm1_n1_-1_rho100_summary.json` に保存されている。

## 1. 数値実験の再実行

```bash
python run_paper4_charged_extension_v1.py \
  --mass-ratio 1 \
  --nA 1 \
  --nB -1 \
  --coulomb-to-gravity 100 \
  --scale-orbit-with-alpha \
  --orbit-dilation 10 \
  --angular-mode newtonian \
  --steps-per-radial-orbit 4000 \
  --radial-orbits 3 \
  --out results_rho100_d10_v2
```

出力:

- `results_rho100_d10_v2/charged_v1_qm1_n1_-1_rho100.csv`
- `results_rho100_d10_v2/charged_v1_qm1_n1_-1_rho100_summary.json`

## 2. 放射式の独立照合

```bash
python validate_charged_extension_v1.py
```

固定軌道に対する数値平均放射と解析式を比較する。

## 3. 図の再生成

```bash
python plot_paper4_charged_extension_v1.py
```

デフォルトでは上記基準ケースを読み、`figures_rho100_d10_v2/` に図を作る。

明示指定する場合:

```bash
python plot_paper4_charged_extension_v1.py \
  --input-csv results_rho100_d10_v2/charged_v1_qm1_n1_-1_rho100.csv \
  --summary-json results_rho100_d10_v2/charged_v1_qm1_n1_-1_rho100_summary.json \
  --output-dir figures_rho100_d10_v2
```

## 状態表現

### 加法 / log側

```text
(chi, p, e, phi, t)
```

### 指数 / 乗法

```text
U = exp(i chi)
P = p
E = e
H = exp(i phi)
Q = exp(t/tau0)
```

荷電基準ケースの時計スケールは

```text
tau0 = 2*pi*sqrt(a0^3/alpha)
```

を使う。

## 図一覧

- `figure01_charged_orbit.*` — 荷電二体軌道
- `figure02_log_five_states_panel.*` — 加法/log側5状態
- `figure02a-e_*` — 加法/log側各状態の個別図
- `figure03_exponential_five_states_panel.*` — 指数/乗法5状態
- `figure03a-e_*` — 指数/乗法各状態の個別図
- `figure04_first_charged_experiment_overview.*` — 第一荷電数値実験の全体図
- `figure05_baseline_vs_charged_normalized_orbit.*` — 無荷電基準との規格化軌道比較
- `charged_v1_derived_five_states.csv` — U,H,Q を追加した再現用状態表

## 注意

この v1 は Einstein–Maxwell 系の自己無撞着な厳密解ではない。Newtonian Coulomb 中心力、Paper 4 の状態/読み出し構造、leading GW/EM multipole radiation reaction を組み合わせた診断用数値実験である。
