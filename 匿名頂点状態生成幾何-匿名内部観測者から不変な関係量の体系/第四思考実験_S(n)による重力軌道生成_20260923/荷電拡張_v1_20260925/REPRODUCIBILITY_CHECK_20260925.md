# Paper 4 荷電拡張 v1 — 再現性確認 2026-09-25

## 対象

基準ケース:

- mass ratio `mA/mB = 1`
- charge integers `(nA,nB)=(+1,-1)`
- `|F_C/F_G| = 100`
- `Gamma = 25`
- `alpha = 101`
- orbit dilation `10`
- angular mode `newtonian`
- GW RR ON / EM RR ON
- 3 radial orbits
- 4000 steps/orbit

## 図の再生成確認

次を実行した。

```bash
cd paper4_charged_extension_v1
python plot_paper4_charged_extension_v1.py --output-dir figures_repro_check
```

正常終了し、軌道図、log側5状態、指数側5状態、第一荷電実験 overview、無荷電基準との規格化軌道比較を再生成した。

出力された時計スケール:

```text
tau0 = 3297487.24417781
```

## Google Drive 保存図とのバイト同一性確認

代表として Figure 01 と Figure 05 を Google Drive から再取得し、再生成PNGと SHA-256 を比較した。

### Figure 01

```text
ee6f784bb830b7d14e19513eca1e1b7e87c0571bc76ffbea1b8c5b995b1cf239
```

Google Drive 保存版と再生成版が一致した。

### Figure 05

```text
f29d58275c79c42de62b01be76ef601ecce5238d6a95b87b6dbd78412bd6c214
```

Google Drive 保存版と再生成版が一致した。

## 保存構成

Google Drive の `荷電拡張_v1_20260925` 直下に、実行プログラム3本、基準結果、無荷電基準、図、解析文書、ZIPを保存した。

また `raw_runs/` に以下の診断系列を保存した。

- original baseline
- Newtonian baseline
- rho=100 initial run
- rho=100, orbit dilation=10
- rho=100, EM RR off
- rho sweep: 1, 10, 100, 1000
- convergence: 2000, 4000, 8000 steps/orbit

## 判定

少なくとも保存済み基準CSV/JSONから、保存済み plotting program によって Figure 01–05 を再生成できることを確認した。代表図 Figure 01/05 は Google Drive 保存版と SHA-256 が一致した。
