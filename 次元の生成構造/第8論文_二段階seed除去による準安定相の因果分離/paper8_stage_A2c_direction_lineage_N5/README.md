# Stage A2c: N=5 方向基底の系譜追跡

Stage A0でbitwise再現された第7論文 `N=5` / `float64` / 明示seedあり軌道を `step=0..5000` だけ再実行し、既存コード内部の方向基底・射影行列・既存横摂動方向を保存して比較する。

物理軌道、方向の物理的定義、q5〜q8、flux、eta_noise、高精度、Delta/seed振幅掃引、`N=40`、`N=300`は追加しない。既存原本とStage A0成果物は読み取り専用である。

## 実行順序

```bash
python3 verify_sources.py
python3 replay_and_extract_bases.py
python3 reconstruct_transverse_directions.py
python3 analyze_direction_lineage.py
python3 make_figures.py
python3 make_report.py
```

各工程は前工程の成功記録を確認する。不一致時は自動修正せず停止する。

## 保存する正本

- `B0`: 初期親平面
- `Bdom(t)`: 瞬時支配平面
- `D34(t) = s4_new_dirs(B0, Bdom(t))`
- `P34(t) = D34(t)D34(t)^T`: D34比較の正本
- `S4(t) = s4_basis(...)` と `PS4(t)`
- `Tperp`: 第7論文横摂動コードの複素 `eta`。比較用実2次元部分空間は `span(eta.real, eta.imag)` とする。

基底列の符号・交換・内部回転と、2次元部分空間自体のambient空間内回転を別々に記録する。

## 停止範囲

報告書完成後に停止し、高精度、Delta掃引、`N=40`、`N=300`、Stage B/Cへ進まない。
