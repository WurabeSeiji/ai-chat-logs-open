# 第3章 複素平面読出し図 — step0・終了時・凝縮中心拡大の3グリッド

（N=3..40 段1+2+3 スイープ再現論文・第3章。総括論文と Concept DOI
10.5281/zenodo.22317635 を共有する。Version DOI 10.5281/zenodo.22317636。
式番号は第1章・第2章（式1〜式25）からの連番）

## 1. 目的

第2章のスイープが保存した状態 npz（Δτ=2π/N 走行）から、各辺の複素波 z_e ∈ ℂ を
複素平面に読み出した3枚のグリッド図——(1) step0、(2) 終了時（step500）、
(3) 終了時の最大角クラスター拡大——を生成する処理を記述する。3図の設計意図
（インフレーション図の始点・終点・終点内部を配置側から裏づける役割分担、および
「図の複素平面 ≠ 親平面 Π」の注意と橋渡しの恒等式）は**第2章 §2.6（式25）**で
与えられており、本章はその実装・実行・観察を固定する。データは読み出しのみで、
力学・保存データには一切触れない。

## 2. 理論的背景

状態と測定の数学は第2章（式22〜式25）に依る。本章の図固有の規約を式26〜式28 として
定義する。

**(式26) 重複計数（縮退の表記）** — 描画対象の複素値集合 {z_e} を、実部・虚部の
**小数第12位丸め**で同値類に分け、要素数 c > 1 の類の位置に「xc」を注記する:

    class(w) = ( round(Re w, 12), round(Im w, 12) )

12桁丸めは倍精度の実効桁（〜16桁）より粗く、丸め誤差程度の差を同一視しつつ
物理的に有意な分離（本系列では相対 10⁻⁴ 以上）は別扱いにする閾値である。
拡大図では丸めを第15位（式28）に締め、機械精度レベルの一致のみを「xn」と数える。

**(式27) パネルスケール規約** — 各パネルの軸範囲は、そのパネルの全体振幅
r = max_e |z_e| を用いて [−1.15r, +1.15r] の正方（equal aspect）とし、
**目盛数値は実値のまま**表示する（正規化した値に貼り替えない）。原点から各点への
線分は波の偏角と振幅を同時に読むための補助線である。

**(式28) 最大角クラスター抽出（拡大図の算法）** — 終了時状態を全体振幅で規格化した
座標で**粗解像度 1/100** に丸めて群化し、最大要素数の群を拡大対象とする:

    key(w) = ( round(Re w / amp, 2), round(Im w / amp, 2) ),   amp = max_e |z_e|
    C* = argmax_{key} |{ w : key(w) = key }|
    中心 c = mean(C*),  割れ幅 spread = max_{w∈C*} |w − c|,  窓幅 = 1.4 × spread

パネル題に本数 |C*|・spread・amp を印字する。spread は「クラスターが厳密な一点
（厳密縮退）か有限幅の束か」を定量化する（第2章 §2.6 の表・拡大図の行）。

## 3. 実装方法

- 図化プログラム `plot_complex_plane_N3_N40_stage123_v1.py`（本パッケージ同梱、
  読み出しのみ）。描画様式は本シリーズの既存図化プログラム
  （`complex_plane_readout_step0_step2000_20260904/plot_complex_plane_step0_step2000.py`
  のグリッド様式、および `自発的分裂予備実験_v1/N40_state_readout_20260904/
  plot_complex_plane_N40_v1.py` のクラスター拡大算法）を N=3..40 の 8×5 グリッドに
  拡張したもので、算法自体は同一である（様式の系譜を変えない）。
- 入力は第2章の成果物 `results/hm_N{N}_den_{N}_states_500.npz` のみ。**den=N
  （Δτ=2π/N）の走行を各 N の代表として用いる**（step0 は分母に依らず同一なので
  代表選択は終了時にのみ効く。他の分母の終了時配置は保存済み npz から同様に
  読み出せる）。
- 出力は PNG 3枚。データの書き換え・再生成はない。

## 4. 詳細設計

### 4.1 全体フロー

```
[初期化]   入力フォルダの決定（results/）
[ループ処理] draw_grid(step=0)  : N=3..40 の各パネルに step0 状態を描画（式26・式27）
           draw_grid(step=500): 同、終了時状態
           拡大グリッド      : N=3..40 の各パネルに最大角クラスター（式28）
[終了処理] 各図を savefig、未使用2パネルを off
```

### 4.2 全体データフロー

- **入力**: `results/hm_N{N}_den_{N}_states_500.npz` × 38（第2章成果物・読み出しのみ）。
  使用フィールド: `Z`（501×M）、`denominator`、`steps`（整合性 assert 用）
- **パラメータ**（すべてプログラム内定数）:
  - 丸め桁: 12（式26）／15（拡大図内、式28）／粗解像度 2桁（式28の 1/100）
  - スケール係数 1.15（式27）、窓係数 1.4（式28）
  - グリッド 8×5（38使用・2 off）、dpi=180
- **出力**:
  - `fig_complex_plane_step0_N3_N40_stage123.png`
  - `fig_complex_plane_final_N3_N40_stage123.png`
  - `fig_complex_plane_final_zoom_N3_N40_stage123.png`

### 4.3 個別処理

#### 4.3.1 初期化（読み込みと整合性検査）

```python
    18	def load(N, step):
    19	    d = np.load(os.path.join(IN, f'hm_N{N}_den_{N}_states_500.npz'))
    20	    assert int(d['denominator']) == N and int(d['steps']) == 500
    21	    return np.asarray(d['Z'][step], dtype=np.complex128)
```

- 20行の assert は「意図したファイル（den=N・500 step）を読んでいる」ことの検査で、
  不一致なら AssertionError で停止する（描画は行われない）。

#### 4.3.2 ループ処理（1）: グリッド描画（式26・式27）

```python
    26	    for k, N in enumerate(range(3, 41)):
    27	        ax = axs[k]
    28	        z = load(N, step)
    29	        M = N * (N - 1) // 2
    30	        assert z.size == M
    31	        for w in z:
    32	            ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
    33	        ax.plot(z.real, z.imag, 'o', ms=2.5, color='tab:red', alpha=0.85, linestyle='none')
    34	        cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
    35	        for (a, b), c in cnt.items():
    36	            if c > 1:
    37	                ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3),
    38	                            fontsize=5, color='black')
    39	        r = float(np.abs(z).max())
    40	        lim = r * 1.15 if r > 0 else 1.0
    41	        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
```

- 31-32行: 原点からの線分、33行: 点、34-38行: 式26の重複注記、39-41行: 式27の
  スケール規約（目盛は実値のまま。46行の `ticklabel_format` は表記の指数化のみ）。

#### 4.3.3 ループ処理（2）: 最大角クラスター拡大（式28）

```python
    67	    z = load(N, 500)
    68	    amp = float(np.abs(z).max())
    69	    coarse = {}
    70	    for w in z:
    71	        key = (round(float(w.real) / amp, 2), round(float(w.imag) / amp, 2))
    72	        coarse.setdefault(key, []).append(w)
    73	    mem = max(coarse.values(), key=len)
    74	    zz = np.array(mem)
    75	    c = zz.mean()
    76	    dev = np.abs(zz - c)
    77	    spread = float(dev.max())
    78	    win = spread * 1.4 if spread > 0 else amp * 1e-12
...
    90	    ax.set_title(f'N={N}: {len(zz)} waves, dev={spread:.2e} (|z|max={amp:.2e})', fontsize=7)
```

- 70-73行: 式28の群化と最大クラスター選択、75-78行: 中心・割れ幅・窓幅、
  90行: パネル題への本数・spread・amp の印字（§6 の観察値の出所はこの印字である）。

#### 4.3.4 終了処理・例外的挙動（数式化しない要素）

- 40行 `lim = r*1.15 if r > 0 else 1.0`: 全点が原点の退化ケースへのフォールバック
  （本データでは発生しない）。
- 78行 `win = spread*1.4 if spread > 0 else amp*1e-12`: クラスターが1点のみ
  （spread=0）の場合の窓幅フォールバック（N=3 で発生: 最大クラスターが1波）。
- 50-51行・91-92行: 8×5=40 パネル中、未使用の 39・40 番目を `axis('off')`。
- 例外処理は assert（20・30行）のみで、それ以外の分岐はない。

## 5. 実行結果

### 5.1 再現コマンド

```bash
cd N3_N40_stage123_sweep_20260905
python3 plot_complex_plane_N3_N40_stage123_v1.py
# または ./run_all.sh の最終段として実行される
```

### 5.2 実行環境

- Python 3.9.6（`.venv/bin/python3`）、numpy 2.0.2、matplotlib（Agg 描画）
- macOS 26.3.1（arm64）

### 5.3 実行時間

3図の生成合計で**約1分以内**（支配項は 38×2 回の npz 読込と 780 本×38 パネルの線分描画）。

### 5.4 検証ゲート

| ゲート | 合格条件 | 実測 | 合否 |
|---|---|---|---|
| G1: 入力整合 | 各読込で denominator==N・steps==500・z.size==M の assert 通過 | 全38N×3図で通過（AssertionError なし） | **PASS** |
| G2: 完走 | `ALL DONE` 出力・PNG 3枚生成 | 確認 | **PASS** |

（入力 npz 自体の bit 系譜は第2章 G1〔全228走行の Z[0] が静的親と bit 一致〕で担保済み）

### 5.5 データ

| 項目 | 内容 |
|---|---|
| 入力 | `results/hm_N{N}_den_{N}_states_500.npz` × 38（第2章成果物、読み出しのみ） |
| 出力先 | パッケージ直下 |
| step0 図 | `fig_complex_plane_step0_N3_N40_stage123.png`（632,949 bytes） |
| 終了時図 | `fig_complex_plane_final_N3_N40_stage123.png`（4,582,367 bytes） |
| 拡大図 | `fig_complex_plane_final_zoom_N3_N40_stage123.png`（494,906 bytes） |
| SHA256 | 同梱 `SHA256SUMS.txt` を正本とする |

### 5.6 図化

上記3枚（8×5 グリッド、各パネル equal aspect・実値目盛）。図の読み方と
インフレーション図との対応は第2章 §2.6 の表の通り。

## 6. 実行分析（客観的報告と観察のみ。数値の出所は図中印字＝本プログラムの出力）

1. **step0 図**: 全 N=3..40 のパネルで、対蹠2ペア（4束）の星型が現れる。束は動径方向に
   振幅の広がりを持つ。式26の重複注記は小 N の一部（N=4 の x4、N=5〜9 の x2〜x3 など）
   に現れ、N が大きくなると厳密重複は消える。
2. **終了時図**: 全パネルで星型は消失し、ほぼ等振幅のリング状配置（小 N では疎な
   スポーク）になる。式26（12桁丸め）の重複注記は小 N の少数（N=4, 5 の x2）を除き
   現れない。
3. **拡大図**: 最大角クラスターの本数は 1〜10（N=3: 1波、N=39, 40: 10波）、割れ幅
   spread は N≥6 で 10⁻⁷〜10⁻⁴ のオーダー（例: N=40 で dev=2.04e-04、|z|max=3.64e-02、
   相対 ~5.6×10⁻³）。N=4, 5 のクラスター（2波）のみ dev が 10⁻¹⁶〜10⁻¹⁰ と機械精度
   レベルで一致する。
4. 観察のまとめ（事実のみ）: 終了時の角クラスターは、N≥6 の全域で「厳密な同一複素数
   への凝縮」ではなく有限幅の束である。第2章 §2.6 の役割分担の通り、これは H⊥/H が
   測らない終状態の微細組織の記録である。

---
（第3章おわり。総括論文は別途）
