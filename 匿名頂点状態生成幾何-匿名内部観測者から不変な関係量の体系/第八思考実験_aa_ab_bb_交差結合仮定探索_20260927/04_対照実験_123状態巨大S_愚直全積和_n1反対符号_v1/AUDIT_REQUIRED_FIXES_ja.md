# 要修正項目：第八思考実験 対照実験「123状態・巨大S・愚直全積和」v1 コード監査記録

**監査日:** 2026-09-27
**監査者:** Claude（Anthropic）
**方法:** ソースの読解と算術。監査対象のプログラムは実行していない。
**改訂:** 2026-09-27 ビット単位の一致に関する記述を削除。項目 1 を書き換え、項目 10 を追加（前回と異なる実装を、前回と同じ実装に戻す指示）。

---

## 監査対象

| ファイル | 行数 | SHA-256 |
|---|---:|---|
| `run_paper8_control_full123_bruteforce_n1_opposite_v1.py` | 750 | `2d3a3b7662b27c1bb95301060eddef8a54bc206f8f05609ac2bd9adf158d663f` |
| `show_progress_full123_v1.py` | 55 | `7c1d257779f719bd3240c5227198884a16ba606822e61b58f030de70d3bd931e` |
| `README_ja.md` | 79 | `874afcbd38663edbb7438ce7ea10eaa688127398f164f2f9077159201677c239` |

以下、行番号 `L…` は断りがない限り `run_paper8_control_full123_bruteforce_n1_opposite_v1.py` の行番号。

## 照合に使った参照

| ファイル | 場所 | SHA-256 |
|---|---|---|
| `run_paper7_self_interaction.py` | 第七思考実験 `01_独立自己相互作用数値実験_v1/` | `9f421cb2673ff4f012d9900189d6e3f32912a559036d96f963cce920fc523de3` |
| `paper7_reference_transition.py`（上と同一内容） | 第八思考実験 `01_最小共面位相差結合_8A_v1/` | `9f421cb2673ff4f012d9900189d6e3f32912a559036d96f963cce920fc523de3` |
| `run_strict_charged8_5cases.py` | 第七思考実験 `01_独立自己相互作用数値実験_v1/` | `78e3cf5fa37ff880fb3a983f9eca31f6dfa0546378a699b5053ca0d1c38b931b` |
| `paper6_init_G_q0_c1.py` | 第六思考実験 `04_補講_Gq0c1_初期化同値検証_20260927/` | `bf5435c9b753158f283f695a2b85d0576acb7e4732bb2cb121030c365ac907fe` |
| `run_paper8_control_block_diagonal.py` | 第八思考実験 `03_対照実験_交差項ゼロ_block_diagonal_v1/` | `95c5abd3ba7d707cb39f6a45f9397654a45e5fe1518ccfe7228dfb094b5318b4` |

## 判定の基準にした指示

1. 状態は `psi[123]` 1 本。毎 microstep、`S[123,123]` の 15,129 セルを全部作り、積和 15,129 項を省略なく計算する。
2. microstep の状態と microstep ごとの情報は保存しない。初期状態と全 macrostep の `psi[123]` は間引かず全部保存する。
3. 条件は n=1 反対符号の 1 条件。
4. 初期化は、第六思考実験 補遺の規格化 G=q0=c=1 と、補遺の初期化コードに準じる。
5. 進捗をバックグラウンドジョブから読み出せる。
6. 実行はコード監査の後。
7. 対照実験である。指示 1〜5 で決めた箇所のほかは、前回（第六・第七思考実験）のプログラムと同じ実装にする。

---

## 要修正項目 一覧

| 番号 | 項目 | 種別 |
|---:|---|---|
| 1 | 停止条件が前回と違う | 前回と異なる実装 |
| 2 | 初期化が補遺（G=q0=c=1）に準じていない | 指示と不一致 |
| 3 | S の構成が実装前の説明と違う | 要決定 |
| 4 | 積項カウンタが実測でなく計算値 | 表示と実体の不一致 |
| 5 | docstring「状態ベクトルは 1 本だけ」と実装の不一致 | 記述の誤り |
| 6 | 再実行すると raw を無警告で上書きする。再開もできない | データ保全 |
| 7 | 既定の出力先が Google Drive の中 | データ保全 |
| 8 | 強制終了時に、flush 前の macrostep 行を失う | データ保全 |
| 9 | 例外停止時の進捗が `micro=0` と表示される | 表示の誤り |
| 10 | 前回にない停止 IF 文が追加されている | 前回と異なる実装 |

項目 1 と項目 10 は、前回と異なる実装になっている箇所である。このままでは対照実験にならない。前回と同じ実装に戻す。

---

## 1. 停止条件が前回と違う

**該当行:** L420-425（`all_channels_reached_target`）、L648（`while not all_channels_reached_target(psi):`）、README L30-33

**前回の実装**

第六思考実験 基準ケース `run_strict_charged8_final.py` L173-178。

```python
while True:
    t,x,y=readout(z)
    w.writerow([nstep,z[P],...])
    if z[P] <= 20.0: break
    prev=z.copy(); z=macrostep(z); nstep += 1
```

第六思考実験 5 条件 `run_strict_charged8_5cases.py`（L136-139、L200）と、第七思考実験 `run_paper7_self_interaction.py`（L19、L136-138）も同じ停止文である。

**今回の実装**

```python
def all_channels_reached_target(psi):
    for ch in range(NCH):
        base = CHANNEL_BASES[ch]
        if float(psi[gidx(base, P)]) > TARGET_P:
            return False
    return True
...
        while not all_channels_reached_target(psi):
```

**違い**

| | 前回 | 今回 |
|---|---|---|
| 止まる条件 | P≤20 になったら止まる | aa・ab・bb の 3 本すべてが P≤20 になるまで止まらない |
| P≤20 になった後 | 走らない | ab は P≤20 の後も更新が続く |

**違いの結果**

L130 の `rate_dp` の式を P について求積した値。実測ではない。

| 出来事 | macrostep | その時の aa/bb の P |
|---|---:|---:|
| ab が P=20 を通過 | 488,344 | 43.32 |
| ab が P=0 に到達 | 約 573,108 | 41.99 |
| aa/bb が P=20 に到達 | 1,456,530 | 20 |

- 同じ求積は、実測値 488,345（第六思考実験 n1_attractive）と 1,456,531（第七思考実験 self_n1）を 1 step 差で再現する。
- ab は P≤20 の後も走り続け、約 573,108 macrostep で P≤0 になる。そこで項目 10 の IF 文が例外を出し、`status=failed` で終了する。aa/bb は P=20 に達しない。
- 前回の実験データでは、ab は 488,345 macrostep、aa/bb の条件は 1,456,531 macrostep で、どちらも P≤20 になった最初の行で止まっている。P=20 より下を走ったデータはない。

**必要な修正**

1. `all_channels_reached_target`（L420-425）と、L648 の `while not all_channels_reached_target(psi):` を削除する。
2. 前回と同じ停止文 `if P <= 20.0: break` に戻す。
3. README の「停止条件」の節（L30-33）を、修正後の停止文に合わせて書き直す。

**未決定の点**

- 123 状態には P が 3 つある（aa、ab、bb）。前回の停止文をどの P に適用するかは木原が決める。
- 最初に P≤20 になるのは ab で、488,345 macrostep。その時点の aa/bb の P は 43.32。

---

## 2. 初期化が補遺（G=q0=c=1）に準じていない

**該当行:** L105（`CHANNEL_CD`）、L79-81（`P0`、`TARGET_P`、`N0`）、L150-172（`initialize_psi`）、L448-454（HDF5 属性）、L612-629（`run_info`）

**事実**

| 補遺の初期化コード `paper6_init_G_q0_c1.py` にあるもの | 今回のプログラム |
|---|---|
| 規格化の宣言 `G = 1.0`、`Q0 = 1.0`、`C_LIGHT = 1.0`（L19-21） | なし |
| `M_OVER_Q0_EXACT = Fraction(20, 3)`、`M_TOTAL`、`M_A`、`M_B`（L26-30） | なし |
| 入力が電荷倍数 n と符号 `(n, sign_a, sign_b)`（L33-39） | なし |
| λ を有理数で作る（L50-51） | なし |
| `C = 1 − λ_A λ_B`、`D = (λ_A − λ_B)²` を計算し、最後に float へ変換（L52-54） | なし。完成値を 10 進で直書き |
| 状態ベクトルへの格納（UR、P、HR、Q、N、C、D、PS、PM、QB） | 同じ成分に同じ値 |

- 今回の初期化は `CHANNEL_CD = ((0.91, 0.0), (1.09, 0.36), (0.91, 0.0))` の直書きだけ。
- n、符号、λ、q0、M/q0 はコードのどこにも出てこない。
- README、`run_info.json`、HDF5 属性にも G=q0=c=1 の記載がない。
- 完成した C、D を表から与える形は、補遺より前の第六思考実験のやり方である。

**初期値**

補遺と同じ式で作った値と、直書き値は同じである。

| チャネル | 補遺の式への入力 | C | D |
|---|---|---:|---:|
| aa | n=1、(+, +) | 0.91 | 0.0 |
| ab | n=1、(+, −) | 1.09 | 0.36 |
| bb | n=1、(−, −) | 0.91 | 0.0 |

N は 0.25 で同じ。

**必要な修正**

1. L105 の直書きを削除する。
2. 規格化 G=q0=c=1、M/q0=20/3 と、入力 `(n, s_A, s_B) = (1, +1, −1)` を明示する。
3. 補遺の `case_relations` と同じ手順で、有理数のまま λ、C、D を作り、float への変換は最後に 1 回だけ行う。
4. チャネルごとの符号の組は次のとおり。自己関係は第七思考実験 `self_relations`（L143-147）の `C = 1 − λ²`、`D = 0` と同じ値になる。

| チャネル | 符号の組 |
|---|---|
| aa | (s_A, s_A) |
| ab | (s_A, s_B) |
| bb | (s_B, s_B) |

5. G、q0、c、M/q0、n、符号、λ_A、λ_B を `run_info.json` と HDF5 属性に記録する。

---

## 3. S の構成が実装前の説明と違う

**該当行:** L336-359（`interaction_entry`）、L181-328（`candidate_value`）

**事実**

- 0 以外になり得るセルは、1 チャネルあたり 341 個。内訳は「行 0-29 × q 列 30-40」の 330 個と、q の巡回の 11 個。
- 3 チャネルで 1,023 個。残りの 14,106 個は、状態の値によらず常に 0。
- q 以外の状態成分（P、U、H、Q、N、C、D、work register）の列は、すべての行で 0。S がこれらの成分に 0 以外を掛けることはない。
- 力学の計算（U の回転、P への増分の加算、RK4 の各段）は、すべて `candidate_value` の中のスカラー演算で行われる。積和がするのは、候補値に q を掛けて足すことだけである。
- 実装前の説明は「セルごとに 11 候補 `q0 F0 + … + q10 F10` を計算し、そのセルの係数を決めて `S_ij` に入れる」だった。実装は「候補値を q 列に 1 つずつ置く」で、構成が違う。
- 第六・第七思考実験の `transition(z)` の one-hot 積和 `out[i] = Σ_k q[k] · cand[i,k]` とは、同じ構成である。

**必要な対応**

- 現状の構成を認めるか、状態成分の列に係数を置く構成へ作り直すかを決める。決めるのは木原。

---

## 4. 積項カウンタが実測でなく計算値

**該当行:** L541、L551

**事実**

- `logical_product_terms_visited` は `完了 microstep 数 × 15,129` の計算値である。
- 積和のループ（L377-381）の中では何も数えていない。
- 終了時に期待値と照合して FAIL にする処理もない。
- この値は、積を実際に計算したことの証拠にならない。

**必要な対応**

- 次のどちらかにする。どちらにするかは木原が決める。
  - 積を 1 回計算するごとに 1 を加える実測カウンタに替え、終了時に `microstep 数 × 15,129` と照合して不一致なら FAIL にする。
  - この項目を `progress.json` と stdout から削除する。

---

## 5. docstring「状態ベクトルは 1 本だけ」と実装の不一致

**該当行:** L18-21（docstring）、L604（`psi_next`）

**事実**

- docstring は "One persistent state vector only" と書いている。
- 実装には `psi[123]` のほかに `psi_next[123]` がある。同期更新（旧状態だけから新状態を作る）のための作業配列で、毎 microstep 全成分が上書きされる。

**必要な修正**

- docstring の規則 1 を「状態ベクトルは `psi[123]` 1 本。同期更新のための作業配列 `psi_next[123]` を 1 本持つ」に書き換える。README の該当箇所も同じに直す。

---

## 6. 再実行すると raw を無警告で上書きする。再開もできない

**該当行:** L447（`h5py.File(path, "w")`）、L603-608（毎回初期状態から開始）

**事実**

- 出力フォルダに `raw_macro_part000.h5` がある状態で再実行すると、警告なしに上書きする。
- 途中から再開する機能がない。止まった場合は macrostep 0 からやり直しになる。

**必要な修正**

1. 出力フォルダに `raw_macro_part*.h5` が既にある場合は、開始せずにエラーで止める。
2. 再開機能を付けるかどうかは木原が決める。全 macrostep の `psi[123]` が保存されているので、最後に保存された行から再開する材料はある。

---

## 7. 既定の出力先が Google Drive の中

**該当行:** L729-733

**事実**

- `--out` の既定値はプログラムと同じフォルダの `results_full123_bruteforce_n1_opposite_v1`。
- プログラムは Google Drive の同期フォルダにあるので、実行中の HDF5（1 part 約 99 MB）と、30 秒または 100 macrostep ごとに置き換わる `progress.json` が同期の対象になる。

**必要な修正**

- 実行時は `--out` で Google Drive の外のローカルフォルダを指定する。完了後に結果をコピーする。

---

## 8. 強制終了時に、flush 前の macrostep 行を失う

**該当行:** L677-678（100 macrostep ごとの flush）、L108-111（定数）

**事実**

- 通常の例外では `finally`（L710-711）でファイルを閉じるので、行は失われない。
- 強制終了や電源断では、最後の flush より後の最大 99 行が失われる。書き込み中の part は読めなくなる場合がある。

**必要な修正**

- flush の間隔を決め直す。毎 macrostep にするか、現状の 100 のままにするかは木原が決める。

---

## 9. 例外停止時の進捗が `micro=0` と表示される

**該当行:** L699-707

**事実**

- macrostep の途中で例外が出た場合も、失敗時の進捗は `micro_in_macro=0` で出力される。止まった microstep の位置が記録に残らない。

**必要な修正**

- 直近に完了した microstep の番号を変数に保持し、失敗時の進捗にその値を渡す。

---

## 10. 前回にない停止 IF 文が追加されている

**該当行:** L121-123（`_require_rate_domain`）、L127、L138（その呼び出し）、L403-417（`validate_macro_boundary`）、L671（その呼び出し）

**前回の実装**

第七思考実験 `run_paper7_self_interaction.py` L56-60。第六思考実験の 2 本（`run_strict_charged8_final.py` L46、`run_strict_charged8_5cases.py` L42）も同じ。

```python
def _rates(p,n,c,d):
    rp=math.sqrt(p); rc=math.sqrt(c)
    dp=-(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
    dt=p*rp/rc
    return dp,0.0,dt
```

- 値を調べる IF 文はない。
- 前回の 3 本で、P の値によって走行を止める文は、項目 1 の停止文だけである。

**今回の実装**

追加 1。変化率の計算のたびに呼ばれる。

```python
def _require_rate_domain(p, c, label):
    if not (math.isfinite(p) and math.isfinite(c) and p > 0.0 and c > 0.0):
        raise ValueError(f"rates domain failure at {label}: p={p!r}, c={c!r}")

def rate_dp(p, n, c, d):
    _require_rate_domain(p, c, "rate_dp")
    ...
def rate_dt(p, _n, c, _d):
    _require_rate_domain(p, c, "rate_dt")
    ...
```

追加 2。macrostep ごとに呼ばれる。

```python
def validate_macro_boundary(psi):
    for i in range(NTOT):
        if not math.isfinite(float(psi[i])):
            raise FloatingPointError(...)
    for ch in range(NCH):
        ...
        if p <= 0.0 or c <= 0.0 or ps <= 0.0:
            raise ValueError(...)
        ...
        if qs != 1.0 or qnz != 1:
            raise ValueError(...)
```

**必要な修正**

1. `_require_rate_domain`（L121-123）を削除する。L127 と L138 の呼び出しも削除する。
2. `validate_macro_boundary`（L403-417）を削除する。L671 の呼び出しも削除する。
3. `rate_dp`・`rate_de`・`rate_dt` の式は前回の `_rates` と同じなので、変更しない。

---

## 実行前の注意（コードの修正ではない）

- 実行時間は未計測。1 microstep で Python の関数呼び出しが 15,129 回、積が 15,129 回ある。概算で 1 macrostep 約 0.2 秒。488,345 macrostep で約 1 日、1,456,531 macrostep で約 3 日になる。

---

## 報告どおりだった項目

| 項目 | 根拠 |
|---|---|
| `psi[123]` を 1 回だけ確保 | L603 |
| `S[123,123]` を 1 回だけ確保 | L605 |
| 41 状態部分ベクトル、41×41 小行列、reshape、concatenate なし | スライスは保存処理の L481 だけ |
| argmax なし。11 候補を q 列に構成 | L350-354 |
| 毎 microstep、S の 15,129 セルに代入 | L366-368 |
| 積和 15,129 項、値によるスキップなし | L377-381 |
| 1 macrostep = 11 microstep | L650 |
| microstate を保存しない | — |
| 初期状態と全 macrostep の `psi[123]` を保存 | L633、L674 |
| 参照 raw を読まない | 開くのは出力の HDF5 と JSON だけ（L447、L568） |
| stdout の flush と、原子的に更新する `progress.json` | L566-590 |
| 30 秒経過で macro 途中も出力 | L657-668 |
| 禁止語（reshape、concatenate、argmax、np.dot、`S @ psi` など）がない | 出現はコメントと docstring だけ |
| 比較プログラムは未作成、本体は未実行 | フォルダは 3 ファイル、出力フォルダなし |

- 生成経路の `if` は、行番号・位相番号・チャネル番号による分岐だけである。値を見るのは L122 の定義域チェックだけである（項目 10）。
- 第七思考実験の `transition(z)` と、30 行 × 11 位相の候補、q の巡回、初期状態の格納先がすべて一致した。

---

## 関連：補遺の初期化コードで確認した点

対象は `paper6_init_G_q0_c1.py`（公開済み。Version DOI 10.5281/zenodo.22985300）。

| 該当行 | 事実 |
|---|---|
| L50-51 | λ は `Fraction(sign*3*n, 10)` と直接書かれている。`M_OVER_Q0_EXACT` からは計算していない |
| L23-25 | M/q0=20/3 から λ=3n/10 への換算は、コメントの中だけにある |
| L42-44 | `charge_to_mass_from_multiple` は定義だけで、どこからも呼ばれていない |
| L19-21、L28-30 | `G`、`Q0`、`C_LIGHT`、`M_TOTAL` は記録用の `meta` に入るだけで、C、D の計算に使われていない |

項目 2 の修正で補遺の `case_relations` をそのまま写すと、この構成も引き継ぐ。λ を M/q0 から計算する形にするかどうかは木原が決める。

---

## 参考：旧対照（03）で確認した箇所

対象は `run_paper8_control_block_diagonal.py`。停止条件を調べるために読んだ該当行だけの記録である。

| 該当行 | 事実 |
|---|---|
| L234 | 積和が `build_S123(psi)@psi` で、numpy の行列積を使っている |
| L235-236 | 生成ループの中で `independent_microstep(ref)` を並走させ、毎 microstep 比較している |
| L239-244 | 保存は macro 1・10・100・1000・4000 の 5 点の P と、最終状態だけである |
| L232 | 走行長は 4000 macrostep 固定 |

---

## 再現用コード

項目 1 の数値を出した算術。監査対象のプログラムを読み込まない。

### 項目 1：ab が P=0 に達する macrostep

```python
import math
H = 2.0*math.pi/4000
def steps(c, d, p_hi, p_lo, n=0.25, m=2_000_000):
    a = (4.0/3.0)*n*d*math.sqrt(c); b = (64.0/5.0)*n*c*math.sqrt(c)
    f = lambda p: p**1.5/(a*p + b)
    h = (p_hi - p_lo)/m
    s = f(p_lo) + f(p_hi) + sum((4 if k % 2 else 2)*f(p_lo + k*h) for k in range(1, m))
    return s*h/3.0/H
ab_50_20 = steps(1.09, 0.36, 50, 20)
ab_20_0  = steps(1.09, 0.36, 20, 0)
aa_50_20 = steps(0.91, 0.0, 50, 20)
t = ab_50_20 + ab_20_0
b = (64.0/5.0)*0.25*0.91*math.sqrt(0.91)
print(round(ab_50_20), round(aa_50_20), round(ab_20_0), round(t))
print((50**2.5 - 2.5*b*t*H)**0.4, (50**2.5 - 2.5*b*ab_50_20*H)**0.4, t/aa_50_20)
```

出力：`488344 1456530 84764 573108`、`41.99…  43.32…  0.393…`
