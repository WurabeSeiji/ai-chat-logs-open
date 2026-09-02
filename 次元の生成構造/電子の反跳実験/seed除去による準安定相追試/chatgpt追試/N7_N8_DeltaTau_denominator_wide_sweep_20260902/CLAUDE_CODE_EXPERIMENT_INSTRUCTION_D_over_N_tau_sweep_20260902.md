\# Claude Code 実験指図：N=7,8 における Δτ \= 2π/D の広域分母スイープ

作成日: 2026-09-02    
目的: 頂点数 N と位相刻み Δτ \= 2π/D の関係を、D を広範囲に独立走査して検証する。    
重要: 本実験では辺数に既に M=N(N-1)/2 を使っているため、位相刻みの分母は \*\*D\*\* と表記する。D と辺数 M を混同しないこと。

\---

\#\# 0\. この実験で判定したいこと

既存の N=3..16 分母制御実験では、D=N-2,N-1,N,N+1,N+2 と固定 D=124 を比較したところ、

\- D が小さいほど Hperp/H の巨視化が早い  
\- D=124 は 500 step 内ではほぼ数値床に留まる  
\- しかし step は D ごとに異なる物理的／位相的時間を表すため、D=124 が本当に特殊なのか、単に Δτ が小さく観測時間が足りなかっただけなのか未判定

である。

今回の実験では、N を固定して D を広範囲にスイープし、次を切り分ける。

1\. \*\*単なる時間刻み効果か\*\*  
   \- τ \= step × 2π/D で曲線が重なるなら、D 効果の大部分は時間再尺度化。

2\. \*\*D/N が支配変数か\*\*  
   \- N=7 と N=8 の結果が D/N で重なる、または共通特徴を持つか。

3\. \*\*D=N に固有の力学的特徴があるか\*\*  
   \- D=N を事前に特別扱いせず、広域スイープ後にピーク・谷・屈曲・安定窓などが D/N=1 に現れるかを判定する。

4\. \*\*124 が特殊数なのか、単に N から遠いのか\*\*  
   \- D=124 を同一 τ 観測窓まで走らせてもフラットか。  
   \- 同一 τ で他の大きな D と同様なら、124 固有性を主張しない。

5\. \*\*奇数 N / 偶数 N の差\*\*  
   \- N=7（奇数・距離クラス）と N=8（偶数・1-factor）を最小対照として比較する。

\#\#\# 禁止する先入観

\- D=N が特別だと仮定して結果を選別しない。  
\- D=124 が特別だと仮定しない。  
\- 逆に「D は単なる時計」と最初から決めつけない。  
\- 結果を既存仮説に合わせるための平滑化、正規化、クリップ、seed 追加をしない。

\---

\#\# 1\. 作業場所

この指図ファイルが置かれた新規フォルダを実験ルートとする。

既存の正本データ・正本プログラムは\*\*読み取り専用\*\*で参照し、上書きしないこと。

親状態の参照元は、実験ルートから見て原則として以下。

\`\`\`text  
../干渉保存力学\_資格審査とシード無し系列\_20260831/  
    hm\_mp\_free\_N3\_N40\_20260901/  
        data/hm\_N7/parent\_v.npz  
        data/hm\_N8/parent\_v.npz  
\`\`\`

参照実装として、以下の分母制御実験および N=40 実験も読むこと。

\`\`\`text  
../干渉保存力学\_資格審査とシード無し系列\_20260831/  
    ChatGPT\_denominator\_controls\_64bit\_with124\_20260901/

../干渉保存力学\_資格審査とシード無し系列\_20260831/  
    N40\_L40000\_DeltaTau\_2pi\_N\_20260901/  
\`\`\`

\*\*注意:\*\* 相対パスが実環境で異なる場合は、同名正本を検索して実体を確認し、README/RUN\_METADATA に実際の参照パスと SHA256 を記録する。推測で別ファイルを代用しない。

\---

\#\# 2\. 固定する力学

既存の「干渉保存力学」から変更しない。

辺数:

\\\[  
M\_{\\rm edge}=\\frac{N(N-1)}{2}  
\\\]

隣接行列:

\- 2 本の辺が頂点を共有する場合 A\_ef=1  
\- 同一辺は 0  
\- その他は 0

干渉行列:

\\\[  
H\_{ef}(z)=A\_{ef}\\,\\overline{z\_e}z\_f  
\\\]

1 step:

\\\[  
z'=\\exp\\left\[-i\\left(\\frac{2\\pi}{D}\\right)H(z)\\right\]z  
\\\]

実装条件:

\- \`numpy.float64\`  
\- \`numpy.complex128\`  
\- \`numpy.linalg.eigh(H)\`  
\- 1 step ごとに現在の z から H を再構成  
\- seed 追加なし  
\- clipping なし  
\- renormalization なし  
\- 人工ノイズなし  
\- IF による状態変更なし  
\- 高精度演算へ勝手に変更しない  
\- 高速化のために力学を書き換えない

\`H\_total \= z†z\` は記録する。丸め誤差を隠すための再正規化は禁止。

\---

\#\# 3\. 初期状態

主実験では、保存済み hm 親をそのまま使用する。

\- N=7: \`hm\_N7/parent\_v.npz\["v"\]\`  
\- N=8: \`hm\_N8/parent\_v.npz\["v"\]\`

各 N について初期状態 z0 から初期実 2 平面を一度だけ定義する。

\`\`\`python  
p \= z0.real / ||z0.real||  
q \= z0.imag \- dot(z0.imag, p) \* p  
q \= q / ||q||  
\`\`\`

各時刻で

\\\[  
z\_\\perp=z-p(p\\cdot z)-q(q\\cdot z)  
\\\]

\\\[  
f\_\\perp=\\frac{H\_\\perp}{H}  
       \=\\frac{z\_\\perp^\\dagger z\_\\perp}{z^\\dagger z}  
\\\]

を主観測量とする。

\#\#\# 初期状態監査

走行前に各 N について最低限以下を記録する。

\- N  
\- M\_edge  
\- norm  
\- mean\_amp2 \= ||z||²/M\_edge  
\- H\_total  
\- Hperp/H at step 0  
\- abs(z@z)/H  
\- parent file SHA256  
\- dtype

既存親を再生成しないこと。

\---

\#\# 4\. D スイープ

\#\#\# Stage A: 広域・密スイープ（探索）

各 N=7,8 について、

\\\[  
D=2,3,4,\\ldots,256  
\\\]

の \*\*全整数\*\*を走らせる。

各 D:

\- 500 step  
\- step 0 を含め 501 点記録  
\- この Stage A は広域構造・局所的な窓・異常点の探索用

さらに遠方確認として、

\`\`\`text  
D \= 320, 384, 512  
\`\`\`

も 500 step 走らせる。

\*\*Stage A の 500 step の結果だけで D が安定だと判定しないこと。\*\*  
D が大きいほど同じ step 数で進む τ が短いためである。

\#\#\# Stage B: 同一 τ 観測窓（本比較）

次の D を重点集合とする。

\`\`\`text  
D \= 2..32 の全整数  
    40, 48, 64, 80, 96, 112, 124, 128,  
    160, 192, 224, 256, 320, 384, 512  
\`\`\`

各 N ごとに、基準 D=N を 500 step 走らせた累積位相時間と同じところまで走らせる。

\\\[  
\\tau\_{\\max}(N)=500\\frac{2\\pi}{N}  
\\\]

各 D の必要 step 数は

\\\[  
S\_{\\max}(N,D)=  
\\left\\lceil  
500\\frac{D}{N}  
\\right\\rceil  
\\\]

とする。

これにより D=124 も「500 step」ではなく、N=7,8 の D=N 系と\*\*同一 τ\*\*まで観察する。

\#\#\# Stage C: 自動追試

Stage A/B から以下が見つかった場合のみ追加する。

\- τ\_onset の局所極小・極大  
\- γ\_tau の局所極小・極大  
\- D/N=1 近傍の不連続・鋭い屈曲  
\- 飽和値の異常窓  
\- 124 だけが近傍 D と明瞭に異なる

整数 D のため、異常 D の前後 ±8 程度を再走行し、環境情報・入力 SHA を固定して再現確認する。

\---

\#\# 5\. 横軸の定義

\#\#\# 主図で step を横軸にしない

主たる時間軸は

\\\[  
\\boxed{  
\\tau=\\frac{2\\pi}{D}\\,{\\rm step}  
}  
\\\]

とする。

同時に以下も CSV に記録する。

\\\[  
C=\\frac{{\\rm step}}{D}  
\\\]

（周期数、\\(\\tau=2\\pi C\\)）

N 間比較用の分母軸:

\\\[  
\\boxed{x\_D=\\frac{D}{N}}  
\\\]

さらに相互作用スケールを含む補助無次元時間として、

\\\[  
r\_N^2=\\frac{||z\_0||^2}{M\_{\\rm edge}}  
\\\]

\\\[  
\\boxed{  
\\chi=  
2r\_N^2(N-2)\\tau  
}  
\\\]

も計算する。

\`step\` は再現性監査用の列として保持するが、主図の横軸には使わない。

\---

\#\# 6\. 各 run で保存する時系列

最低限:

\`\`\`text  
N  
D  
D\_over\_N  
step  
tau  
cycles  
chi  
Hperp\_frac  
H\_parallel\_frac  
H\_total  
H\_total\_rel\_drift  
global\_closure \= abs(z@z)/H  
PR  
PR\_over\_M  
amp\_min  
amp\_max  
amp\_std  
finite  
\`\`\`

可能なら初期平面外の rank / singular-value 指標も別解析で追加してよいが、  
\*\*主力学を変更してはならない。\*\*

\---

\#\# 7\. onset と成長率

\#\#\# Primary onset

既存実験との連続性のため、

\\\[  
f\_\\perp=H\_\\perp/H \> 0.05  
\\\]

を primary macroscopic onset とする。

保存する値:

\`\`\`text  
onset\_step\_0p05  
onset\_tau\_0p05  
onset\_cycles\_0p05  
onset\_chi\_0p05  
\`\`\`

\#\#\# 補助閾値

閾値依存性確認として、

\`\`\`text  
1e-6  
1e-3  
0.05  
\`\`\`

の 3 種を記録する。

\#\#\# 初期指数成長率

固定窓

\\\[  
10^{-12}\\le f\_\\perp\\le10^{-4}  
\\\]

に入る点を使い、

\\\[  
\\ln f\_\\perp \= a \+ \\gamma\_{\\rm step}\\,{\\rm step}  
\\\]

および

\\\[  
\\ln f\_\\perp \= a\_\\tau \+ \\gamma\_\\tau\\,\\tau  
\\\]

を最小二乗で求める。

必ず以下を記録する。

\- slope  
\- intercept  
\- R²  
\- 使用点数  
\- fit の最小/最大 step  
\- fit の最小/最大 τ

点数 5 未満なら \`NA\`。無理に fit しない。

理論上は

\\\[  
\\gamma\_\\tau=  
\\gamma\_{\\rm step}\\frac{D}{2\\pi}  
\\\]

なので数値的にもクロスチェックする。

\---

\#\# 8\. 飽和後の統計

ユーザーは飽和後の安定範囲を重視しているため、必ず定量化する。

primary onset が存在する run について、走行末尾 20% を tail window とし、

\`\`\`text  
sat\_mean\_Hperp\_frac  
sat\_std\_Hperp\_frac  
sat\_min\_Hperp\_frac  
sat\_max\_Hperp\_frac  
sat\_q05\_Hperp\_frac  
sat\_q95\_Hperp\_frac  
\`\`\`

を記録する。

ただし tail window 開始時点が onset より前の場合は \`NA\` とし、追加走行が必要と判定する。

飽和を「1 に近いから固定点」と解釈しない。  
Hperp/H が一定でも、closure、PR、振幅分布が内部で運動している可能性があるため、tail window でそれらも mean/std/min/max を保存する。

\---

\#\# 9\. 必須図

\#\#\# Fig 1: N=7 選択 D の Hperp/H vs τ

\- x: τ  
\- y: Hperp/H（log）  
\- selected D:

\`\`\`text  
4,5,6,7,8,9,10,12,16,24,32,64,124,256  
\`\`\`

\#\#\# Fig 2: N=8 選択 D の Hperp/H vs τ

同じ D 集合。

\#\#\# Fig 3: τ\_onset vs D/N

\- x: D/N  
\- y: onset\_tau\_0p05  
\- N=7 と N=8 を重ねる  
\- D/N=1 に細い参照線を置いてよい  
\- \*\*D/N=1 が特別だという注釈は結果を見てから\*\*

\#\#\# Fig 4: γ\_tau vs D/N

\- N=7 / N=8 overlay  
\- R² が低い fit は別マーカーまたは除外して明記

\#\#\# Fig 5: 飽和値 vs D/N

\- y: sat\_mean\_Hperp\_frac  
\- q05-q95 または min-max の帯  
\- N=7 / N=8

\#\#\# Fig 6: χ 軸での collapse test

selected D について

\- x: χ  
\- y: Hperp/H（log）

N=7 と N=8 の両方を比較し、τ では残った N 差が χ で消えるかを見る。

\#\#\# Fig 7: 広域 heatmap

各 N を別図として、

\- x: D/N  
\- y: τ  
\- color: log10(Hperp/H)

Stage A/B のデータから作る。  
補間で偽の細構造を作らない。原データ格子を尊重する。

\#\#\# 監査用 step 図

必要なら \`figures/audit\_step\_axis/\` に保存してよい。  
\*\*論理解釈の主図に step 軸を使わない。\*\*

\---

\#\# 10\. 最終集約 CSV

最低限以下を 1 行 1 (N,D) で出す。

\`results/sweep\_summary.csv\`

\`\`\`text  
N  
M\_edge  
D  
D\_over\_N  
r2bar  
stage  
steps\_run  
tau\_max  
initial\_Hperp\_frac

onset\_step\_1e-6  
onset\_tau\_1e-6  
onset\_step\_1e-3  
onset\_tau\_1e-3  
onset\_step\_0p05  
onset\_tau\_0p05  
onset\_cycles\_0p05  
onset\_chi\_0p05

gamma\_step  
gamma\_tau  
fit\_R2  
fit\_n

sat\_mean\_Hperp\_frac  
sat\_std\_Hperp\_frac  
sat\_min\_Hperp\_frac  
sat\_max\_Hperp\_frac  
sat\_q05\_Hperp\_frac  
sat\_q95\_Hperp\_frac

Htotal\_max\_rel\_drift  
closure\_tail\_mean  
closure\_tail\_std  
PR\_over\_M\_tail\_mean  
PR\_over\_M\_tail\_std

numpy\_version  
python\_version  
platform  
input\_sha256  
program\_sha256  
\`\`\`

Stage A と Stage B が同じ (N,D) を含む場合、  
summary では両方を区別できる \`stage\` 列を必ず持つ。

\---

\#\# 11\. 期待結果を「仮説」としてのみ記録する

以下は検証対象であり、結果として先に書かない。

\#\#\# H0: 単純な時間再尺度化

D を変えた違いが τ 軸でほぼ消える。

この場合:

\- 124 固有性なし  
\- 500 step でフラットだったのは短い τ しか見ていなかったため  
\- γ\_tau は D にあまり依存しない

\#\#\# H1: D/N 支配

N=7 と N=8 が D/N 軸で共通の構造を持つ。

\#\#\# H2: D=N 特異／準特異

D/N=1 に、

\- onset の局所極小・極大  
\- γ\_tau のピーク・谷  
\- 飽和値変化  
\- 安定窓境界

などが再現可能に現れる。

\#\#\# H3: 124 固有性

同一 τ まで走らせても D=124 のみ近傍の大きな D と異なる。

\*\*H3 は特に厳しく判定する。\*\*  
124 を単独比較せず、96,112,124,128,160 などの近傍・周辺と比較する。

\---

\#\# 12\. 数値監査

各 run で:

\- \`H-H†\` の相対ノルム  
\- \`H\_total\` drift  
\- finite  
\- dtype

を記録。

norm drift を補正しない。

推奨:

\- warning: relative norm drift \> 1e-10  
\- hard failure: non-finite  
\- relative norm drift \> 1e-7 は run failure として記録し停止可

ただし failure を隠して再正規化して続行しない。

同一入力・同一プログラムでも環境差で軌道が分岐し得るため、

\- Python version  
\- NumPy version  
\- OS/platform  
\- BLAS/LAPACK 情報（取得可能なら）  
\- CPU architecture

を \`RUN\_METADATA.json\` に保存する。

\---

\#\# 13\. 出力フォルダ構成

\`\`\`text  
.  
├── CLAUDE\_CODE\_EXPERIMENT\_INSTRUCTION\_D\_over\_N\_tau\_sweep\_20260902.md  
├── README.md  
├── RUN\_METADATA.json  
├── SHA256SUMS.txt  
├── ANALYSIS.md  
├── program/  
│   ├── run\_sweep.py  
│   ├── analyze\_sweep.py  
│   └── plot\_sweep.py  
├── data/  
│   ├── N7/  
│   │   └── Dxxx/  
│   └── N8/  
│       └── Dxxx/  
├── results/  
│   ├── sweep\_summary.csv  
│   ├── stageA\_summary.csv  
│   ├── stageB\_summary.csv  
│   └── anomaly\_followups.csv  
└── figures/  
    ├── fig01\_N7\_tau\_curves.png  
    ├── fig02\_N8\_tau\_curves.png  
    ├── fig03\_onset\_tau\_vs\_D\_over\_N.png  
    ├── fig04\_growth\_gamma\_tau\_vs\_D\_over\_N.png  
    ├── fig05\_saturation\_vs\_D\_over\_N.png  
    ├── fig06\_chi\_collapse.png  
    ├── fig07\_N7\_heatmap.png  
    ├── fig08\_N8\_heatmap.png  
    └── audit\_step\_axis/  
\`\`\`

ファイル名は同内容ならこの命名を優先する。  
既存フォルダへの出力は禁止。

\---

\#\# 14\. ANALYSIS.md に必ず回答する問い

結果を見た後、データに基づいて次を順番に回答する。

1\. τ 軸に直すと D ごとの曲線はどこまで重なるか。  
2\. D=124 は同一 τ でもフラットか。  
3\. D=124 は 96,112,128,160 と比べて統計的／視覚的に異常か。  
4\. N=7 と N=8 は D/N 軸で重なるか。  
5\. D/N=1 に局所特徴は存在するか。  
6\. γ\_tau は D に依存するか、D/N に依存するか、ほぼ一定か。  
7\. 飽和後 Hperp/H の平均・範囲は D で変わるか。  
8\. N=7 と N=8 の差は奇偶差と呼べるほど再現的か。それとも N 個別差に留まるか。  
9\. χ 軸にすると N=7/8 の差は縮小するか。  
10\. 現データだけで「D=N は自然な時間刻み」と言えるか。言えないなら何が不足か。

\#\#\# 表現上の注意

\- 「証明した」「必然」と書かない。  
\- 数値的に観測したことと、解釈を分ける。  
\- onset=0.05 は巨視化閾値であり、指数不安定性の開始時刻そのものではない。  
\- saturation は Hperp/H の飽和であり、状態 z 自体の固定点とは限らない。

\---

\#\# 15\. 完了条件

以下が揃うまで「完了」としない。

\- N=7,8 の Stage A 全 D が走行済み  
\- Stage B の重点 D が同一 τ まで走行済み  
\- \`sweep\_summary.csv\` が存在  
\- 主図が step ではなく τ / D/N / χ 軸で作成済み  
\- saturation range が定量化済み  
\- D=124 が「短時間だっただけか」を同一 τ 比較で回答済み  
\- D=N の特殊性について肯定・否定どちらにも先入観なく結論を記載  
\- SHA256 と環境情報を保存  
\- 既存正本を一切上書きしていない

\---

\#\# 16\. Claude Code への最終指示

まず既存正本の場所と入力 SHA256 を確認し、\*\*実行前に短い実験計画と推定計算量を表示すること\*\*。  
その後は質問待ちで止まらず、Stage A → 集約 → Stage B → 集約 → 図化 → ANALYSIS.md の順で実行する。

ただし、実行時間が想定を大きく超える場合は、力学を変更して高速化してはならない。  
どの Stage / D まで完了したかを保存し、再開可能な checkpoint 方式にすること。

最重要の判定は次の一文である。

\> \*\*D=124 がフラットだった理由は、124 自身の特殊性なのか、それとも N=7,8 に対して Δτ=2π/124 が小さく、従来の 500 step 観測窓が短すぎただけなのか。さらに D=N に、単純な時間再尺度化では消えない固有構造が存在するか。\*\*  
