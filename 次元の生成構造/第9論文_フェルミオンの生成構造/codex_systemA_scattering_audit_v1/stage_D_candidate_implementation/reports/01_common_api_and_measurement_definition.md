# 共通APIと二層測定定義

## モデル定義

### 物理散乱層

一次元カーネル \(a(u),b(u)\) を、搬送波と正規化ηモードによって

\[
D_X[f](u,\eta)
=
f(u)e^{iq_Xp_0u}
\frac{e^{im_X\eta}}{\sqrt{N_\eta}}
\]

と全状態へ埋め込む。実装値は
\((q_A,m_A)=(+1,1)\)、\((q_B,m_B)=(-1,2)\)、
\((N_u,N_\eta)=(512,16)\) である。

散乱の正本は

\[
\Psi_A^{\rm raw}=r_{\rm eff}\Psi_A+t_{\rm eff}\Psi_B,\qquad
\Psi_B^{\rm raw}=t_{\rm eff}\Psi_A+r_{\rm eff}\Psi_B
\]

であり、経路振幅、干渉項、rawノルムはすべてこの全状態で計算した。
`*_norm` は一貫して二乗L2ノルム \(\lVert\cdot\rVert^2\) を表す。

### パリティ読出し層

入力ごとに対応するηモードへ射影し、搬送波を逆シフトして
\(D_X^{-1}\Psi_X\) を得る。半周期移動 \(P_u f(u)=f(u+\pi)\) に対して

\[
C_\pi[f]=\langle f,P_uf\rangle,\qquad
c_\pi[f]=\frac{\Re C_\pi[f]}{\lVert f\rVert^2},
\]

\[
p_B=\frac{\lVert(f+P_uf)/2\rVert^2}{\lVert f\rVert^2},\qquad
p_F=\frac{\lVert(f-P_uf)/2\rVert^2}{\lVert f\rVert^2}
\]

を測定した。

Candidate 1・3の応答は、それぞれ

\[
F_1=\frac{c_\pi[a]+c_\pi[b]}{2},\qquad
F_3=c_\pi[ab]
\]

である。Candidate 3では
\(\lVert ab\rVert^2\le10^{-14}\lVert a\rVert^2\lVert b\rVert^2\)
を固定エラー条件とし、暗黙のゼロ返却を禁止した。

### raw出力パリティ

たとえばA出力は、全和を単一搬送波で復調せず、
\(r\Psi_A\) と \(t\Psi_B\) をそれぞれ元の \(D_A,D_B\) で復調する。
直交する由来空間の直和として、ノルム、\(C_\pi\)、偶奇成分ノルムの
分子を加算した。この定義では由来間の人工的な交差項を作らない。

### 半周期同変性

\[
P_A=D_AP_uD_A^{-1},\qquad P_B=D_BP_uD_B^{-1}
\]

を使い、入力と各出力由来経路を対応する作用素で移動した。
状態依存写像全体について

\[
\varepsilon_P=
\frac{\lVert
S(P_A\Psi_A,P_B\Psi_B)
-(P_A\oplus P_B)S(\Psi_A,\Psi_B)
\rVert}
{\lVert S(\Psi_A,\Psi_B)\rVert}
\]

を保存した。

## コード上の事実

- `ScatteringResult` は指示された全フィールドを保持し、復調残差、
  クリップ前角度、由来別再構成残差、交換残差用データを追加した。
- raw出力を物理的正本とし、チャネル別正規化後出力を別配列として
  保持した。
- A/B交換試験では、配列だけでなく復調仕様と期待カーネルも一緒に
  literal に交換した。
- `B_to_A_transfer` はStage Dでは計算していない。Stage Bの同名量は
  初期Bスペクトルに対するA出力の余弦類似度であり経路量ではない。
  Stage Dの経路量は `path_b_to_a_norm` である。

## 数学的帰結

- \(t=e^{i\theta}\cos\theta\)、
  \(r=-ie^{i\theta}\sin\theta\) から
  \(|r|^2+|t|^2=1\) と \(r^*t+t^*r=0\) が自動的に成立する。
- \(c_\pi[P_uf]=c_\pi[f]\)、
  \(c_\pi[(P_ua)(P_ub)]=c_\pi[ab]\) なので、C1・C3は指定した
  同時半周期移動の下で応答値を保つ。
- ηモードの直交により、現設定の二由来経路干渉は理論上ゼロである。

## 数値観測

入力復調残差の全量最大値は \(1.589\times10^{-16}\)、復調後の
再変調往復残差は \(1.627\times10^{-16}\)、raw出力由来別再構成残差は
\(2.154\times10^{-16}\) だった。

## 作業仮説

パリティ指標を散乱角へ結び付けること自体は、Stage Dで採用した
候補モデルであって既存System Aからの導出ではない。

## 未導出

復調読出しが実在する相互作用機構に対応するか、κと包絡ρの力学的起源、
反復散乱でraw状態と正規化状態のどちらを次状態にするかは未導出である。
