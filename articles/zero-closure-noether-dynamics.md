---
title: "離散零閉包からの Noether 保存則と関係位相動力学 ── 公理を壊さない自己写像としての動力学と標準模型一世代表現"
emoji: "🔁"
type: "idea"
topics: ["物理学", "幾何学", "ゲージ理論", "標準模型", "仮説"]
published: true
---

## 概要

「次元の生成構造」シリーズ、[閉包公理からの対称性導出](https://zenn.dev/noriaki_kihara/articles/closure-axioms-symmetry-derivation) の続編です（2026-08-21 公開、公開版 v1.0）。

- Version DOI: https://doi.org/10.5281/zenodo.22040736
- Concept DOI: https://doi.org/10.5281/zenodo.22040735
- Zenodo: https://zenodo.org/record/22040736

前論文では、複素零閉包 $\sum X_a^2=0$・有限回帰 $U^N=I$・simplex 閉包・自己無撞着という少数の条件から多数の対称構造を導きましたが、二つが残っていました。**Noether 型の保存則**と、**次状態を決める動力学**です。本論文はこの二つを扱い、さらに局所ゲージ幾何・標準模型一世代表現・hypercharge・anomaly cancellation・chirality 選択までの理論接続を与えます。

## 中心主張：状態の書換えを「動力学」として暗黙に許さない

通常の場の理論では、状態 $\Phi$ から $\Phi'$ への更新則をまず置き、それを時間発展と呼びます。本公理系ではそれができません。自己無撞着が基礎条件なので、**更新後の状態が零閉包などの許容条件から外れるなら、その更新則はこの理論の動力学ではない**からです。

問うべきは「次状態を作れるか」ではなく、

$$
\phi^{(n)}\in\mathcal Z_N\ \Longrightarrow\ \phi^{(n+1)}\in\mathcal Z_N
$$

を満たす**許容状態空間の内部の自己写像**を構成できるか、です。

## 関係位相から current と作用が出る

等振幅 $X_i=Ae^{i\phi_i}$ の下で、隣接関係の虚部と実部から

$$
J_{ij}=A^2\sin(\phi_j-\phi_i),\qquad
S_N=-A^2\sum_{\langle ij\rangle}\cos(\phi_j-\phi_i)
$$

が出ます。$\partial S_N/\partial\phi_i=-\sum_j J_{ij}$ なので、停留条件は離散 continuity equation $\sum_{j\sim i}J_{ij}=0$ です。これが有限 $N$ で**厳密に**成り立つ離散 Noether 保存則で、連続の $\partial_\mu J^\mu=0$ はその $N\to\infty$ 近似です。

## 零閉包を壊さない自己写像

零閉包 $\sum_i e^{2i\phi_i}=0$ は実2条件 $C_R=C_I=0$。未制約の力 $F_i=\sum_j\sin(\phi_j-\phi_i)$ をその接空間へ射影し、有限反復では法線方向の retraction で多様体上へ戻します：

$$
\phi^{(n+1)}=R_{\phi^{(n)}}\!\left(\eta\,P_{\phi^{(n)}}F(\phi^{(n)})\right),
\qquad
\mathcal F_N:\mathcal Z_N\to\mathcal Z_N .
$$

各反復で $\sum_i X_i^2=0$ が厳密保存されます。「公理を一度破って後から拘束し直す」のではなく、**公理を保つ写像として動力学が実現できる**──これが本論文の動力学に関する結論です。

## 自己写像パラメータ $s$ は物理時間ではない

ここは誤解されやすい点です。$s$ は自己無撞着解へ到達するための**構成・選択パラメータ**で、物理時間 $t$ ではありません。$t$ は前論文で導いた通り、観測不能な複素軸 $it$ の Lorentz 読出しとして現れます。自己写像の固定点は**時間軸を含む全配置**に対する条件で、それを Lorentz 読出しで読んだものが

$$
\Box\phi=0
$$

という場の方程式です。時間を基礎で特権化せずに、離散自己写像と連続場方程式が両立します。

## 局所ゲージ幾何

各頂点の位相原点を独立に選び直せる（無名性）とすると、頂点間の比較には辺に補償量 $\theta_{ij}$ が必要になります。これが $U(1)$ connection の離散ゲージ変換則を満たし、simplex の**面**の閉路和が gauge 不変な curvature になります。

$$
\text{頂点}\to\text{位相場},\quad
\text{辺}\to\text{connection},\quad
\text{面}\to\text{curvature}
$$

連続極限で $D_\mu=\partial_\mu-igA_\mu$、$F_{\mu\nu}$、Maxwell 型作用が現れ、頂点状態を $\mathbb C^r$ にすると閉路積から Yang–Mills の非可換項 $[A_\mu,A_\nu]$ が出ます。

## 標準模型一世代表現

前論文の複素5自由度 $V=V_3\oplus V_2$ と $S(U(3)\times U(2))$ に接続すると：

- $U(1)$ generator の trace-zero 条件 $3y_3+2y_2=0$ から **hypercharge 比** $-1/3:1/2$ が固定
- 一体側の双対 $V^*$ と、向き付き二体関係 $\Lambda^2V$ を合わせた $V^*\oplus\Lambda^2V$（$5+10=15$ 成分）が

$$
d^c\oplus L\oplus u^c\oplus Q\oplus e^c
$$

に分解──標準模型一世代の左手 Weyl 15 成分と一致（右手ニュートリノなし）

- $SU(3)^3$・$SU(3)^2U(1)$・$SU(2)^2U(1)$・$U(1)^3$・重力–$U(1)$ の全 anomaly と $SU(2)$ global anomaly の相殺を直接検算

$SU(5)$ は仮定していません。$\overline{\mathbf5}\oplus\mathbf{10}$ という既知の内容に、6軸零閉包 → 5自由度 → $3\oplus2$ の順で**到達した**形です。

## chirality 選択

残っていた共役二 Weyl sector の選択を、過去論文の A/B 二チャネル選択系と同定しました。mirror-odd な内部相関 $J=\operatorname{Im}(B^{*2}C)$ と既存の非線形選択項から、最小 normal form

$$
\dot S_\chi=\lambda J+gS_\chi(1-S_\chi^2)
$$

を得ます。$g>0$ で $S_\chi=0$ は不安定、$\pm1$ が安定。基礎方程式は左右対称のまま、$J$ の符号を種として一方が自発選択されます。

## 数値検証仕様と「公理保存監査」

後続の数値検証では、既存の A/B Fermi 型実験系の力学を変えず、観測量の追加と mirror run だけで上の normal form を fit します。最重要の判定は**選択が起きること**ではなく、

$$
\text{観測上の選択が成立}\ \land\ \text{公理保存が成立}
$$

の両方です。もっともらしい出力を出しても零閉包や有限回帰を壊す更新則は、この公理系の動力学として採用しません。

## 何が残っているか

本論文で「閉じた」と言うのは gauge group・局所 gauge geometry・一世代内部表現・hypercharge・anomaly cancellation・chirality 選択の鎖についてです。Higgs 動径モード、Yukawa、世代数、質量階層、量子補正は別の問題として残ります。

## 参照

- 前論文（閉包公理からの対称性導出 v1.0）: https://doi.org/10.5281/zenodo.22028072
- 論文一覧（GitHub Pages）: https://wurabeseiji.github.io/ai-chat-logs-open/
