# 倍音干渉による周期波内部の有限履歴読出し
## ― ベース位相を基準とする時間読出しゲージの理論解析 ―

**木原範昭**  
Version 1.0 / 2026-09-29

---

## 要旨

離散的な状態更新を考えるとき、次状態が直前の一状態だけから決まるとは限らない。放射、遅延相互作用、自己無撞着な場と運動の更新などを考えると、有限時間窓 \(\Delta T\) 内の複数の過去状態を区別して参照する必要が生じ得る。本稿では、その履歴を外部メモリに保存するのではなく、周期波そのものに保持させることが可能かを理論的に検討する。

単一周波数の波は、一周期内の位相を与えることはできるが、独立な時間構造を任意に形成する自由度に乏しい。これに対して、位相関係を保った複数の倍音を重ねると、Fourier 合成によって一周期内部により細かな時間構造を形成できる。ベース波 \(U=e^{i\omega_0\tau}\) に対する \(m\) 次倍音の相対位相は \(U^{m-1}\) となり、これを周期内部の位置を識別する「時間読出しゲージ」と解釈できる。

本稿は、Nyquist–Shannon の標本化理論、Fourier パルス合成、optical frequency comb、および delay-line memory に関する既知の結果を基礎として、倍音帯域の増加が周期内部の時間分解能を増加させることを整理する。その上で、一周期の位相構造を有限時間窓内の履歴セルとして読む物理解釈を提案する。この解釈は既知の信号理論そのものを変更するものではなく、波動系を「現在値」ではなく「有限長の循環履歴を保持する状態媒体」として扱うための写像を与えるものである。

---

## 1. 研究動機

これまでの離散状態更新の考え方では、現在の状態から次状態を生成するために、直前の離散時間幅 \(\Delta\tau\) に対応する状態を読み出せればよいと仮定できる。

しかし、重力波や電磁波を含む放射系、あるいは場と運動を反作用まで含めて自己無撞着に更新する系を考えると、次状態の決定に単一の \(\Delta\tau\) 前の状態だけでは不足し、有限時間窓

\[
\Delta T=N\Delta\tau
\]

に含まれる

\[
X(\tau-\Delta\tau),\quad
X(\tau-2\Delta\tau),\quad\ldots,\quad
X(\tau-N\Delta\tau)
\]

を区別して参照する必要がある可能性が生じる。

ここで本質的な問題は、これらの過去状態を**何が保持するのか**である。

外部の履歴配列や補助メモリを導入すれば計算上は容易である。しかし、閉じた物理系の状態生成として考えるなら、履歴を保持する物理的媒体自身を状態の中に求めたい。

そこで本稿では、状態保持媒体を周期波と仮定する。

単色波

\[
\psi_1(\tau)=A e^{i\omega_0\tau}
\]

だけでは、一周期内部を多数の独立な時間位置として識別する自由度が不足する。一方、複数の倍音を含む

\[
\Psi(\tau)=\sum_{m=1}^{M}A_m e^{im\omega_0\tau}
\]

では、倍音間の干渉によって一周期内部に細かな時間構造を形成できる。

したがって本稿の中心的な問いは、

\[
\boxed{
\text{倍音を持つ周期波は、有限時間窓内部の状態列を識別する媒体となり得るか}
}
\]

である。

この問題は信号理論、Fourier 合成、周波数コム、時間分解能、および遅延線記憶と密接に関係する。本稿では、それらの先行研究を基礎に実現可能性を理論解析する。

---

## 2. 単色波の限界

ベース周期を

\[
T_0=\frac{2\pi}{\omega_0}
\]

とする。

単色波

\[
\psi_1(\tau)=A e^{i\omega_0\tau}
\]

は、位相

\[
\phi=\omega_0\tau\pmod{2\pi}
\]

によって周期内部の位置を表すことができる。

しかし、その波形は振幅 \(A\) と位相 \(\phi\) が与えられれば決まり、一周期内部に任意の複数構造を形成する独立なスペクトル自由度を持たない。

したがって、単色波の一周期を

\[
\tau_0,\tau_1,\ldots,\tau_{N-1}
\]

という複数の履歴セルとして物理的に識別したい場合、位相を読み出すための追加構造が必要になる。

本稿では、その追加構造として倍音を考える。

![図1 単色波と倍音波の時間分解能比較](figure01_monochromatic_vs_harmonic_time_resolution.png)

**図1.** 単色波と複数倍音を含む波の周期内部構造の比較。倍音を加えることで、一つのベース周期内部により細かな干渉構造を形成できる。

---

## 3. 倍音による周期内部構造

有限個の整数倍音を

\[
\Psi_M(\tau)=\sum_{m=1}^{M}A_m e^{im\omega_0\tau}
\]

とする。

最も単純に \(A_m=1\) とすれば、位相 \(\phi=\omega_0\tau\) を用いて

\[
S_M(\phi)=\sum_{m=1}^{M}e^{im\phi}
\]

となる。

等比級数から

\[
S_M(\phi)
=e^{i(M+1)\phi/2}
\frac{\sin(M\phi/2)}{\sin(\phi/2)}
\]

を得る。

したがって振幅包絡は

\[
|S_M(\phi)|
=
\left|
\frac{\sin(M\phi/2)}{\sin(\phi/2)}
\right|.
\]

主ピークの代表的な幅は \(M\) の増加とともにおおむね

\[
\delta\phi\sim \frac{2\pi}{M}
\]

の尺度で狭くなる。

時間に戻せば

\[
\delta\tau
\sim
\frac{T_0}{M}.
\]

これは新しい標本化定理ではない。最高周波数が \(M\omega_0\) まで増えるため、利用可能な帯域が増加し、それに応じて短い時間構造を表現できるという通常の Fourier/time-bandwidth の結果である。

重要なのは、本稿ではこの既知の性質を**周期内部の状態読出し**として解釈する点である。

---

## 4. ベース波を基準とする時間読出しゲージ

ベース波を

\[
U(\tau)=e^{i\omega_0\tau}=e^{i\phi}
\]

とする。

\(m\) 次倍音は

\[
U^m=e^{im\phi}
\]

である。

ベース波との相対位相を積によって読み出せば、

\[
U^mU^{-1}=U^{m-1}=e^{i(m-1)\phi}
\]

となる。

したがって、絶対位相を外部時計に対して測るのではなく、ベース波自身を基準として各倍音の相対位相を読むことができる。

本稿では

\[
G_m(\phi)=U^{m-1}
\]

を **時間読出しゲージ (temporal readout gauge)** と呼ぶ。

\(m\) が増えるほど \(G_m\) は一つのベース周期内でより多く回転するため、周期内部のより細かな位相位置を区別できる。

この意味で、

\[
\boxed{
\text{倍音次数}\uparrow
\quad\Longrightarrow\quad
\text{周期内部の読出し解像度}\uparrow
}
\]

となる。

ただし、これは帯域を増やさずに情報量を無制限に増やす主張ではない。高次倍音を追加すること自体が帯域の増加であり、Nyquist–Shannon 理論と矛盾しない。

---

## 5. 周波数領域と時間領域

optical frequency comb では、等間隔の多数の周波数成分が位相関係を保って並び、時間領域では周期的なパルス列として表現される。Fortier and Baumann は、frequency comb を高安定基準から多数の optical tones へ位相・周波数情報を精密に移す optical synthesizer として整理している [3]。

また Hyodo らは、3本の連続波半導体レーザーを位相同期して Fourier 合成するだけでも、9.6 GHz の optical pulse train を実験的に生成できることを示した [4]。その後、3本の独立レーザーによる 1.8 THz pulse train の合成も報告されている [5]。

したがって、

\[
\text{離散周波数成分}
\longleftrightarrow
\text{周期的時間構造}
\]

という対応そのものは十分に確立された信号・光学上の事実である。

![図3 周波数コムと時間領域分解能](figure03_frequency_comb_to_time_domain_resolution.png)

**図2.** 周波数成分数を増やしたときの frequency-domain comb と time-domain waveform の対応。高次倍音まで位相整合して加えるほど、ベース周期内部の構造は細かくなる。

---

## 6. 一周期を有限履歴レジスタとして読む

ここからが本稿の物理解釈である。

ベース周期 \(T_0\) を有限時間窓

\[
\Delta T=T_0
\]

と対応させ、その内部を

\[
\Delta T=N\delta\tau
\]

と分割する。

周期内部の位相位置

\[
\phi_j=\frac{2\pi j}{N},
\qquad j=0,1,\ldots,N-1
\]

を、

\[
\phi_j
\longleftrightarrow
X(\tau-j\delta\tau)
\]

と対応させれば、一周期の波形を有限長の履歴レジスタとして読むことができる。

概念的には

\[
\boxed{
\text{memory cell}
\quad\longrightarrow\quad
\text{phase-resolved cell}
}
\]

という写像である。

これは、情報を静的な半導体セルに保持するのではなく、媒体を伝搬するパルスとして循環保持した初期の delay-line memory と構造的に似ている。

EDSAC などで用いられた mercury delay-line memory では、電気信号を音響パルスへ変換し、媒体中を伝搬させ、受信後に再増幅して循環させることで情報を保持した [6,7]。情報は空間的なランダムアクセスセルではなく、循環する波の時間順序として保持された。

本稿の提案は同一の装置原理ではないが、記憶概念としては近い。違いは、単なる伝搬遅延位置ではなく、**ベース波と倍音との干渉位相をアドレスとして用いる**点にある。

![図2 波周期FIFO](figure02_wave_period_fifo_phase_delay_memory.png)

**図3.** 通常の FIFO / delay-line memory と、周期波内部の位相セルとして履歴を読む概念の比較。波形自身が循環履歴媒体となる。

---

## 7. 離散読出しと循環シフト

一周期を \(N\) 点で読むとする。

\[
r_j=\Psi(\phi-j\delta\phi),
\qquad
\delta\phi=\frac{2\pi}{N}.
\]

一つの microstep だけ位相が進めば、理想化した循環表現では

\[
r_j' = r_{j-1}
\]

となる。

添字を \(N\) を法として扱えば

\[
r_0' = r_{N-1}
\]

であり、\(N\) 回のシフト後に

\[
R^N=I
\]

となる。

ここで \(R\) は循環シフト作用素である。

したがって有限周期波は、適切な位相分解能を持つなら、数学的には circular buffer / cyclic FIFO と同型の構造を持ち得る。

ただし重要なのは、\(N\) 個のセルを宣言するだけでは \(N\) 個の独立情報を保持したことにはならない点である。独立に復元できる自由度は、利用可能なスペクトル帯域、振幅・位相自由度、雑音、観測演算によって制限される。

---

## 8. 奇数倍音と偶数倍音の再解釈

倍音をベース波との相対位相として読むと、奇数倍音と偶数倍音の意味も変わる。

整数倍音では、

\[
m=2k+1
\]

の奇数倍音に対して

\[
U^{2k+1}U^{-1}=U^{2k}
\]

となり、ベース波との相対読出しは偶数次数となる。

一方、偶数倍音

\[
m=2k
\]

では

\[
U^{2k}U^{-1}=U^{2k-1}
\]

となる。

したがって「波そのものの倍音次数」と「ベース波に対する読出しゲージ次数」は一つずれる。

さらに二重被覆変数

\[
V=e^{i\phi/2},
\qquad
U=V^2=e^{i\phi}
\]

を導入すると、奇数半整数系列は

\[
V^{2k+1}=VU^k
\]

と分解できる。

偶数系列は

\[
V^{2k}=U^k
\]

である。

したがって奇数系列では、全モードに共通する \(V\) と、周期内部を細分する整数回転 \(U^k\) を分離して読むことができる。

\[
\boxed{
V^{2k+1}
=
\underbrace{V}_{4\pi\text{ carrier}}
\underbrace{U^k}_{\text{intra-period readout gauge}}
}
\]

この分解は、4\(\pi\) 周期を持つ位相構造と、2\(\pi\) 周期内部の時間読出し構造を混同しないために有用である。

![図4 奇数偶数倍音の役割](figure04_odd_even_harmonic_roles_and_readout_gauge.png)

**図4.** \(V=e^{i\phi/2}\), \(U=V^2\) による奇数・偶数系列の分解。奇数系列は共通の 4\(\pi\) carrier \(V\) と整数位相ゲージ \(U^k\) に分けられる。

なお、この 4\(\pi\) 周期性だけから量子力学的 spin-1/2 の完全な変換則を導いたとは言えない。本稿では位相構造上の二重被覆としてのみ扱う。

---

## 9. 情報量と分解能についての注意

本稿の提案は、単一の有限帯域に Shannon/Nyquist 限界を超える任意量の情報を埋め込めるという主張ではない。

Nyquist は有限帯域通信における信号速度とパルス間干渉の条件を体系化し [1]、Shannon は帯域制限された信号の標本化・通信に関する基本構造を示した [2]。

本稿では高次倍音を追加するため、最高周波数

\[
\omega_{\max}=M\omega_0
\]

も増加する。

したがって

\[
M\uparrow
\quad\Rightarrow\quad
B\uparrow
\quad\Rightarrow\quad
\delta\tau\downarrow
\]

という関係は、通常の time-bandwidth tradeoff と整合する。

本稿の新しい点は情報容量の上限を変更することではなく、**増加した時間分解能を、有限履歴の物理的読出しに対応させること**である。

---

## 10. 現代の遅延系との接点

波動・遅延系を履歴依存計算に利用する考え方は、現代の physical reservoir computing にも現れている。

特に delay-based photonic reservoir computing では、物理的な遅延と時間多重化された virtual nodes を用いて時系列情報を処理する。2026年には thin-film lithium niobate 上で time–wavelength-coupled virtual nodes を用いる delay-based photonic reservoir が報告されている [9]。

この種の研究は、本稿の「一つの物理波動系の時間構造を複数の有効状態として読む」という発想に近い。ただし reservoir computing は入力履歴を計算資源として利用することを目的とするのに対し、本稿では、**物理系自身の過去状態を同じ波の内部位相構造として保持できるか**を問題にしている。

したがって両者は技術目的が異なるが、物理的な時間多重化・遅延・波動記憶という観点で重要な先行研究となる。

---

## 11. 遅延相互作用への物理解釈

有限履歴レジスタとしての波を導入する動機は、遅延を含む物理相互作用にある。

一般に、有限伝播速度を持つ相互作用では、現在の相互作用が単一の現在値だけでなく、過去の source state と関係する。

通常の計算では、

\[
\{X(\tau-j\delta\tau)\}_{j=1}^{N}
\]

を履歴配列として外部に保存できる。

しかし、閉じた状態生成系として考えるなら、

\[
\boxed{
\text{history array}
\longrightarrow
\text{phase-resolved wave state}
}
\]

と写像できれば、履歴を外部メモリとして追加する必要がなくなる可能性がある。

このとき一つの波は単なる瞬間値

\[
\psi(\tau)
\]

ではなく、

\[
\boxed{
\text{有限時間窓を内部位相構造として保持する状態媒体}
}
\]

として扱われる。

これは、次状態を作るための必要情報が現在の物理状態の外部に隠れていてはならない、という閉じた状態更新の要求に対して一つの候補を与える。

ただし本稿では、具体的な Maxwell 方程式または Einstein 方程式にこの写像を組み込み、自己無撞着な放射反作用解を構成するところまでは扱わない。

---

## 12. 成立に必要な条件

この解釈が物理的な履歴保持として成立するためには、少なくとも次の条件が必要である。

### 12.1 位相コヒーレンス

各倍音の相対位相が履歴読出しに十分な期間保持されなければならない。

\[
\phi_m-m\phi_1
\]

が無制御に拡散すれば、時間アドレスとしての干渉構造は失われる。

### 12.2 十分な帯域

必要な時間分解能 \(\delta\tau\) に対して、十分な最高倍音次数 \(M\) が必要である。

概略的には

\[
\delta\tau\gtrsim \frac{T_0}{M}
\]

が一つの尺度となる。

### 12.3 読出し可能性

波形に構造が存在するだけでは不十分であり、その構造を系内部の相互作用によって区別して読み出せなければならない。

したがって必要なのは

\[
\text{storage}
+
\text{addressing}
+
\text{readout}
\]

の三つである。

本稿で倍音干渉が与えるのは主として addressing/readout gauge の候補である。

### 12.4 更新後の履歴シフト

次状態生成後には、履歴が

\[
X_j\rightarrow X_{j+1}
\]

に相当する形で波内部に更新される必要がある。

したがって完全な物理メモリモデルには、波形生成だけでなく、履歴を書き込みながら循環させる動力学が必要である。

---

## 13. 本稿で示したことと示していないこと

本稿で示したことは以下である。

1. 複数の位相整合した倍音が、一つのベース周期内部に単色波より細かな時間構造を形成できること。
2. ベース波との相対位相 \(U^{m-1}\) を周期内部の読出しゲージとして定義できること。
3. 利用帯域の増加に伴って周期内部の時間分解能が向上することが、既知の Fourier・sampling 理論と整合すること。
4. 一周期を有限時間窓に対応させれば、位相分解された波を circular buffer / delay-line memory に類似した有限履歴媒体として解釈できること。
5. \(V^{2k+1}=VU^k\) により、4\(\pi\) carrier と周期内部の整数位相読出しを分離できること。

一方、本稿では次のことを示していない。

1. 自然界のすべての波が実際に過去状態をこの方法で保存していること。
2. \(N\) 個の位相セルが無条件に \(N\) 個の独立情報自由度を持つこと。
3. 重力波または電磁波の完全な自己無撞着遅延方程式がこの表現だけで解けること。
4. 4\(\pi\) carrier が量子力学的 spin-1/2 の完全な導出になっていること。
5. Shannon/Nyquist の情報容量または標本化限界を超えること。

この区別は、本稿の提案を信号理論上の既知事実と物理解釈に明確に分離するために重要である。

---

## 14. 結論

本稿では、有限時間窓の履歴を外部配列として保持する代わりに、倍音を持つ周期波そのものに保持させる可能性を検討した。

位相整合した倍音群

\[
\Psi(\tau)=\sum_{m=1}^{M}A_m e^{im\omega_0\tau}
\]

は、一つのベース周期内部に細かな時間構造を形成する。ベース波

\[
U=e^{i\omega_0\tau}
\]

に対する各倍音の相対位相

\[
G_m=U^{m-1}
\]

を時間読出しゲージとすれば、倍音次数の増加を周期内部の読出し解像度の増加として解釈できる。

この構造は Fourier synthesis、optical frequency comb、sampling theory と整合し、記憶概念としては古典的 delay-line memory や現代の delay-based physical computing と接点を持つ。

本稿で提案した追加の物理解釈は、

\[
\boxed{
\text{一周期の位相構造}
\quad\leftrightarrow\quad
\text{有限時間窓の履歴}
}
\]

という対応である。

この対応が成立するなら、波は単なる瞬間的な状態量ではなく、有限長の過去を内部に保持する循環状態媒体として扱える。

次の課題は、この読出し構造を具体的な遅延相互作用に組み込み、履歴の書込み・循環・読出しを同一の閉じた状態更新則から生成できるかを調べることである。

---

## 参考文献

[1] H. Nyquist, “Certain Topics in Telegraph Transmission Theory,” *Transactions of the American Institute of Electrical Engineers*, vol. 47, no. 2, pp. 617–644, 1928. DOI: 10.1109/T-AIEE.1928.5055024.

[2] C. E. Shannon, “Communication in the Presence of Noise,” *Proceedings of the IRE*, vol. 37, no. 1, pp. 10–21, 1949. DOI: 10.1109/JRPROC.1949.232969.

[3] T. Fortier and E. Baumann, “20 years of developments in optical frequency comb technology and applications,” *Communications Physics*, vol. 2, Article 153, 2019. DOI: 10.1038/s42005-019-0249-y.

[4] M. Hyodo, N. Onodera, and K. S. Abedin, “Fourier synthesis of 9.6-GHz optical-pulse trains by phase locking of three continuous-wave semiconductor lasers,” *Optics Letters*, vol. 24, no. 5, pp. 303–305, 1999. DOI: 10.1364/OL.24.000303.

[5] M. Hyodo, K. S. Abedin, and N. Onodera, “Fourier synthesis of 1.8-THz optical-pulse trains by phase locking of three independent semiconductor lasers,” *Optics Letters*, vol. 26, no. 6, pp. 340–342, 2001. DOI: 10.1364/OL.26.000340.

[6] J. P. Eckert Jr., “A Survey of Digital Computer Memory Systems,” *Proceedings of the IRE*, 1953. Historical context and delay-line implementation are also documented by the Computer History Museum, “1949: EDSAC computer employs delay-line storage.”

[7] Computer History Museum, “Delay Lines,” *Revolution: The First 2000 Years of Computing*. Mercury and magnetostrictive delay lines stored serial information as circulating acoustic pulses; accessed 2026-09-29.

[8] L. Chang, S. Liu, and J. E. Bowers, “Integrated optical frequency comb technologies,” *Nature Photonics*, vol. 16, pp. 95–108, 2022. DOI: 10.1038/s41566-021-00945-1.

[9] D. Kong et al., “Delay-based photonic reservoir computing on thin-film lithium niobate with time–wavelength-coupled virtual nodes,” *npj Unconventional Computing*, vol. 3, Article 36, 2026. DOI: 10.1038/s44335-026-00081-5.

---

## 図の再現性

本稿の図1–4は外部画像の転載ではなく、解析式および概念構造から独自に作成した。再生成用プログラム、SVG、PNG、および説明ファイルは本稿と同じ研究系列の

`倍音干渉による有限履歴読出し_説明図_v1/`

に保存した。

再生成プログラム：

`generate_harmonic_history_figures_v1.py`

図ファイル：

- `figure01_monochromatic_vs_harmonic_time_resolution.svg/.png`
- `figure02_wave_period_fifo_phase_delay_memory.svg/.png`
- `figure03_frequency_comb_to_time_domain_resolution.svg/.png`
- `figure04_odd_even_harmonic_roles_and_readout_gauge.svg/.png`

