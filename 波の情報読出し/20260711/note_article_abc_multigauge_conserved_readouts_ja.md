# 質量や運動量やエネルギーは、どこから読まれているのか

新しい論文を公開しました。

今回のテーマは、質量、運動量、エネルギーです。

と言っても、標準物理で使われている質量、運動量、エネルギーそのものを、いきなり導出したという話ではありません。

今回調べたのは、その一歩手前です。

背景座標を先に置かない閉じた位相系の中で、質量的、運動量的、エネルギー的に見える保存量を、干渉から読み出せるのか。

この問いを、ABC 閉鎖位相系という最小モデルで数値的に検証しました。

正式版の論文、英訳、TeX/PDF、実行スクリプト、結果データは Zenodo に公開しています。

・Concept DOI（常に最新版）: https://doi.org/10.5281/zenodo.21308049

・この版: https://doi.org/10.5281/zenodo.21308050

・Zenn 記事: https://zenn.dev/noriaki_kihara/articles/abc-multigauge-conserved-readouts

## 背景座標を先に置かない

これまでの波の情報読出しシリーズでは、空間や時間を、最初から背景として置くのではなく、閉じた位相系からの読出しとして扱ってきました。

つまり、x, y, z, t が先にあるのではなく、位相差、干渉、参照波、読出し窓の関係から、空間的・時間的な量が読まれる、という立場です。

今回の問いは、その次です。

空間や時間が読出しなら、質量、運動量、エネルギーも、読出しとして扱えるのではないか。

この問いは、かなり危ない問いでもあります。

なぜなら、質量、運動量、エネルギーは、物理学の中では非常に基本的な量だからです。

だからこそ、今回の論文では、標準物理量そのものを導出した、とは言っていません。

まずは、閉鎖位相系の内部で、それらに見える保存読出しが構成できるかを調べました。

## ABC 閉鎖位相系

今回使ったのは、局所波 A、局所波 B、観測機 C からなる ABC 系です。

A と B は、衝突する局所波です。

C は外部観測者ではありません。

同じ閉鎖位相系の中にある参照波であり、A や B との干渉相関を読むための観測機です。

ここで重要なのは、単一のゲージで得られた値を、そのまま測定値とはみなさないことです。

一つの見方だけで得た値は、まだ読出しとして不十分です。

複数の読出しゲージ、複数の参照波、複数の読出し窓を通して、同じ量が安定に再構成されることを要求します。

これを、多ゲージ干渉読出しと呼んでいます。

## 何を読んだのか

今回の実験では、三つの読出し量を使いました。

空間位相方向の相関勾配を、

p_read

として読みます。

これは運動量そのものではありません。

運動量的に振る舞う、空間位相勾配の読出しです。

時間位相方向の相関勾配を、

E_read

として読みます。

これもエネルギーそのものではありません。

エネルギー的に振る舞う、時間位相勾配の読出しです。

そして、複数ゲージを変えても安定して残る振幅二乗残差を、

R_read

として読みます。

これは質量そのものではありません。

質量的に振る舞う、安定残差の読出しです。

## R は測りにくい

今回、かなり面白かったのは R の扱いです。

空間位相勾配 p は、反転や相対勾配の変化として比較的読みやすい。

時間位相勾配 E も、時間方向の読出し窓を変えれば追いやすい。

しかし R は違います。

R は、変動が小さいからこそ安定している量です。

でも、変動が小さいということは、測りにくいということでもあります。

つまり、質量的なものは、そもそも読みにくい。

これは直感的にも面白い点です。

大きく変動している軸を t と読み、安定して残っている軸を R と読む。

そう考えると、t と R の区別も、最初から固定された名前ではなく、読出しの安定性から決まっていることになります。

今回の実験では、この t/R 分離も、多ゲージ読出しの中で確認しました。

## 対称衝突では何が起きたか

まず、等しい振幅を持つ A と B の単回衝突を調べました。

結果として、p_read、E_read、R_read は、複数ゲージから非常に高い精度で再構成されました。

最大誤差は、

・p_read: 2.5202062658991053e-14

・E_read: 2.2315482794965646e-14

・R_read: 4.440892098500626e-16

でした。

この条件では、p は反転し、E と R は保存されました。

さらに、8回の反復衝突でも、p 反転、E と R の保存、識別振動の保存、補償付き閉鎖が維持されました。

つまり、単回の偶然ではなく、反復しても読出し構造が壊れないことを確認しました。

## 非対称 R では、単純反転が壊れる

次に、A と B の R が異なる場合を調べました。

等振幅条件では、単純に

q → -q

と反転させれば、反射のように見えます。

しかし、R_A と R_B が異なる場合、この単純反転は保存量を壊します。

特に、

R_A p_A + R_B p_B

が保存されません。

これはモデルの失敗ではありません。

むしろ、R が質量的な量として読まれるなら、保存写像も R で重み付けされた形に一般化される必要がある、という診断です。

この点が今回の重要な分岐でした。

## R 重み付きの保存写像

そこで、R_A p_A + R_B p_B と、

R_A p_A^2 + R_B p_B^2

を保存する一般化衝突写像を構成しました。

これは、標準物理の運動量保存やエネルギー保存を、そのまま仮定したという意味ではありません。

閉鎖位相系の中で読まれた R と p に対して、R*p と R*p^2 が保存されるかを調べたものです。

8種類の非対称振幅ケースでは、

・R*p 保存誤差: 2.3803181647963356e-13

・R*p^2 保存誤差: 1.4086509736443986e-12

の精度で保存が確認されました。

さらに、非単位・非対称の位相勾配、同方向の追いつき衝突、複数回衝突、読出しノイズ、極端な R 比も調べました。

R_B/R_A は、0.015625 から 64.0 まで掃引しました。

統合サマリーでは、9本の実験すべてが valid になりました。

しかも、単一ゲージだけで判定した実験は一つもありません。

## 何が見えてきたのか

今回の結果から見えてきたのは、質量、運動量、エネルギーを、最初から外部にある実体量として置かなくてもよい可能性です。

少なくとも、この ABC 閉鎖位相系の中では、

・空間位相勾配が、運動量的な読出しになる

・時間位相勾配が、エネルギー的な読出しになる

・安定振幅二乗残差が、質量的な読出しになる

・R*p が、運動量保存のように振る舞う

・R*p^2 が、二乗量保存のように振る舞う

という構造を、多ゲージ干渉から一貫して構成できました。

これは、かなり重要な結果だと思っています。

なぜなら、これまで空間や時間を位相の読出しとして扱ってきた流れが、質量的・運動量的・エネルギー的な保存量にも拡張できる可能性が出てきたからです。

## これは何を主張していないか

誤解を避けるために、はっきり書いておきます。

今回の論文は、標準物理の質量、運動量、エネルギーを完全に導出したものではありません。

標準力学を再導出したものでもありません。

実在粒子の衝突を定量予言するものでもありません。

今回確認したのは、ABC 閉鎖位相系の内部で、質量的・運動量的・エネルギー的に見える保存読出しが、多ゲージ干渉から構成できる、ということです。

標準物理量との対応写像は、次の課題です。

ただし、その対応写像を作るための読出し側の土台は、今回かなり明確になりました。

## まとめ

今回の数値実験では、背景座標を先に置かない ABC 閉鎖位相系から、p_read、E_read、R_read を多ゲージ干渉により再構成しました。

対称衝突では、p が反転し、E と R が保存されました。

非対称 R 条件では、単純な q → -q 反転が R*p 保存を破ることを確認しました。

そのうえで、R*p と R*p^2 を保存する一般化衝突写像を構成し、非対称振幅、非対称速度、複数回衝突、ノイズ、極端 R 比に対しても保存読出しが維持されることを確認しました。

質量、運動量、エネルギーを、先験的な実体量として置くのではなく、閉鎖位相系からの保存読出しとして扱う。

今回の論文は、そのための最初の数値構成実験です。

――――

正式版の論文・英訳・TeX/PDF・実行スクリプト・結果データは Zenodo に公開しています。

・Concept DOI: https://doi.org/10.5281/zenodo.21308049

・Version DOI: https://doi.org/10.5281/zenodo.21308050

・Zenn 記事: https://zenn.dev/noriaki_kihara/articles/abc-multigauge-conserved-readouts

<!-- pdf-links -->
論文本体の PDF は、公開リポジトリから直接ダウンロードいただけます。

- ab_two_body_harmonic_readout_c1_area_sweep_preliminary_summary_en.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/ab_two_body_harmonic_readout_c1_area_sweep_preliminary_summary_en.pdf
- ab_two_body_harmonic_readout_c1_area_sweep_preliminary_summary_ja.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/ab_two_body_harmonic_readout_c1_area_sweep_preliminary_summary_ja.pdf
- abc_closed_phase_independent_c_gauge_relation_decomposition_distance_exponent_preliminary_summary_en.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/abc_closed_phase_independent_c_gauge_relation_decomposition_distance_exponent_preliminary_summary_en.pdf
- abc_closed_phase_independent_c_gauge_relation_decomposition_distance_exponent_preliminary_summary_ja.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/abc_closed_phase_independent_c_gauge_relation_decomposition_distance_exponent_preliminary_summary_ja.pdf
- abc_closed_phase_system_multigauge_conserved_readouts_en.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/abc_closed_phase_system_multigauge_conserved_readouts_en.pdf
- abc_closed_phase_system_multigauge_conserved_readouts_ja.pdf
  https://raw.githubusercontent.com/WurabeSeiji/ai-chat-logs-open/main/波の情報読出し/20260711/abc_closed_phase_system_multigauge_conserved_readouts_ja.pdf

リポジトリ: https://github.com/WurabeSeiji/ai-chat-logs-open

#物理学 #量子力学 #複素数 #波動 #干渉 #保存量 #質量 #運動量 #エネルギー #数値実験 #シミュレーション #独立研究 #Zenodo #サイエンス
