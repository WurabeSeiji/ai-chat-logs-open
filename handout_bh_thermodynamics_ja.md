# BH熱力学プログラム：6篇の研究プログラム概要

**著者：** 木原 範昭（Noriaki Kihara）  
**所属：** WF System Co., Ltd. / 大阪大学 基礎工学部 卒業  
**ORCID：** [0009-0004-6753-4020](https://orcid.org/0009-0004-6753-4020)

---

## 3つの主張

**主張 1：4次元球への単位立方体充填の漸近係数 $c = 8/(3\pi)$**  
半径 $R = 2k+1$ の4次元球内に整数中心の単位立方体を完全包含で詰めた個数 $N(k)$ について、体積不足 $\Delta(R) = V_4(R) - N(k)$ の漸近係数 $c$ を包除原理から閉形式で導出。$c = 8/(3\pi) \approx 0.84883$。（論文3）

**主張 2：4次元 Tangherlini ブラックホール熱力学との完全一致**  
ブラックホールのエントロピーは事象の地平面の表面積に比例するので、直径 $R$ の4次元球では $R^3$ に比例する。中心投影由来計量の自己整合的決定から $R \propto M^{1/2}$、Hawking 温度 $T_H \propto M^{-1/2}$、寿命 $\tau \propto M^{3/2}$、蒸発レート $dM/dt \propto -M^{-1/2}$ を導出。Bekenstein–Hawking エントロピーとの比 $\Delta(R)/S_{BH} \to 32/(3\pi) \approx 3.40$。（論文4・5）

**主張 3：ホログラフィー非仮定の独立な経路**  
't Hooft–Susskind ホログラフィー原理、Maldacena AdS/CFT、特定の量子重力プログラム（LQG/CDT/弦理論）を**仮定しない**。中心投影と離散性という最小限の幾何学的・組合せ論的仮定のみから上記結果を導出する独立経路。自己引用ゼロ。（論文1〜6）

---

## 研究の全体構成

**論文1：俯瞰・仮説提示**  
仮定 A（BH エントロピーは幾何学的）・B（蒸発は地平面局在）・C（時空は Planck スケールで離散）から、仮説1（射影仮説）と仮説2（離散ドリフト仮説）を抽出。

**論文2：中心投影の数学的基盤**  
$S^4(R) \subset \mathbb{R}^5$ から4次元接超平面へのグノモン投影。誘導計量 $g_{\mu\nu} = (R^2/\ell^2)(\delta_{\mu\nu} - x_\mu x_\nu/\ell^2)$ は $G_{\mu\nu} + (3/R^2)g_{\mu\nu} = 0$ を満たし、Lorentz 連続化で de Sitter、質量導入で Schwarzschild–de Sitter 形。

**論文3：離散構造の数学**  
充填条件 $\sum (|c_i|+1/2)^2 \le R^2$、$N(k)$ を $k = 0..60$ で正確計算（$N(60) = 1{,}028{,}515{,}513$）。$\Delta(R) = (16\pi/3)R^3 - 6\pi R^2 + O(R)$。Lagrange–Jacobi 四平方表現数 $r_4(N) = 8\sigma(N)$ と接続。

**論文4・5：物理接続と動力学**  
$\Delta(R)$ と4次元 Bekenstein–Hawking エントロピーの定数比 $32/(3\pi)$。Friedmann 類似 $E_{\mathrm{total}} \propto R^2$ から Tangherlini スケーリング全項目を回復。

**論文6：統合と未解決問題**  
4次元での結果が観測 3+1 BH 物理学にどう削減されるかを主要未解決問題として精密に表明。3候補機構（KK・Randall–Sundrum/ADD・dS-CFT）。

---

## 6篇のDOI一覧

| 論文 | DOI |
|------|-----|
| 1（俯瞰）  | [10.5281/zenodo.19837587](https://doi.org/10.5281/zenodo.19837587) |
| 2（射影） | [10.5281/zenodo.19837589](https://doi.org/10.5281/zenodo.19837589) |
| 3（充填） | [10.5281/zenodo.19837591](https://doi.org/10.5281/zenodo.19837591) |
| 4（エントロピー） | [10.5281/zenodo.19837593](https://doi.org/10.5281/zenodo.19837593) |
| 5（動力学） | [10.5281/zenodo.19837595](https://doi.org/10.5281/zenodo.19837595) |
| 6（統合） | [10.5281/zenodo.19837597](https://doi.org/10.5281/zenodo.19837597) |

## 4次元 Tangherlini との対応

| 量 | 中心投影構成 | 4次元 Tangherlini |
|------|--------------|-----------------|
| 質量–半径 | $R \propto M^{1/2}$ | $r_h \propto M^{1/2}$ |
| Hawking 温度 | $T_H \propto M^{-1/2}$ | 同上 |
| 寿命 | $\tau \propto M^{3/2}$ | 同上 |
| エントロピー | $\Delta(R) \propto R^3$ | $S_{BH} \propto R^3$ |
| 比 | $\Delta/S_{BH} \to 32/(3\pi)$ | – |

---

## 本研究が主張しないこと

- 物理的時空が 4 次元であること
- 高次元時空の物理的実在
- Bekenstein–Hawking エントロピーの第一原理導出
- ホログラフィック・量子重力プログラムの代替

蒸発レート指数の「不一致」は 3+1 と 4 規約の混同であり、4 次元では完全一致。真の未解決問題は 4 → 3+1 削減機構の同定。

---

**Zenn 記事**: https://zenn.dev/noriaki_kihara/articles/bh-thermodynamics-projection  
**note 記事（日本語）**: https://note.com/kiharanoriaki/n/n4524006ef175  
**note 記事（英語）**: https://note.com/kiharanoriaki/n/n916eac778c1c  
**ソースコード（GitHub）**: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/BH熱力学プログラム  
全論文のPDF・LaTeX・MarkdownはZenodoに公開済み（日英対訳、CC BY 4.0）
