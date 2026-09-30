# EM 自己無撞着有限履歴モデル v1 — 設計と初回検証

## 目的
現行モデルの leading electric-dipole 放射式を使わず、古典点電荷の標準理論を直接使う。

## 採用した物理
1. 相互場: 任意運動点電荷の retarded Lienard-Wiechert field。
2. 運動: Lorentz force。
3. 自己力: covariant Lorentz-Dirac (LAD) equation。
4. 履歴: retarded time を解くための有限 worldline window。
5. 波内部履歴との接続: 有限窓全体を Fourier/harmonic coefficients に可逆変換し、全モード保持時に履歴を損失なく再構成できることを検査する。

## 主方程式
計量符号を (-,+,+,+)、c=1, 4 pi epsilon0=1 とする。

m a^mu = f_ext^mu + (2/3) q^2 (d a^mu/dtau - a^2 u^mu)

これを

du^mu/dtau = a^mu

da^mu/dtau = (a^mu - f_ext^mu/m)/tau0 + a^2 u^mu

tau0 = 2 q^2/(3m)

として積分する。

相互電磁場は Lienard-Wiechert の velocity term + acceleration/radiation term を retarded time で評価する。

## 重要な点
- dipole radiation power を後付けしない。
- aa / ab を別の放射則として与えない。
- 相手からの retarded field と自己の LAD self-force を同じ worldline history から評価する。
- LL reduction は主系に使用していない。

## 数値条件
LAD は runaway / preacceleration を持ち得るため、単純な大刻み積分は不安定になる。初回試行では dt が tau0 より大きく失敗した。そこで

q1=+0.1, q2=-0.1, m1=m2=1,
separation=1,
dt=1e-4,
history=12000

とし、history window > light travel time を満たした。

## 初回検証 A: radial control
50 steps。u.u=-1 および u.a=0 を数値精度で維持。全履歴の Fourier encode/decode 最大誤差は約 1.1e-15。

## 初回検証 B: transverse motion
v1=+0.05 y, v2=-0.05 y、200 steps。retarded acceleration field が非ゼロになることを確認した。有限履歴の Fourier encode/decode は機械精度で可逆。

## 解釈
これは「解析的閉形式解」を得たという意味ではない。Maxwell の retarded point-charge solution と Lorentz-Dirac point-charge equation をそのまま数値的に結合した参照実装である。

次の課題は、現在配列として保持している finite worldline window を、先に提案した倍音・位相セル表現へ置き換え、同じ retarded field / LAD result を再現できるか確認することである。
