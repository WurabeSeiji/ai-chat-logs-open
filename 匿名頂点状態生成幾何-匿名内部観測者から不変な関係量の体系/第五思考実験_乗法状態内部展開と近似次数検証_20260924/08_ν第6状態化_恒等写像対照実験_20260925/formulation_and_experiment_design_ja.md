# νの第6状態化：定式化と対照実験設計

作成日: 2026-09-25

## 1. 目的

この実験の目的は、新しい物理を導入することではない。第五思考実験で状態外の関数引数として残っていた対称質量比

$$
\nu=\frac{m_A m_B}{(m_A+m_B)^2}=\frac{q}{(1+q)^2}
$$

を、第6の永続状態として内部化し、旧5状態系の出力を変更せずに再現できるかを検証する。

旧系は概念的に

$$
Z_5' = F(Z_5;\nu),\qquad Z_5=(U,P,E,H,Q)
$$

であり、\(\nu\) が状態ベクトル外に存在していた。

新系は

$$
Z_6=(U,P,E,H,Q,\mathcal N)
$$

とし、

$$
Z_6'=\widetilde F(Z_6)
$$

だけで更新する。新系の `micro_six_fast(s)` と `macro_six_fast(s)` は状態 `s` 以外の物理引数を受け取らない。

この対照実験では電荷自由度はまだ導入しない。したがって、今回確認するのは「重力-only Paper 5 における ν の状態化」だけである。

## 2. \(\mathcal N\) の乗法状態表現

\(\nu>0\) なので、正の乗法群要素として

$$
\boxed{\mathcal N_0=e^{\log\nu_0}=\nu_0}
$$

と置ける。

ただし数値実装では `exp(log(nu))` を実際に往復計算しない。旧実装の浮動小数点値を厳密に保持するため、

```text
Nscale := nu
```

とする。これは数学的には上式と同じ表現だが、IEEE-754 のビット列を保存するための実装上の厳密化である。

実際、負の対照として \(q=4\) を調べると

$$
\nu=0.16,
$$

に対し、数値的な

$$
\exp(\log\nu)=0.15999999999999998
$$

となり、1 ULP だけ値が変わる。今回の力学核では12000 step後も軌道変数には差が現れなかったが、状態そのもののビット一致は失われるため、厳密な対照実験では採用しない。

## 3. 恒等写像を同一相互作用の積和として実装

第五思考実験のB方式では11相の one-hot 位相状態

$$
q_m=(q_{0,m},\ldots,q_{10,m}),\qquad
q_{j,m}\in\{0,1\},\qquad
\sum_{j=0}^{10}q_{j,m}=1
$$

を使う。

\(\mathcal N\) の各相における自己係数を

$$
\boxed{c^{(j)}_{\mathcal N\leftarrow\mathcal N}=1\quad(j=0,\ldots,10)}
$$

とし、旧5状態から \(\mathcal N\) への交差係数を

$$
\boxed{c_{\mathcal N\leftarrow U}
=c_{\mathcal N\leftarrow P}
=c_{\mathcal N\leftarrow E}
=c_{\mathcal N\leftarrow H}
=c_{\mathcal N\leftarrow Q}=0}
$$

とする。

すると1 microstep の \(\mathcal N\) 更新は

$$
\boxed{
\mathcal N_{m+1}
=
\sum_{j=0}^{10}
q_{j,m}\,1\,\mathcal N_m
=
\mathcal N_m
}
$$

となる。

コードでは

```python
n['Nscale'] = sum(
    float(o['q'][k]) * (IDENTITY_N_COEFF[k] * nscale)
    for k in range(11)
)
```

とし、`Nscale = old_Nscale` の別代入を恒等更新として用いていない。

## 4. 旧 \(\nu\) 依存式の置換

旧放射反作用核は

$$
C_{\rm RR}=1.6\,\nu\,r^{-3}
=\frac85\nu r^{-3}
$$

を使う。

新系では物理式・演算順序を変えず、

$$
\boxed{
C_{\rm RR}=1.6\,\mathcal N\,r^{-3}
}
$$

とする。

つまり変更点は「数値」ではなく「情報の所在」だけである。

旧：

```text
micro_legacy(s, nu_external)
rr_kernel(U,P,E,nu_external)
```

新：

```text
micro_six_fast(s)
Nscale = s['Nscale']
rr_kernel(U,P,E,Nscale)
```

新 generator の macro/micro 更新には外部 `nu` 引数を残さない。

## 5. 読み出し

Paper 5 の軌道読み出しは

$$
\cos\chi=\frac12(U+U^{-1}),
$$

$$
r=\frac{P}{1+E\cos\chi},
$$

$$
x=r\operatorname{Re}H,\qquad
y=r\operatorname{Im}H,
$$

$$
t=\tau_0\log Q
$$

であり、\(\nu\) を使用しない。

したがって新系でも読み出し式は変更しない。必要な監査時のみ

$$
\nu_{\rm audit}=\mathcal N
$$

と読むことができるが、この値は状態更新へ戻さない。

## 6. 対照実験条件

### 6.1 12ケース

過去の厳密監査と同じケース集合を用いる。

- \(p_0\in\{24.3,60,120\}\)
- \(e_0\in\{0.20,0.45,0.55\}\)
- 上記9ケースは \(q=1\)
- 追加3ケースは \((p_0,e_0)=(60,0.45)\), \(q\in\{2,4,10\}\)
- 4000 macro step
- 1 macro step = 11 synchronous microsteps

比較対象は

$$
(U,P,E,H,Q),\quad (x,y,r,t),\quad\mathcal N
$$

である。

### 6.2 12000-step 長時間対照

次の3ケースを12000 macro stepまで比較する。

1. \((p_0,e_0,q)=(24.3,0.4358898943540673,1)\)
2. \((60,0.45,4)\)
3. \((120,0.55,1)\)

### 6.3 ランダム microstep 監査

seed=20260924 で500初期状態を生成し、各11 microsteps、合計5500 microstepsについて旧系と新系をビット比較する。

### 6.4 保存済み Paper-4 C1 データとの照合

保存済み

`experiment02_sn_2p5pn_rr_C1_N12.csv`

を第三の参照として用いる。旧実装と新6状態実装の両方を同じ初期条件から走らせ、保存済みCSVとの差と、旧対新の差を同時に記録する。

## 7. 合格条件

今回の変更は情報の所在だけを変えるため、単なる許容誤差ではなく、旧対新について

$$
\boxed{(U,P,E,H,Q)_{\rm new}=(U,P,E,H,Q)_{\rm old}}
$$

を浮動小数点ビット列まで一致させることを目標とする。

さらに

$$
\boxed{\mathcal N_n=\mathcal N_0=\nu_0}
$$

が全 step で成立することを要求する。

## 8. この実験で証明しないこと

この対照実験が成功しても、次は証明されない。

- 電荷を追加した Einstein-Maxwell 近似反応系が同じ6状態だけで閉じること。
- 電荷状態や結合スケールが恒等状態でよいこと。
- 既存B方式内部に残る `P+hK` 型加法中間演算が、より強い乗法閉包条件を満たすこと。
- 現在の5状態が数学的に最小であること。

したがって本実験の結論は、νという既知の隠れ状態1個の内部化に限定する。