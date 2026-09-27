# 第八思考実験 対照実験：交差項ゼロの123状態 block-diagonal S(Ψ)

## 目的

第八思考実験で用いる長い状態ベクトル

$$
\Psi=(\Psi_{aa},\Psi_{ab},\Psi_{bb})^T\in\mathbb R^{123}
$$

をそのまま使い、相互作用配列の非対角ブロックを厳密に0とする。

$$
S_{\rm ctrl}(\Psi)=
\begin{pmatrix}
S_{aa}(\Psi_{aa})&0&0\\
0&S_{ab}(\Psi_{ab})&0\\
0&0&S_{bb}(\Psi_{bb})
\end{pmatrix}.
$$

更新は毎 microstep、愚直に

$$
\Psi' = S_{\rm ctrl}(\Psi)\Psi
$$

で行う。

この対照系が Paper 6 / Paper 7 の3分離系と一致すれば、123状態化と大行列化そのものは既存物理を変更しておらず、交差項を導入したときの変化を非対角ブロックへ帰属できる。

## 初期状態

n=1 条件を使用する。

- `aa`: Paper 7 self n1, `(C,D)=(0.91,0)`
- `ab`: Paper 6 n1 attractive, `(C,D)=(1.09,0.36)`
- `bb`: Paper 7 self n1, `(C,D)=(0.91,0)`

その他の41状態初期値は Paper 6 / 7 と同一。

## 実験

1. 41×41 の状態依存行列 `S41(z)` を明示構成し、`S41(z) @ z` が Paper 6/7 `transition(z)` を再現するか監査。
2. 3個を block diagonal に並べた123×123 `S123(Psi)` を構成。
3. 全非対角ブロックが厳密0であることを監査。
4. 1 macrostep = 11 microsteps を123×123 dense matrix productで実行。
5. 4000 macrosteps（1軌道）まで実行し、3つの独立 `transition` と逐次比較。
6. 保存済み Paper 6 / Paper 7 の first-macro raw CSV と比較。

## 再現

```bash
python run_paper8_control_block_diagonal.py
python verify_against_papers6_7.py
```

結果は `results_control/` に出力される。
