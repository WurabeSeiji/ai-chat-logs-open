#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一万能運動関数 K v1 —— 調和閉鎖分散による並進

位置づけ
--------
統一万能関数は F（相互作用）・G（読出し）・D（次元読出し）の三つで、
**運動が無かった**。F の線形部にも媒介頂点にも分散が無く、どちらも
AB チャネル空間の実回転であるため状態に e^{iφ} が掛からない。結果として
パケットは静止したままで、位置の巻数が上がらない
（実測: CR0 v2 で T=2000・位置の総回転 0.005 周）。

運動は調和閉鎖から出る。帯 k に角速度 ω_k = k·ω₁ を与えると、
位置表示では形不変の剛体並進になる（ソリトン的）。これは加速度逆二乗則論文
（Concept DOI 10.5281/zenodo.21441081）§6 の調和閉鎖 |ω_n|Δθ_n = Ω と
同じ構造であり、運動項と逆二乗則が同一の起源を持つことを意味する。

正本
----
`万能相互作用多体接続_v1/run_kinetic_dispersion_demo_v3.py` の KE クラス:

    self.disp = np.exp(1j * np.arange(self.nreg) * omega1)[None, :]
    def step(self):
        super().step()
        self.C = self.C * self.disp

本モジュールはこれを関数として取り出したものである。**振る舞いを変えない。**
資格審査 `run_qualification_kinetic_v1.py`（Q-K1〜Q-K3）で正本との
ビット一致を確認してから使用すること（恒久ルール: 新メンバー追加時は
正本との対照テストを資格審査に追加してから使用）。

符号の規約（正本 v3 の記録より）
--------------------------------
位相因子 e^{+ik·ω₁} は **−方向** への並進である
（c_k → c_k·e^{−2πiks/N} が +s 並進）。したがって

    1步あたりの並進 = −ω₁·N/(2π) セル
    予言位置 = (x₀ − v·t) mod N ,  v = ω₁·N/(2π)

正本 v2 はこの符号を誤り、v3 で修正された。ここでも同じ規約を使う。

分岐と定数（規約 R0b）
----------------------
分岐: **なし**。定数: **なし**。ω₁ は呼出側が宣言する引数である。
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "k_dispersion",
    "k_translate_flat",
    "k_cells_per_step",
    "k_omega1_for_cells",
    "k_predicted_position",
]


def k_dispersion(C, omega1, axis=-1):
    """調和閉鎖分散を、**帯添字がフーリエ側にある**配列へ作用させる。

        C_k  →  C_k · exp(i·k·ω₁)

    多体エンジンのレジスタ行列 C（… × Nn）はこの形なので直接使える。

    引数
      C      : 複素配列。`axis` が帯（倍音）添字
      omega1 : 基本角速度 ω₁（帯 k の角速度は ω_k = k·ω₁）
      axis   : 帯添字の軸（既定 −1）

    戻り
      同じ形の複素配列（新しい配列・入力は破壊しない）
    """
    C = np.asarray(C)
    ks = np.arange(C.shape[axis])
    ph = np.exp(1j * ks * omega1)
    shape = [1] * C.ndim
    shape[axis] = -1
    return C * ph.reshape(shape)


def k_translate_flat(psi, omega1, n_chi, n_eta):
    """**位置表示**の状態（χ×η を平坦化した 1 次元配列）を並進させる。

    二体正本の状態はこの形（χ 512 × η 16 = 8192）である。

    **符号付き周波数を使うこと（v1 の修正点）。** 当初は `k_dispersion` を
    そのまま流用して生添字 0..N−1 で位相を掛けていたが、これは並進ではない。
    FFT では負周波数が後半（N/2..N−1）に折り返して入るため、周波数対
    (+j, −j) は添字 (j, N−j) に置かれる。生添字で位相を掛けると

        ψ̂_j ψ̂_{N−j} → ψ̂_j ψ̂_{N−j} · e^{ijω}·e^{i(N−j)ω} = … × e^{iNω}

    となり、**Σψ² に大域位相 e^{iNω} が乗る**。閉包 Σa²+Σb² は、衝突を1回
    通ると個別の Σa²・Σb² が非零になり和だけが保存される状態になるので、
    ここに大域位相が乗ると和が壊れる（実測: 1步で 1.8e-16 → 7.1e-09）。

    符号付き k なら指数が j+(−j)=0 で消え、Σψ² は厳密不変になる
    （実測: arg(Σψ²) の変化 0.000000・100步累積でも 1.9e-16）。

    引数
      psi    : 長さ n_chi·n_eta の複素配列（位置表示）
      omega1 : 基本角速度（= 1步あたりの回転角［rad］）
      n_chi  : χ 格子点数（円周）
      n_eta  : η 格子点数（毛円）

    戻り
      同じ長さの複素配列
    """
    F = np.fft.fft(np.asarray(psi).reshape(n_chi, n_eta), axis=0)
    k = np.fft.fftfreq(n_chi, d=1.0 / n_chi)      # ★ 符号付き周波数
    F = F * np.exp(1j * k * float(omega1))[:, None]
    return np.fft.ifft(F, axis=0).reshape(-1)


def k_cells_per_step(omega1, n):
    """1步あたりの並進速度［セル/步］。符号は −方向を正とする。"""
    return float(omega1) * float(n) / (2.0 * np.pi)


def k_omega1_for_cells(v_cells, n):
    """指定の並進速度［セル/步］を与える ω₁（k_cells_per_step の逆写像）。"""
    return 2.0 * np.pi * float(v_cells) / float(n)


def k_predicted_position(x0, omega1, t, n):
    """t 步後の予言位置［セル］: (x₀ − v·t) mod n 。"""
    return (float(x0) - k_cells_per_step(omega1, n) * float(t)) % float(n)
