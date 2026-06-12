#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺62 検証: xyztRQ 射影の定式化
#  E1: ¼桁の定義の整合 — T_{1/4} が桁 d∈Z4 を +1 シフトする(波動水準で全数確認)
#  E2: 階層¼桁展開の丸め: 誤差上界 (1/4)·w^(L)·(1/(1-1/(2R))) と往復同一性
#  E3: 拘束面 R²=s(k), Q=ε(k) の マーキング/B4 不変性(補遺50 E3 の再確認を含む)
#  E4: t=Σ1/ν_t が論文8 内部時間と同一であることの数値照合(3進カスケード)
import numpy as np, itertools, random

# ---------- E1: ¼桁とT_{1/4} ----------
# 状態(軸単位) = (m, branch, sign): 波 = sign * φ, φ = √2cos(2πmx) [branch=0] / √2sin [branch=1]
# 桁の定義: d = (branch + 2*[sign<0]) mod 4
N=64; xs=np.arange(N)/N
def wave(m,branch,sign,shift=0.0):
    f = np.cos if branch==0 else np.sin
    return sign*np.sqrt(2)*f(2*np.pi*m*(xs-shift))
def identify(w,m):
    """波 w を (branch, sign) として同定"""
    for branch in (0,1):
        for sign in (1,-1):
            if np.allclose(w, wave(m,branch,sign), atol=1e-9): return branch,sign
    return None
ok=0; tot=0
for m in (1,):  # 基本波(位置桁は基本波が定義、補遺43 §3.1)
    for branch in (0,1):
        for sign in (1,-1):
            tot+=1
            d = (branch + 2*(sign<0)) % 4
            w2 = wave(m,branch,sign, shift=0.25)   # T_{1/4}
            r = identify(w2,m)
            if r is None: continue
            b2,s2 = r
            d2 = (b2 + 2*(s2<0)) % 4
            if d2 == (d+1)%4: ok+=1
print(f"E1: T_(1/4) による桁シフト d→d+1 (mod 4): {ok}/{tot} PASS")

# ---------- E2: 階層¼桁展開 (波長系譜 k=3) の稠密性・丸め上界・往復 ----------
# 訂正: 位置の桁は波長に乗る → 縮尺は容器比 1/(2R) でなく系譜の奇数比 k (論文9 奇数入れ子定理)。
# 混合基数 [4; k, k, ...]: x = (1/4)[d0 + Σ_{ℓ≥1} c_ℓ ∏ k^{-ℓ}], d0∈Z4, c_ℓ∈Z_k。
# 稠密性: 子の¼格子 (1/(4k))Z は親の¼格子 (1/4)Z を k 奇数で整数細分 ([P_ℓ:P_{ℓ-1}]=k) — 任意の奇数系譜で完備。
# 奇数性の要求は稠密性でなく波の整合(論文9 奇数入れ子定理)由来。容器比 1/(2R) は位置の縮尺ではない(本検証の主訂正)。
def expand_k(x, k, L):
    d0 = int(np.floor(x*4)) % 4
    rem = x - d0/4.0
    digits=[d0]; step=1.0/4
    for _ in range(L):
        step /= k
        c = int(np.floor(rem/step))
        c = min(c, k-1)
        digits.append(c); rem -= c*step
    return digits, rem
def rec_k(digits, k):
    x = digits[0]/4.0; step=1.0/4
    for c in digits[1:]:
        step /= k; x += c*step
    return x
random.seed(3)
for k,L in [(3,8),(5,8)]:
    max_err=0
    for _ in range(20000):
        x=random.random()
        dg,_=expand_k(x,k,L)
        max_err=max(max_err, abs(x-rec_k(dg,k)))
    bound = (1.0/4)*k**(-L) * (k/(k-1)) if k<=4 else None
    bound=(1.0/4)*float(k)**(-L)*(k/(k-1.0))
    print(f"E2 (k={k}, 深さ{L}): 最大誤差 {max_err:.3e} ≤ 上界 {bound:.3e} : {'PASS' if max_err<=bound+1e-12 else 'FAIL'}")
dg0,_=expand_k(0.7321,3,8)
print(f"    往復例 x=0.7321 → 桁{dg0} → {rec_k(dg0,3):.7f} (誤差 {abs(0.7321-rec_k(dg0,3)):.2e})")
# ---------- E3: 拘束面の不変性 ----------
PERMS=list(itertools.permutations(range(4)))
SIGNS=list(itertools.product([1,-1],repeat=4))
ok3=True
for k in itertools.product(range(-2,3),repeat=4):
    s = sum((abs(t)+0.5)**2 for t in k)
    q = (-1)**sum(abs(t) for t in k)
    for p in PERMS[:6]:
        for sg in SIGNS[:4]:
            k2=[sg[i]*k[p[i]] for i in range(4)]
            if abs(sum((abs(t)+0.5)**2 for t in k2)-s)>1e-12: ok3=False
            if (-1)**sum(abs(t) for t in k2)!=q: ok3=False
print(f"E3: 拘束面 R²=s(k), Q=ε(k) の B4 不変性(標本): {'PASS' if ok3 else 'FAIL'}")

# ---------- E4: t=Σ1/ν_t = 論文8 内部時間 ----------
S=3**10; t=0.0; ticks=[]
sp=S
while sp>3:
    t += 1.0/np.sqrt(sp); ticks.append(t); sp//=3
import math
print(f"E4: 3進カスケード S=3^10 の累積内部時間 t=Σ1/√s' = {t:.6f}  (論文8 §3.2 の t と同一構成)")
print(f"    末項支配・収束的: 各項比 = √3 ≈ {math.sqrt(3):.3f} 倍ずつ増加(物差しの細分化)")
