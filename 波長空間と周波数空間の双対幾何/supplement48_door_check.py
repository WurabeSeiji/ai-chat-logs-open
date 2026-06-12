# -*- coding: utf-8 -*-
# 補遺48 §1: 二つのハードルが「同じ扉」か。
# 時間逆写像の標的集合 vs 振幅Born重みの集合 を直接照合する。
# nu_t^2 = Sum(|k_i|+1/2)^2 (零点込み). 振幅: |z|^2 = a^2+b^2 = 二平方和.
from sympy import factorint

def is_sum_of_k_squares_classical(n, k):
    """2:Fermat / 3:Legendre / 4:Lagrange"""
    if n < 0: return False
    if k == 2:
        for p, e in factorint(n).items():
            if p % 4 == 3 and e % 2 == 1: return False
        return True
    if k == 3:
        m = n
        while m % 4 == 0 and m > 0: m //= 4
        return m % 8 != 7
    if k == 4: return True
    raise ValueError

K = 25
oblong = sorted({a*(a+1) for a in range(0, K+1)})

def sumset(vals, r, cap):
    s = {0}
    for _ in range(r):
        ns = set()
        for a in s:
            for v in vals:
                if a + v <= cap: ns.add(a + v)
        s = ns
    return s

CAP = 400
s4 = {x + 1 for x in sumset(oblong, 4, CAP)}          # nu_t^2 (4軸,零点込み)
odd_sq = sorted({(2*a+1)**2 for a in range(0, K+1)})
s3_times4 = sumset(odd_sq, 3, CAP*4)                   # 4*nu_t^2 (3軸)

print("nu_t^2 (4axis,zp) <=60:", sorted(x for x in s4 if x <= 60))
print("all odds covered:", set(range(1, CAP, 2)).issubset(s4))
s4l = sorted(s4)
gaps4 = [s4l[i+1]-s4l[i] for i in range(len(s4l)-1) if s4l[i+1] <= CAP]
print("gap max/min:", max(gaps4), min(gaps4))
print("3axis 4*nu_t^2 mod8:", sorted({x % 8 for x in s3_times4}))

N = 2000
c2 = sum(1 for n in range(1, N+1) if is_sum_of_k_squares_classical(n, 2))
c3 = sum(1 for n in range(1, N+1) if is_sum_of_k_squares_classical(n, 3))
print("density 2/3/4-sq in 1..%d:" % N, c2, c3, N)

odds = list(range(1, 200, 2))
inA = [n for n in odds if is_sum_of_k_squares_classical(n, 2)]
notA = [n for n in odds if not is_sum_of_k_squares_classical(n, 2)]
print("odd<200 that are 2-square:", len(inA), "/", len(odds))
print("not 2-square:", notA[:20])
print("137 2-square:", is_sum_of_k_squares_classical(137, 2), "odd:", 137 % 2 == 1)
