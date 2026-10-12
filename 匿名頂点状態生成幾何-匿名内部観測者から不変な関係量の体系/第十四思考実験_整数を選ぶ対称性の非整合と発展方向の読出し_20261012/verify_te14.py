#!/usr/bin/env python3
"""Independent numerical checks for Thought Experiment 14 (2026-10-12).
No physical model is fitted here; this checks only published algebra and counting.
"""
import math
from itertools import product

sites = []
for c in product(range(-3,4), repeat=4):
    farthest_corner_r2 = sum((abs(v) + .5)**2 for v in c)
    if farthest_corner_r2 <= 9.0:
        sites.append(c)

inner = [c for c in sites if sum(v*v for v in c) <= 1]
assert len(sites) == 137
assert len(inner) == 9
assert len(sites)-len(inner) == 128

n, m = 124, 23
assert math.gcd(n,m) == 1
R = math.cos(math.pi*m/n)**2
alpha_inv_root = 4*math.pi/(1-R)**2
alpha_inv_geo = (137 + math.sqrt(137**2 + 2*math.pi**2))/2
alpha_inv_codata22 = 137.035999177

print('4D cells enclosed:', len(sites), '= 1 + 8 + 128')
print('R124_23:', format(R,'.16f'))
print('alpha_inverse_root:',format(alpha_inv_root,'.12f'))
print('alpha_inverse_geo:', format(alpha_inv_geo,'.12f'))
print('alpha_inverse_CODATA2022:',format(alpha_inv_codata22,'.12f'))
print('root - geometric:',format(alpha_inv_root-alpha_inv_geo,'.12f'))
print('geometric - CODATA:',format(alpha_inv_geo-alpha_inv_codata22,'.12f'))
print('root - CODATA:',format(alpha_inv_root-alpha_inv_codata22,'.12f'))
print('NOTE: The numerical proximity does not establish physical identification.')
