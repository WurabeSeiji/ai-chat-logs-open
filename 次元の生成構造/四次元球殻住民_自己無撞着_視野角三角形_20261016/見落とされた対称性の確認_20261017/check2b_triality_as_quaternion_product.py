"""確認 2b：確認 2 で見つけた位数 3 の三重対称性の行列 T が、四元数の掛け算 q → a q b（または a q̄ b）で
書けるか。a, b はフルヴィッツの単位 24 個から探す。"""
import json, itertools
import numpy as np


def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return np.array([a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
                     a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0])


conj = lambda q: q * np.array([1, -1, -1, -1])
T = np.array(json.load(open('results/check2_triality_d4.json'))['triality_matrix'])
units = [np.eye(4)[i] * s for i in range(4) for s in (1, -1)] + [np.array(s) / 2 for s in itertools.product((1, -1), repeat=4)]
E = np.eye(4)
found = []
for a in units:
    for b in units:
        for cj in (False, True):
            M = np.array([qmul(qmul(a, conj(e) if cj else e), b) for e in E]).T
            if np.allclose(M, T):
                found.append(dict(a=a.tolist(), b=b.tolist(), conjugate=cj,
                                  order_a=int(next(k for k in range(1, 13) if np.allclose(np.linalg.matrix_power(np.array([qmul(a, e) for e in E]).T, k), np.eye(4)))),
                                  order_b=int(next(k for k in range(1, 13) if np.allclose(np.linalg.matrix_power(np.array([qmul(e, b) for e in E]).T, k), np.eye(4))))))
r = dict(T=T.tolist(), det=float(np.linalg.det(T)), representations=found, n_found=len(found))
json.dump(r, open('results/check2b_triality_as_quaternion_product.json', 'w'), indent=1, ensure_ascii=False)
print(json.dumps(r, indent=1, ensure_ascii=False))
