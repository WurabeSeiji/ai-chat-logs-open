# -*- coding: utf-8 -*-
"""Are the unnormalized-run parent and the normalized-run parent symmetry-equivalent (vertex permutation x global phase)?
Output: results/parent_symmetry_check.json"""
import os, json, itertools, importlib.util, numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.join(HERE, "results")
def load_engine(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "program", name + ".py")); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); m.progress = lambda s: None; return m
eng_un, eng_n = load_engine("original_engine"), load_engine("normalized_engine"); N = 5
su = eng_un.LowRankSystem(N); vu = eng_un.make_parent(su, np.random.default_rng(40260721 + 1000 * N))[0]
sn = eng_n.LowRankSystem(N); vn = eng_n.make_parent(sn, np.random.default_rng(40260721 + 1000 * N))[0]
vu = vu / np.linalg.norm(vu)
edges = list(zip(su.ea.tolist(), su.eb.tolist())); idx = {e: k for k, e in enumerate(edges)}
best = None
for perm in itertools.permutations(range(N)):
    P = [idx[tuple(sorted((perm[a], perm[b])))] for a, b in edges]
    w = vn[P]
    for conj in (False, True):
        ww = np.conj(w) if conj else w
        # edge sign flips z_e -> -z_e are a symmetry (K' = D K D): choose D to maximize |<D ww, vu>| = sum of |components| when phases align
        prod = np.conj(ww) * vu; ph = np.angle(prod.sum())
        D = np.sign(np.real(prod * np.exp(-1j * ph))); D[D == 0] = 1
        ov = abs(np.vdot(D * ww, vu))
        if best is None or ov > best[0]: best = (float(ov), perm, conj, D.tolist())
out = {"best_overlap_after_vertex_permutation_sign_flips_optional_conjugation": best[0], "perm": list(best[1]), "conjugated": best[2], "sign_flips": best[3],
       "raw_overlap": float(abs(np.vdot(vn, vu))), "sorted_amplitudes_un": sorted(np.round(np.abs(vu), 8).tolist()), "sorted_amplitudes_norm": sorted(np.round(np.abs(vn), 8).tolist())}
json.dump(out, open(os.path.join(R, "parent_symmetry_check.json"), "w"), indent=1)
print(f"raw overlap {out['raw_overlap']:.6f} -> best after vertex permutation + edge sign flips{' + conj' if best[2] else ''}: {best[0]:.12f}, perm={best[1]}, flips={best[3]}")
print("sorted |v|: un", out["sorted_amplitudes_un"]); print("sorted |v|: nm", out["sorted_amplitudes_norm"])
