#!/usr/bin/env python3
"""Analyze Z4 edge phases and the induced 4-sheet covering of a static make_parent network.

Input: parent_static_Nxxxxx_makeparent_*.npz containing array `v` (M complex edges).
Assumption: edge ordering is np.triu_indices(N, k=1).
No time stepping is performed.
"""
import argparse
import itertools
import json
import math
import numpy as np


def z4_labels(v):
    phase = np.angle(v)
    # Estimate common phase modulo pi/2 using fourth powers.
    alpha = np.angle(np.mean(np.exp(4j * phase))) / 4.0
    q = np.mod(np.rint((phase - alpha) / (np.pi / 2.0)).astype(int), 4)
    err = np.angle(np.exp(1j * (phase - alpha - q * np.pi / 2.0)))
    return alpha, q, float(np.max(np.abs(err)))


def edge_map(N, q):
    ii, jj = np.triu_indices(N, k=1)
    return {(int(i), int(j)): int(x) for i, j, x in zip(ii, jj, q)}


def oriented_q(E, i, j):
    if i < j:
        return E[(i, j)]
    return (-E[(j, i)]) % 4


def triangle_fluxes(N, E):
    vals = []
    for i, j, k in itertools.combinations(range(N), 3):
        f = (oriented_q(E, i, j) + oriented_q(E, j, k) + oriented_q(E, k, i)) % 4
        vals.append((i, j, k, f))
    return vals


def fundamental_fluxes(N, E, root=0):
    # Spanning star rooted at root. Each non-tree edge (i,j) closes root-i-j-root.
    vals = []
    others = [x for x in range(N) if x != root]
    for i, j in itertools.combinations(others, 2):
        f = (oriented_q(E, root, i) + oriented_q(E, i, j) + oriented_q(E, j, root)) % 4
        vals.append((root, i, j, f))
    return vals


def generated_subgroup(fluxes):
    gens = {f for *_, f in fluxes}
    reachable = {0}
    changed = True
    while changed:
        changed = False
        for a in list(reachable):
            for g in gens:
                b = (a + g) % 4
                if b not in reachable:
                    reachable.add(b); changed = True
    return sorted(reachable)


def cover_components(N, E):
    # Undirected lifted graph on (vertex, sheet), sheet in Z4.
    adj = [[] for _ in range(4 * N)]
    def idx(i, t): return 4 * i + (t % 4)
    for (i, j), q in E.items():
        for t in range(4):
            a = idx(i, t); b = idx(j, t + q)
            adj[a].append(b); adj[b].append(a)
    seen = set(); sizes = []
    for s in range(4 * N):
        if s in seen: continue
        stack = [s]; seen.add(s); n = 0
        while stack:
            x = stack.pop(); n += 1
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); stack.append(y)
        sizes.append(n)
    return sorted(sizes, reverse=True)


def counts(vals):
    out = [0, 0, 0, 0]
    for *_, f in vals: out[f] += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('npz')
    ap.add_argument('--N', type=int, default=None)
    ap.add_argument('--json-out', default=None)
    args = ap.parse_args()
    d = np.load(args.npz)
    v = d['v']
    N = args.N or (int(d['n']) if 'n' in d else None)
    if N is None:
        M = len(v)
        N = int((1 + math.isqrt(1 + 8 * M)) // 2)
    if N * (N - 1) // 2 != len(v):
        raise ValueError('N and M are inconsistent')
    alpha, q, max_err = z4_labels(v)
    E = edge_map(N, q)
    tri = triangle_fluxes(N, E)
    fund = fundamental_fluxes(N, E)
    subgroup = generated_subgroup(fund)
    comps = cover_components(N, E)
    result = {
        'N': N, 'M': len(v), 'alpha': alpha, 'z4_max_phase_error_rad': max_err,
        'triangle_count': len(tri), 'triangle_flux_counts': counts(tri),
        'independent_cycle_count': len(fund), 'fundamental_flux_counts': counts(fund),
        'generated_subgroup_Z4': subgroup,
        'cover_vertex_count': 4 * N, 'cover_component_sizes': comps,
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
    if args.json_out:
        with open(args.json_out, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)

if __name__ == '__main__':
    main()
