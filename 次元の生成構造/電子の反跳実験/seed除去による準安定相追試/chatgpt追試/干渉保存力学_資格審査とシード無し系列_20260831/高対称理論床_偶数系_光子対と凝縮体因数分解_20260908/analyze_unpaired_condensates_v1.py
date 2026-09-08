#!/usr/bin/env python3
import csv, math, json
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np

HERE=Path(__file__).resolve().parent
from audit_highsym_even_photon_pairs_v2 import build_theoretical_floor, edges_for_n

EVEN_NS=list(range(4,41,2))
TOL=1e-10

def smallest_prime_factor(n):
    if n < 2: return None
    p=2
    while p*p<=n:
        if n%p==0: return p
        p+=1
    return n

def divisors_ge2(n):
    return [d for d in range(2,n+1) if n%d==0]

def key_phase(u, nd=10):
    return (round(float(u.real),nd), round(float(u.imag),nd))

def partition_class_regular_polygons(z, idxs, p):
    # Partition same-amplitude class into p-gons in the squared-wave plane.
    buckets=defaultdict(list)
    unit={}
    for i in idxs:
        w=z[i]*z[i]
        u=w/abs(w)
        k=key_phase(u)
        buckets[k].append(i); unit[k]=u
    groups=[]
    # deterministic repeated extraction
    while any(buckets.values()):
        k0=min((k for k,v in buckets.items() if v), key=lambda x:(math.atan2(x[1],x[0]),x))
        i0=buckets[k0].pop(0)
        u0=(z[i0]*z[i0])/abs(z[i0]*z[i0])
        group=[i0]
        ok=True
        for a in range(1,p):
            target=u0*np.exp(2j*np.pi*a/p)
            kt=key_phase(target)
            if not buckets.get(kt):
                # nearest fallback within numerical tolerance
                candidates=[kk for kk,v in buckets.items() if v]
                if not candidates:
                    ok=False; break
                kk=min(candidates,key=lambda q:abs(complex(q[0],q[1])-target))
                if abs(complex(kk[0],kk[1])-target)>1e-8:
                    ok=False; break
                kt=kk
            group.append(buckets[kt].pop(0))
        if not ok:
            raise RuntimeError(f'failed polygon partition p={p}, remaining class size={len(idxs)}')
        groups.append(group)
    return groups

def brute_n4_triads(z):
    from itertools import combinations
    unused=set(range(len(z))); groups=[]
    while unused:
        found=None
        for comb in combinations(sorted(unused),3):
            den=sum(abs(z[i])**2 for i in comb)
            res=abs(sum(z[i]*z[i] for i in comb))/den
            if res<TOL:
                found=comb; break
        if found is None: break
        groups.append(list(found)); unused-=set(found)
    return groups,sorted(unused)

def write_csv(path,rows):
    if not rows: return
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

def main():
    # Read selected maximum matching from prior audit.
    match_by_N=defaultdict(list)
    pth=HERE/'selected_maximum_matchings.csv'
    with pth.open(encoding='utf-8') as f:
        for r in csv.DictReader(f):
            match_by_N[int(r['N'])].append((int(r['edge_index_1']),int(r['edge_index_2'])))

    summary=[]; groups_out=[]; class_out=[]
    for N in EVEN_NS:
        edges,d,z,lam,H,nd=build_theoretical_floor(N); M=len(z)
        matched=set()
        for a,b in match_by_N[N]: matched.update([a,b])
        unmatched=sorted(set(range(M))-matched)
        total_cl=abs(np.sum(z*z))
        rem_cl=abs(np.sum(z[unmatched]*z[unmatched])) if unmatched else 0.0
        rem_norm=float(np.sum(np.abs(z[unmatched])**2)) if unmatched else 0.0
        rem_rel=rem_cl/max(rem_norm,1e-300) if unmatched else 0.0
        L=N//math.gcd(N,4)
        sp=smallest_prime_factor(L)
        possible=divisors_ge2(L)

        primitive_groups=[]
        leftover=[]
        if N==4:
            primitive_groups,leftover=brute_n4_triads(z)
            primitive_size=3 if not leftover and primitive_groups else None
            mode='cross-distance exceptional triads' if primitive_size else 'unresolved'
        elif unmatched:
            primitive_size=sp
            mode=f'within-distance regular {sp}-gons in squared-wave plane'
            # all unmatched == all waves for standard branch outside N=8k
            for dd in range(1,N//2+1):
                idxs=[i for i in unmatched if int(d[i])==dd]
                if idxs:
                    gs=partition_class_regular_polygons(z,idxs,sp)
                    primitive_groups.extend(gs)
            used=set(i for g in primitive_groups for i in g)
            leftover=sorted(set(unmatched)-used)
        else:
            primitive_size=2
            mode='complete two-wave ±i pairing'
            # represent prior matching as primitive groups
            primitive_groups=[list(x) for x in match_by_N[N]]
            leftover=[]

        # validate groups
        max_gres=0.0
        for gi,g in enumerate(primitive_groups):
            den=sum(abs(z[i])**2 for i in g)
            gres=abs(sum(z[i]*z[i] for i in g))/max(den,1e-300)
            max_gres=max(max_gres,gres)
            ds=sorted(set(int(d[i]) for i in g))
            groups_out.append({
                'N':N,'group_id':gi,'group_size':len(g),'distance_classes':';'.join(map(str,ds)),
                'edge_indices':';'.join(map(str,g)),
                'edges':';'.join(f'{edges[i][0]}-{edges[i][1]}' for i in g),
                'closure_relative':gres,
                'sum_z2_real':float(np.real(sum(z[i]*z[i] for i in g))),
                'sum_z2_imag':float(np.imag(sum(z[i]*z[i] for i in g)))
            })
        for dd in range(1,N//2+1):
            idx=[i for i in range(M) if int(d[i])==dd]
            csum=sum(z[i]*z[i] for i in idx)
            cnorm=sum(abs(z[i])**2 for i in idx)
            class_out.append({'N':N,'distance_class':dd,'wave_count':len(idx),
                              'closure_relative':abs(csum)/max(cnorm,1e-300),
                              'squared_phase_order_L':L,
                              'allowed_regular_polygon_sizes_dividing_L':';'.join(map(str,possible)) if possible else ''})

        summary.append({
            'N':N,'M':M,'L=N/gcd(N,4)':L,'smallest_prime_factor_L':sp if sp else '',
            'prior_pair_matching_pairs':len(match_by_N[N]),'prior_pair_coverage':2*len(match_by_N[N])/M,
            'unmatched_wave_count_after_max_2pair':len(unmatched),
            'unmatched_relative_closure':rem_rel,
            'primitive_condensate_size':primitive_size if primitive_size else '',
            'primitive_group_count':len(primitive_groups),
            'primitive_group_coverage_of_all_waves':sum(len(g) for g in primitive_groups)/M,
            'leftover_after_primitive_factorization':len(leftover),
            'max_primitive_group_closure_relative':max_gres,
            'alternative_regular_polygon_sizes':';'.join(map(str,possible)),
            'factorization_mode':mode
        })

    write_csv(HERE/'unpaired_condensate_summary.csv',summary)
    write_csv(HERE/'primitive_zero_closure_groups.csv',groups_out)
    write_csv(HERE/'distance_class_zero_closure.csv',class_out)
    meta={'even_N':EVEN_NS,'tolerance':TOL,'definition':'factor total high-sym theoretical-floor waves after maximal exact 2-wave pairing; each distance class tested and factored by regular root-of-unity polygons in z^2 plane'}
    (HERE/'CONDENSATE_METADATA.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    print('N  M  paircov unmatched relclosure primitive groups leftover alt')
    for r in summary:
        print(f"{r['N']:2d} {r['M']:3d} {r['prior_pair_coverage']:.1f} {r['unmatched_wave_count_after_max_2pair']:3d} {r['unmatched_relative_closure']:.2e} {str(r['primitive_condensate_size']):>3} {r['primitive_group_count']:3d} {r['leftover_after_primitive_factorization']:3d} {r['alternative_regular_polygon_sizes']}")

if __name__=='__main__': main()
