#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Q軸の証拠テスト: 励起次数 ε(k)=(-1)^Σ|k| の頂点簿記
# 端点保存が R 認可(ε中立)で救えない「第二の1単位不足」が出るかの全数検査
import itertools, numpy as np
from collections import defaultdict

def shell_cells(m):
    out=[]
    for k in itertools.product(range(-4,5),repeat=4):
        if abs(sum((abs(t)+0.5)**2 for t in k)-m)<1e-9: out.append(np.array(k))
    return out
SH={m: shell_cells(float(m)) for m in (1,3,5,7,9,11,13)}
def shape(v): return tuple(sorted(np.abs(v)))
def eps(v): return (-1)**int(np.abs(v).sum())

print("="*72)
print("E1: 殻×軌道の ε と裸ノルム²")
print("="*72)
shapes_by_shell={}
for m in (1,3,5,7,9,11,13):
    d=defaultdict(int)
    for v in SH[m]: d[shape(v)]+=1
    shapes_by_shell[m]=sorted(d)
    for sh,cnt in sorted(d.items()):
        v=np.array(sh)
        print(f"  s={m:2d}  形{sh}  x{cnt:3d}  ε={eps(v):+d}  |k|²={int((v**2).sum())}"
              f"  (q≡Σ|k| mod2 = {int(np.abs(v).sum())%2})")

CHANNELS={7:[(3,3,1)], 9:[(5,3,1),(3,3,3)], 11:[(7,3,1),(5,5,1),(5,3,3)],
          13:[(9,3,1),(7,5,1),(7,3,3),(5,5,3),(3,3,3,3,1)]}

print()
print("="*72)
print("E2: 端点 ε 保存の可否 (親軌道×チャネル×娘軌道割当の全数)")
print("="*72)
def channel_eps_options(ch):
    opts=set()
    for combo in itertools.product(*[ [eps(np.array(sh)) for sh in shapes_by_shell[m]] for m in ch ]):
        opts.add(int(np.prod(combo)))
    return opts
def line_closes(vp, ch):
    """補遺38基準: v_parent ∈ {±(va±vb)} (娘の対、origin=0含む)"""
    pools=[SH[m] for m in ch]
    # 全配置は重いので: 対ごとの和差の可能集合で判定 (originは±v_d)
    tgt={tuple(vp), tuple(-vp)}
    ms=list(ch)
    for i,j in itertools.combinations(range(len(ms)),2):
        for va in pools[i]:
            for vb in pools[j]:
                if ms[i]==1 and ms[j]==1: continue
                for x in (va+vb, va-vb):
                    if tuple(x) in tgt: return True
    # origin との対: 線 = ±v_d → ノルム一致が必要
    if 1 in ms:
        for i,m in enumerate(ms):
            if m==1: continue
            for vd in pools[i]:
                if tuple(vd) in tgt: return True
    return False
print("  親s 親形             ε親  チャネル        ε終集合     ε保存可? 線閉鎖(R可視)?")
stuck_eps=[]; stuck_line=[]
for s in (7,9,11,13):
    for psh in shapes_by_shell[s]:
        vp=np.array(psh); ep=eps(vp)
        for ch in CHANNELS[s]:
            opts=channel_eps_options(ch)
            ok_eps = ep in opts
            ok_line = line_closes(vp, ch) if len(ch)<=3 else line_closes(vp, ch)
            print(f"  {s:3d} {str(psh):16s} {ep:+d}  {str(ch):14s} {str(sorted(opts)):10s}"
                  f"  {'YES' if ok_eps else 'NO '}      {'YES' if ok_line else 'NO '}")
            if not ok_eps: stuck_eps.append((s,psh,ch))
            if not ok_line: stuck_line.append((s,psh,ch))
print()
print("  ε 保存不能 (全チャネルで): ")
allstuck=defaultdict(list)
for s in (7,9,11,13):
    for psh in shapes_by_shell[s]:
        if all((s,psh,ch) in [(a,b,c) for a,b,c in stuck_eps] for ch in CHANNELS[s]):
            allstuck[s].append(psh)
for s,v in allstuck.items(): print(f"    s={s}: {v}")
print()
print("="*72)
print("E3: R 取引は ε を救えない (借入単位は origin 型: ε=+1 中立)")
print("="*72)
print("  検算: ε保存則の同値形 — Σ|k_d|² ≡ |k_p|² (mod 2) (裸ノルム²パリティ保存)")
for s in (9,):
    for psh in shapes_by_shell[s]:
        vp=np.array(psh)
        for ch in CHANNELS[s]:
            # ε保存可能な娘軌道割当が存在するときの |k|² パリティ照合
            for combo in itertools.product(*[shapes_by_shell[m] for m in ch]):
                e=int(np.prod([eps(np.array(c)) for c in combo]))
                k2=sum(int((np.array(c)**2).sum()) for c in combo)
                if e==eps(vp):
                    print(f"   s=9 親{psh}(|k|²={int((vp**2).sum())}) ch={ch}: ε一致割当 {combo}"
                          f" Σ|k_d|²={k2} → パリティ{'一致✓' if k2%2==int((vp**2).sum())%2 else '不一致✗'}")
                    break
            else:
                print(f"   s=9 親{psh} ch={ch}: ε一致割当なし — R(ε中立)では端点不足を吸収不能")
