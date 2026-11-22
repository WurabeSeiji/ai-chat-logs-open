#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
電子の候補(ベース波だけ)の検算。仕様は思考実験2 の §23.1〜§23.2。
octonion_pulse_map_O2_to_O_v1.py と同じフォルダで実行する: python3 octonion_electron_candidates_base_v1.py
格子点数は 256 に固定。
"""
import sys, numpy as np
sys.argv=[sys.argv[0],'256']
sys.path.insert(0,'.')
import octonion_pulse_map_O2_to_O_v1 as m
NG=m.NG; psi=m.psi; mul=m.mul; conj=m.conj
E=np.eye(8)
def zc(X):  # 4 つの複素面 (1,e7),(e1,-e6),(e2,e5),(e3,e4): 左から e^{e7 t} で全部 e^{it}
    return np.stack([X[:,0]+1j*X[:,7], X[:,1]-1j*X[:,6], X[:,2]+1j*X[:,5], X[:,3]+1j*X[:,4]],axis=1)
def wind(X,tol=1e-9):  # 各面に入っている巻数とその振幅
    Z=np.fft.fft(zc(X),axis=0)/NG; out=[]
    for k in range(4):
        idx=np.where(np.abs(Z[:,k])>tol)[0]
        out.append({int(i if i<NG//2 else i-NG):round(float(abs(Z[i,k])),4) for i in idx})
    return out
def rotC(w):  # (1,e7) 面を w 回巻く単位の波
    X=np.zeros((NG,8)); X[:,0]=np.cos(w*psi); X[:,7]=np.sin(w*psi); return X
def left(w,x):  return mul(rotC(w),np.broadcast_to(x,(NG,8)))
# 候補: C 部分 c0·e^{e7 w psi}, C^3 部分は 3 枚の面を -1 回巻く
def electron(w,c0,v):
    X=c0*rotC(w)
    V=np.zeros(8); V[1:4]=v
    return X+left(-1,V)
v=np.array([1,1,1])/np.sqrt(3)
for w in (-1,0,1):
    for c0 in (0.0,np.sqrt(0.5)):
        e=electron(w,c0,v*np.sqrt(1-c0**2))
        print("\n== C 面の巻数 w=%+d, c0^2=%.2f =="%(w,c0**2))
        print(" ノルム一定:",np.allclose((e**2).sum(1),1)," 面ごとの巻数{巻数:振幅} (1,e7),Q1面,Q2面,Q3面:",wind(e))
        re=e[:,0]; print(" 1+7: 実部の範囲 %.3f..%.3f"%(re.min(),re.max()))
        q=(e[:,:4]**2).sum(1); print(" (1+3)(1+3): 前半のパワーの範囲 %.3f..%.3f"%(q.min(),q.max()))
        print(" 自分との積 e·e :",wind(mul(e,e)))
        p=electron(-w,c0,v*np.sqrt(1-c0**2)); p=conj(e)
        print(" 共役との積 e·conj(e):",wind(mul(e,conj(e))))
        # 巻数が逆の相手(陽電子の候補): 面を +1 回
        V=np.zeros(8); V[1:4]=v*np.sqrt(1-c0**2); pos=c0*rotC(-w)+left(+1,V)
        print(" 巻数が逆の相手との積 e·p:",wind(mul(e,pos)))

print("\n\n===== 追加: 面 1 枚ずつの状態と、その積 =====")
def plane(k,wn):   # Q_k 面だけを wn 回巻く単位の波 (k=1,2,3)
    V=np.zeros(8); V[k]=1; return left(wn,V)
d1,d2,d3=plane(1,-1),plane(2,-1),plane(3,-1)
print("d1:",wind(d1)); 
print("d1·d2      :",wind(mul(d1,d2)))
print("d2·d1      :",wind(mul(d2,d1)))
print("(d1·d2)·d3 :",wind(mul(mul(d1,d2),d3)))
print("d1·(d2·d3) :",wind(mul(d1,mul(d2,d3))))
print("(d2·d1)·d3 :",wind(mul(mul(d2,d1),d3)))
print("d1·d1      :",wind(mul(d1,d1)))
print("d1·conj(d1):",wind(mul(d1,conj(d1))))
print("d1·plane(1,+1):",wind(mul(d1,plane(1,1))))
print("d1·plane(2,+1):",wind(mul(d1,plane(2,1))))
print("\n===== 候補 B: (1,e7) 面を -3 回 =====")
eB=rotC(-3); pB=rotC(3)
print("eB:",wind(eB)," eB·eB:",wind(mul(eB,eB))," eB·pB:",wind(mul(eB,pB)))
print("eB·d1:",wind(mul(eB,d1)),"  d1·eB:",wind(mul(d1,eB)))
print("\n===== 掛け方ごとの回転 (a=e^{e7 t}, t=0.3) : 各面の回転角/t =====")
t=0.3; a=np.zeros(8); a[0]=np.cos(t); a[7]=np.sin(t)
x=m.rng.normal(size=8)
def ang(y):
    z0=np.array([x[0]+1j*x[7],x[1]-1j*x[6],x[2]+1j*x[5],x[3]+1j*x[4]])
    z1=np.array([y[0]+1j*y[7],y[1]-1j*y[6],y[2]+1j*y[5],y[3]+1j*y[4]])
    return np.round(np.angle(z1/z0)/t,6), np.allclose(abs(z1),abs(z0))
print(" a(a x)a :",ang(mul(mul(a,mul(a,x)),a)))
print(" a x     :",ang(mul(a,x)),"  x a :",ang(mul(x,a)))
