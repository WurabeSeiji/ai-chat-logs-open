# 補遺 A・B の数値の検算（元プログラムと同じ S10・f_fast・window_sum を使う）
import numpy as np
def S10(x): return np.sum([np.cos(n*x) for n in range(1,11)],axis=0)
phi_deg=np.linspace(-180,180,6000); phi=np.deg2rad(phi_deg)
f=S10(10*phi); a=np.sqrt(0.5)
def Y2(w,amp=1.0,p=phi):
    r=np.deg2rad(w); return (S10(10*p)*(S10(p-r)+amp*S10(p+r)))**2
def peak_at(w,amp,loc):  # 解析点での値
    p=np.deg2rad(np.array([loc])); return Y2(w,amp,p)[0]
print('1) 図のラベルの値（ピーク位置での厳密値）')
for w in [-180,-144,-36,0,36,144,180]:
    print(f'  w={w:5d}  対称: 左={peak_at(w,1,w):9.1f} 右={peak_at(w,1,-w):9.1f} | 非対称: φ=w {peak_at(w,a,w):9.1f} φ=-w {peak_at(w,a,-w):9.1f}')
print('2) 対称版で w と -w の Y² が一致するか（最大差）', max(np.abs(Y2(w)-Y2(-w)).max() for w in range(0,181,36)))
print('   非対称版で Y²(φ,w) と Y²(-φ,-w) が一致するか', max(np.abs(Y2(w,a)-Y2(-w,a)[::-1]).max() for w in range(0,181,36)))
print('3) 36° の倍数以外の位置での山の高さ（φ=w の値、対称版）')
for w in [-36,-30,-27,-24,-18,-12,-9,-6,-3,0]:
    print(f'  w={w:4d}  f_fast(φ=w)={S10(10*np.deg2rad(w)):7.3f}  Y²(φ=w)={peak_at(w,1,w):10.1f}  全体の最大={Y2(w).max():10.1f}')
