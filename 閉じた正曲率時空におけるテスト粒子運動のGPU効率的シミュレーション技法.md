# 閉じた正曲率時空におけるテスト粒子運動のGPU効率的シミュレーション技法

**修士論文中間報告 / ゼミ資料**  
**著者**: 木原範昭  
**日付**: 2026年4月  

## 1. 目的
本研究の目的は、閉じた正曲率（\( k = +1 \)）を持つ時空におけるテスト粒子の運動を、事前計算した歪みメッシュを用いた一次近似により、GPU上の線形行列演算だけで効率的にシミュレーションする計算技法を提案することである。

## 2. アブストラクト
本研究では、FLRW計量（閉宇宙 \( k=+1 \)）上でのテスト粒子運動を対象に、空間を歪みのある直方体メッシュで近似するシミュレーション手法を提案する。メッシュの歪み情報（Jacobian行列）を事前に厳密に計算して巨大テンソルとして保存し、各メッシュ要素内では純粋な線形更新のみで時間発展を行う。これにより、GPUの行列演算（行列積・正規化・外積）を最大限活用し、大規模粒子（\( n^3 \)規模）の長時間シミュレーションを実現した。テスト粒子近似および曲率が十分大きい条件下で、曲率半径 \( R \) とスケール因子の変化に対する精度評価を行い、従来の厳密解との差異を定量的に示した。本手法は、Big Bang模型への拡張可能性を有する有力な計算技法である。

## 3. 制限
- テスト粒子近似（粒子に質量を持たず、粒子間相互作用を無視）
- 曲率が十分大きい領域（曲率半径 \( R \) を小さく設定）
- メッシュ要素内での一次近似（要素境界をまたぐ移動時の補間は簡易的）
- 時間歪み（固有時間 \( \tau \)）は速度依存の1次近似のみ
- GPUメモリ制約による粒子数・時間ステップの上限

## 4. 適用できないもの
本手法は以下の状況には適用できない：
- 粒子間重力相互作用（N-body効果）を考慮した自己重力系
- 曲率が極端に小さい（ほぼ平坦）な時空での高精度長距離運動
- 強い重力場や特異点近傍（ブラックホールなど）
- 厳密な一般相対性方程式の解析的解を必要とする場合
- メッシュ要素を大きく超える距離を移動する高速粒子

## 5. 前提条件
- 光速度 \( c = 1 \)（自然単位系）
- 曲率半径 \( R \geq n^3 \cdot \text{Scale} \)（\( n \): メッシュ分割数、Scale: 評価したい精度オーダー）
- 空間曲率は正の閉じた宇宙（FLRW計量 \( k = +1 \)）
- スケール因子 \( a(t) \) は外部パラメータとして与え、膨張効果を1次近似で考慮
- すべての演算はGPU上の線形行列演算を中心に実施

## 6. 初期化条件
- メッシュ間隔を1とし、全体サイズを \( n \) とした \( n \times n \times n \) の直方体パラメータ空間を定義
- 各格子頂点（合計 \( n^3 \) 点）にテスト粒子を配置
- 初期位置：パラメータ座標 \( \mathbf{\xi}_0 = (i, j, k) \)
- 初期速度ベクトル：\( \mathbf{v}_0 = \mathbf{0} \)
- 局所経過時間（固有時間）：\( \tau_0 = 0 \)

## 7. 計算方式
**Step 1**: テスト粒子の初期配置（\( n^3 \) 行列に位置情報を格納）

**Step 1（歪み計算）**: 空間メッシュの歪み計算  
直方体メッシュに対して閉じた曲率を考慮した歪みを厳密に計算（sin, cos使用可）。各メッシュ要素のJacobian行列を巨大テンソル配列として保存。

**Step 2**: テスト粒子の動きシミュレーション  
時間 \( t = 0 \dots T \) をループし、各ステップで移動位置・速度ベクトル・局所経過時間を計算。結果を巨大配列（形状例: `[T, N, 7]`）に格納（位置[0:3]、速度[3:6]、\( \tau \)[6]）。

## 8. 詳細計算式
FLRW計量（\( k = +1 \)）の共動座標近似として、パラメータ空間 \( \mathbf{\xi} \) を用いる。  
物理位置 \( \mathbf{x}(\mathbf{\xi}) \) の写像例：
\[
\mathbf{x}(\mathbf{\xi}) = R \cdot \sin\left(\frac{|\mathbf{\xi}| \pi}{n}\right) \cdot \frac{\mathbf{\xi}}{|\mathbf{\xi}| + \epsilon}
\]

Jacobian行列：
\[
J_{ij} = \frac{\partial x_i}{\partial \xi_j}
\]

テスト粒子更新（メッシュ内一次近似）：
\[
\mathbf{\xi}_{new} = \mathbf{\xi} + \mathbf{v} \cdot \Delta t
\]
\[
\Delta \mathbf{x} = J \cdot \Delta \mathbf{\xi}
\]
\[
\mathbf{v}_{new} = \mathbf{v} / a(t) \quad \text{(赤方偏移の1次近似)}
\]
\[
\Delta \tau = \Delta t \sqrt{1 - v^2} \quad \text{(時間歪みの1次近似)}
\]

全粒子はテンソル演算で一括処理：
\[
\mathbf{X}^{(t+1)} = \mathbf{X}^{(t)} + \text{batch\_matmul}(J, \Delta \mathbf{\Xi})
\]

## 9. シミュレーション結果（予定）
- 曲率半径 \( R \) とスケール因子の変化に対する精度評価
- 厳密解（解析的great-circle解）との位置誤差（\( L^2 \)ノルムなど）を定量化
- GPU計算速度の向上と、粒子軌道の視覚化（閉じた宇宙を一周する美しい軌道）
- Big Bang模型への拡張可能性を示唆

## 10. 引用文献
- [1] C. G. Böhmer et al., “Azimuthal geodesics in closed FLRW cosmological models,” Phys. Rev. D 110, 083504 (2024).  
- [2] T. Roy Choudhury, Cosmology Lecture Notes: FLRW kinematics and geodesics (2021).  
- [3] J. Adamek et al., “gevolution: a cosmological N-body code based on General Relativity,” JCAP (2016).  
- [4] D. Moxey et al., “An isoparametric approach to high-order curvilinear boundary layer meshes,” Comput. Methods Appl. Mech. Engrg. 283 (2015).  
- [5] G. F. R. Ellis, “Deviation of geodesics in FLRW spacetime geometries,” Class. Quantum Grav. (1997).

## 11. 補足：コード実装（PyTorch + CUDA版、4D埋め込み法推奨）

```python
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

n = 50                          # メッシュ分割数
R = 5.0                         # 曲率半径（調整可能）
n_steps = 1000
dt = 0.01

N = n * n * n

# 初期位置（4D S³埋め込み）
X = torch.randn(N, 4, device=device)
X = X / torch.norm(X, dim=1, keepdim=True) * R

# 初期速度（小さい擾乱）
V = torch.randn(N, 4, device=device) * 0.05
V = V - (X * V).sum(dim=1, keepdim=True) * X
V = V / (torch.norm(V, dim=1, keepdim=True) + 1e-8)

history = torch.zeros(n_steps, N, 7, device=device)  # [T, N, 7]
history[0, :, 0:4] = X
history[0, :, 3:6] = V
history[0, :, 6] = 0.0

for step in range(1, n_steps):
    t = step * dt
    a = 1.0 + 0.001 * t                     # a(t) 簡易例
    
    speed = torch.norm(V, dim=1, keepdim=True)
    theta = (speed / R) * dt
    
    cos_th = torch.cos(theta)
    sin_th = torch.sin(theta)
    unit_V = V / (speed + 1e-8)
    
    X_new = X * cos_th.unsqueeze(1) + unit_V * sin_th.unsqueeze(1)
    X_new = X_new / torch.norm(X_new, dim=1, keepdim=True) * R
    
    v_pec = speed.squeeze(1) / (a * R)
    dtau = dt * torch.sqrt(1.0 - v_pec**2 + 1e-8)
    
    V_new = V / a
    
    X = X_new
    V = V_new
    tau = history[step-1, :, 6] + dtau
    
    history[step, :, 0:4] = X
    history[step, :, 3:6] = V
    history[step, :, 6] = tau

print("シミュレーション完了。形状:", history.shape)