# 初期化仕様 v2.0 見直し監査

**対象:** `universal_initialization_einstein_maxwell_xcts_spec_ja_v2.0_20260930.md`  
**上位:** `universal_interaction_operator_complete_spec_ja_v2.0_20260930.md`  
**日付:** 2026-09-30

## 1. 見直し理由

旧v1.0は「superposed Kerr–Newman + XCTS + electromagnetic constraints」という方針は定めていたが、実装に必要な未知量、free data、境界条件、solver、target matching、ADM→GH変換の一部が未定義だった。このため同じ仕様から複数の異なるinitializerが作れた。

## 2. interaction v2.0 から初期化へ逆流した必須条件

1. `D_i,V_*,P_*` は隠れ外部情報でなく `ψ_grid` state。
2. `Phi_iab=D_i g_ab` を同じ operator で作る。
3. `Pi_ij` は interaction の `K_ij(g,Pi,Phi)` と一致させる。
4. `Pi_0a` を任意の lapse/shift time derivative で決めず、`H+Gamma=0` から一意に決める。
5. `H=F`、`Theta=-beta^k D_k H` を固定して gauge-driver と整合。
6. Maxwell sign/orientation は interaction v2.0 と同一。
7. initial stress-energy と interaction reconstruction を相互監査。
8. identity state, unit state, time state を `ψ0` に全て格納。

## 3. 旧仕様の穴を埋めた項目

- boosted Kerr–Newman の明示式
- attenuation rule `w_A=zeta M_A/M_T d`
- `zeta=2` common baseline
- XCTS unknowns `psi, alpha psi, beta^i`
- EM unknown `phi`
- coupled XCTS + EM solve
- outer/inner BC
- zero-charge BC
- Newton-Raphson / GMRES / block-Jacobi
- forward-AD Jacobian action
- Armijo line search
- ADM momentum control
- spin control
- mass seed matching
- extremal no-substitution rule
- ADM→GH conversion
- interaction compatibility audit

## 4. 重要なルール修正

### extremal

`Q/M=1` でも MP へ切り替えない。同じ initializer が収束しなければ `INITIALIZATION_FAILED_EXTREMAL`。

### superextremal

regular KN gate を超えるものは reject。別モデルへ救済しない。

### zero charge

charge-fixing Neumann BC の `Qd/Qsp` が 0/0 になる neutral case は、total normal flux zero の連続極限を仕様化。

### grid

initializerはapproved `ψ_grid` を受け取り、同じ `D_i,V_*,P_*` を使用する。別 derivative operator を持たない。

## 5. source-derived / project-choice の分離

Mukherjee et al. 2022 から直接採ったものと、本研究で自由度を消すために固定したものを本文第43節で明示的に分離した。

この区別を曖昧にしない。

## 6. まだ別仕様で必要なもの

interaction v2.0 と同様、`ψ_grid` 自体の domain-map / interface / closed semi-discrete operator の**生成法**は別 `grid/operator construction spec` が必要。

本初期化仕様は `ψ_grid` を明示的入力とすることで隠れ自由度にはしていないが、正式実装開始には grid/operator spec の承認も必要。

## 7. 結論

旧v1.0の初期化部をそのまま使ってはならない。

正式な初期化正本候補は v2.0 とし、interaction v2.0 と一対で監査する。
