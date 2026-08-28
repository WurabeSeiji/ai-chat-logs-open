# -*- coding: utf-8 -*-
"""fixed_equimodular/ の全エンジン複製と followup 自前エンジンの make_parent を、等モジュラー自己無撞着版（3 段階、N5〜N20 で検証済）に置換する。
 段階1: 旧アルゴリズム（位相のみ K の固有モード反復、正規化なし）→ rank-2 親 v0（既存の make_parent を _make_parent_phase_only として保持）
 段階2: v0 から位相のみ K の線形回転 exp(ANGLE·K_phase) で |z|² 相対幅 < equi_tol まで時間発展（旧力学の等モジュラー・ヌル単体）
 段階3: 振幅込み K_amp の λ=−iσ_max モードへの混合反復で残差 < tol まで磨く（等モジュラーなら K_amp = r²K_phase）
 正規化なし。不収束は RuntimeError（実験を走らせない）。段階3 後に等モジュラーでなくなった解（縮退 rank-2 解への逸脱）は不採用。
 呼び出し側の tol/restarts（旧アルゴリズム用）は段階2/3 に使わず、N5〜N20 で検証した定数（L_pre=124, restarts=20, tol=1e-10）を用いる。差分は results/equimodular_parent_patches.diff"""
import os, glob, re
HERE=os.path.dirname(os.path.abspath(__file__)); FX=os.path.join(HERE,"fixed_equimodular"); log=[]
HELPERS='''
def _adjacency(sys_lr):
    """隣接辺マスク A（頂点を共有する辺対に 1）。"""
    m = sys_lr.m
    A = np.zeros((m, m))
    for i in range(m):
        share = (sys_lr.ea == sys_lr.ea[i]) | (sys_lr.ea == sys_lr.eb[i]) | (sys_lr.eb == sys_lr.ea[i]) | (sys_lr.eb == sys_lr.eb[i])
        share[i] = False
        A[i, share] = 1.0
    return A


def _K_amplitude_aware(A, v):
    K = A * np.imag(np.conj(v)[:, None] * v[None, :])
    np.fill_diagonal(K, 0.0)
    return K


def _K_phase_only(A, v):
    th = np.angle(v)
    K = A * np.sin(th[None, :] - th[:, None])
    np.fill_diagonal(K, 0.0)
    return K


def _selfconsistency_residual(K, v):
    kv = 1j * (K @ v)
    mu = np.vdot(v, kv) / np.vdot(v, v)
    return float(np.linalg.norm(kv - mu * v) / np.linalg.norm(v))


def _exp_step(K, z, angle):
    w, V = np.linalg.eigh(1j * K)
    return V @ (np.exp(-1j * angle * w) * (V.conj().T @ z))


EQUI_L_PRE = 124        # 段階2 の位相のみ線形回転の分割数（N5〜N20 検証時の値）
EQUI_RESTARTS = 20      # 段階1〜3 のリスタート上限（検証時の値）
EQUI_TOL = 1e-10        # 段階3 のスケール不変残差の収束判定（検証時の値。N=5 で 3.8e-11 到達）
EQUI_PRE_STEPS = 40000  # 段階2 の最大 step 数（N=3 は約 30000 step 必要。N≥4 は 1300〜5300 step で到達）
EQUI_SPREAD_TOL = 1e-9  # 段階2 の等モジュラー判定 |z|² 相対幅
EQUI_SPREAD_FINAL = 1e-8  # 段階3 後も等モジュラーであることの受理条件（縮退解への逸脱を不採用にする）


def make_parent(sys_lr, rng, iters=2000, beta=0.5, tol=None, restarts=None, pre_steps=EQUI_PRE_STEPS, equi_tol=EQUI_SPREAD_TOL):
    """等モジュラー自己無撞着親（EQUIMODULAR PARENT, 2026-08-29）。
    段階1: 旧アルゴリズムで rank-2 親 v0（_make_parent_phase_only、正規化なし）。
    段階2: v0 から位相のみ K の線形回転 exp(ANGLE·K_phase) で |z|² 相対幅 < equi_tol（等モジュラー）まで時間発展。
    段階3: 振幅込み K_amp(v) の λ=−iσ_max モードへの混合反復で残差 < tol。等モジュラーなら K_amp = r²K_phase（恒等式）。
    振幅の正規化なし（全体スケールは段階1の親のまま）。全リスタート失敗で RuntimeError。
    戻り値は旧 make_parent と同じ (v, residual, sigma_spectrum)。residual は振幅込み K に対するスケール不変残差。"""
    # 呼び出し側の tol/restarts は旧（位相のみ）アルゴリズム用の値なので段階2/3 には使わず、検証済みの定数を用いる
    tol = EQUI_TOL
    restarts = EQUI_RESTARTS
    A = _adjacency(sys_lr)
    m = sys_lr.m
    log = []
    for r in range(restarts):
        v0, res0, sig0 = _make_parent_phase_only(sys_lr, rng, iters=max(iters, 400), beta=beta, tol=1e-8, restarts=1)
        z = v0.copy()
        spread = np.inf
        for t in range(pre_steps):
            z = _exp_step(_K_phase_only(A, z), z, 2.0 * math.pi / EQUI_L_PRE)
            if t % 100 == 99:
                r2 = np.abs(z) ** 2
                spread = float((r2.max() - r2.min()) / r2.mean())
                if spread < equi_tol:
                    break
        progress(f"親構成(等モジュラー) restart={r+1} 段階2: {t+1} step で |z|² 相対幅={spread:.2e}")
        if spread >= equi_tol:
            log.append((r + 1, "stage2", spread))
            continue
        v = z
        res = _selfconsistency_residual(_K_amplitude_aware(A, v), v)
        for it in range(iters):
            if res < tol:
                break
            K = _K_amplitude_aware(A, v)
            w, U = np.linalg.eigh(1j * K)
            top = np.where(w <= w.min() + 1e-12 * max(1.0, abs(w.min())))[0]
            ov = [abs(np.vdot(U[:, jj], v)) for jj in top]
            u = U[:, top[int(np.argmax(ov))]]
            ph = np.vdot(u, v)
            if abs(ph) > 0:
                u = u * (ph / abs(ph))
            u = u * (np.linalg.norm(v) / np.linalg.norm(u))
            v = (1.0 - beta) * v + beta * u
            res = _selfconsistency_residual(_K_amplitude_aware(A, v), v)
        K = _K_amplitude_aware(A, v)
        res = _selfconsistency_residual(K, v)
        sig = np.sort(np.linalg.eigvalsh(1j * K))[::-1]
        sig = sig[sig > 1e-12]
        r2 = np.abs(v) ** 2
        spread = float((r2.max() - r2.min()) / r2.mean())
        log.append((r + 1, "stage3", res if spread < EQUI_SPREAD_FINAL else spread))
        if res < tol and spread < EQUI_SPREAD_FINAL:
            sys_lr.set_state(v)
            progress(f"親構成(等モジュラー) 収束 restart={r+1} 残差={res:.2e} |v|={np.linalg.norm(v):.6f} |z|²相対幅={(r2.max()-r2.min())/r2.mean():.2e} 非零辺={int(np.sum(np.abs(v)>1e-8))}/{m} σ_amp={sig[:3]}")
            return v, res, sig
    msg = "make_parent: 等モジュラー自己無撞着解が得られず（restarts=%d）: " % restarts + ", ".join(f"r{r}:{st}:{x:.2e}" for r, st, x in log)
    progress(msg)
    raise RuntimeError(msg)

'''
ENG=sorted(p for p in glob.glob(os.path.join(FX,"**","run_n_scaling_lowrank_v1_*.py"),recursive=True)+glob.glob(os.path.join(FX,"**","original_RAW_K_reference.py"),recursive=True) if "NORMALIZED_ORIGINAL" not in p)
for e in ENG:
    s=open(e,encoding="utf-8").read()
    assert s.count("def make_parent(sys_lr, rng, iters=400, beta=0.5, tol=1e-8, restarts=3):")==1, e
    s=s.replace("def make_parent(sys_lr, rng, iters=400, beta=0.5, tol=1e-8, restarts=3):","def _make_parent_phase_only(sys_lr, rng, iters=400, beta=0.5, tol=1e-8, restarts=3):")
    i=s.index("def zero_closure_kernel_seed(sys_lr, rng):"); s=s[:i]+HELPERS.lstrip("\n")+"\n"+s[i:]
    open(e,"w",encoding="utf-8").write(s); log.append(os.path.relpath(e,FX))
# followup 自前エンジン（1 行形式）
f=os.path.join(FX,"N5_dynamics_followup_theorems_and_stability_20260826","followup_dynamics_20260826","run_followup_experiments.py"); s=open(f,encoding="utf-8").read()
assert s.count("def make_parent(sys,rng,iters=400,beta=0.5,tol=1e-8,restarts=3):")==1
s=s.replace("def make_parent(sys,rng,iters=400,beta=0.5,tol=1e-8,restarts=3):","def _make_parent_phase_only(sys,rng,iters=400,beta=0.5,tol=1e-8,restarts=3):")
FU=HELPERS.replace("progress(","print(")
FU=FU.replace("        v0, res0, sig0 = _make_parent_phase_only(sys_lr, rng, iters=max(iters, 400), beta=beta, tol=1e-8, restarts=1)","        v0, res0, sig0, nit0 = _make_parent_phase_only(sys_lr, rng, iters=max(iters, 400), beta=beta, tol=1e-8, restarts=1)")
FU=FU.replace("            return v, res, sig\n","            return v, res, sig, it\n")
i=s.index("def f_relative_to_parent(Z,v):"); s=s[:i]+FU.lstrip("\n")+"\n"+s[i:]; open(f,"w",encoding="utf-8").write(s); log.append(os.path.relpath(f,FX))
open(os.path.join(HERE,"results","apply_equimodular_parent.log"),"w").write("\n".join(log)+"\n"); print("patched:",len(log),"files")
