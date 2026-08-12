#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N 掃引 1→20・三系（中性／電子型／シードなし）を同一測定系で走らせる v2

v1 系（3フォルダ・破棄済み）からの変更:
  1. **力学を F v1 に戻した**。恒久ルール「力学は unified_interaction_v1.py の
     UnifiedEngine／build_standard_universe のみ」に従う。v1 系が使っていた
     F v2 は仕様書 §1.8 で不採用（真空に 9.66e-15 の物質が湧く自己維持
     アーティファクト）。F v1 は資格審査 Q1–Q5 ALL PASS・G2 短絡あり。
  2. **r（混合率＝反射率）の記録を追加**。r = sin²θ = P奇/(P奇+P偶) は新柱6 の
     一行定式そのもので、α 根 R_α = cos²(23π/124) = 0.697177928 において
     透過率 1−R_α = sin²(23π/124) = √(4πα) ＝ 素電荷、α⁻¹ = 4π/(1−R_α)² = 137.036
     となる（本スクリプト起動時に検算して表示する）。
  3. **全16帯のパワー**を毎步記録。偶数帯（ボゾン帯）が何段目で立つか——3次選択則
     k* = k₁+k₂−k₃ から 1段目=k3（奇・δ²）／2段目=k4（偶・δ⁴）という予言の検証用。
  4. 三系を一つのスクリプトのモード切替にした（測定系の同一性を構造的に保証）。
  5. **モードごとに図を分けて保存**。無改変 import の `fig_one`/`fig_summary`/
     `fig_matrix` は固定名（`fig_nsweep_*_v1.png`）で出力しモード間で上書きされる
     ため、生成直後に `fig_{mode}_4panel_N{n}_v2.png`／`fig_{mode}_summary_v2.png`／
     `fig_{mode}_birth_matrix_v2.png` へ改名する（関数は書き換えない）。

**無改変 import の方針（シリーズ規約）**
中性掃引スクリプト run_tb_nsweep_1to20_v1.py を本フォルダへ無改変コピー
（md5 bfa5d854b637fe97c33a7148be9c7f86）し module として import して
`run_one` / `summarize` / `fig_one` / `fig_summary` / `fig_matrix` をそのまま使う。
走行ロジック・記録項目・判定・図は一行も書き換えない。コピーの HERE は本フォルダ
なので出力は本フォルダに出る（原本フォルダは無傷）。

**モード（唯一の物理的差）**
| mode            | δ      | primary (k_seed, m_seed) | primary (k*, m*) | 種 |
|-----------------|--------|--------------------------|------------------|----|
| neutral         | 1e-2   | (1, 0)                   | (3, 0)           | 中性フェルミオン（ν型） |
| electron        | 1e-2   | (1, 3)                   | (3, −3)          | 電子型 |
| vacuum          | 0      | (1, 3) 宣言のみ          | —                | ボゾン真空（シードなし） |
| fermion_family  | 1e-2×5 | k=1, η∈{0,1,3,5,6}    | k*=3, 5巻き    | フェルミオン積荷のみ |
| boson_family    | 1e-2×3 | k=6, η∈{0,3,5}        | k*=14, 3巻き   | ボゾン積荷のみ |

ポンプは全モードで k=2。primary seed は既存4モードと
fermion_family で k=1（相棒 k*=3）、boson_family で k=6（相棒 k*=14）。

条件（三系で完全同一・すべて宣言値）:
  F=unified_interaction_v1.py（★v1）/ D=unified_dimension_v1.py /
  G=unified_readout_v3.py / S=selection_v1.py
  Nn=16・Nη=8・T=4000・親 seed=2 固定・cell=(2,0)・order=6・窓[2000,4000]・N=1..20
  各 N で真空対照（δ=0）も走行。mode=vacuum では2条件が同一になるため決定論検査に転用。

事前登録（実行前固定）:
  (V1) 構成できない N は三系で同一。シードの巻き・振幅は make_parent／build_init に
       影響しないため。不一致なら異常。
  (V2) 空間: crossing（f₂>0.05）が全 N で起きる。τ_space を三系で比較する。
  (V3) 物質: neutral/electron は f_seed が立つ。vacuum は全奇数帯パワーが全步で厳密 0。
  (V4) 電子型: 相棒帯 k=3 の巻き集中度 ≥0.9・優勢巻き m̂=−3・Q̂=−1・可読率 1。
  (V5) **r の挙動**: r（node 平均）の全步時系列と α 根との距離を記録する。
       vacuum は P奇=0 ゆえ r=0（真空固定点）が厳密に成立するはず。
       neutral/electron の r が α 根 0.697177928 に向かうか・通り過ぎるかは
       **予言を置かず記録する**（本実験の主目的）。
  (V6) **偶数帯の立ち上がり**: 帯別パワーから k=4 等の偶数帯が立つ時刻を読む。
       立つなら δ⁴ 抑制のオーダーか（δ 掃引は次段）。

**mode="mixed"（追加）**——実際の時空に近い多種共存。周期表の 62 状態は (k,η) 住所と
しては 7〜8 個に縮約される（u,c,t は全て m=+2 等・世代は別軸）ので、8 セルでほぼ全ての
生成可能な住所を覆う。和則の逆解き k_seed = 2k_pump − k*、m_seed = 2m_pump − m*:
  シードA（点火＋フェルミオン積荷）k=1・η∈{0,1,3,5,6}・各δ=1e−2
    → 相棒帯 k*=3（奇）・m*∈{0,−1,−3,+3,+2} ＝ ν型・d型・e型・陽電子型・u型
  シードB（ボゾン積荷）           k=6・η∈{0,3,5}・各δ=1e−2
    → 相棒帯 k*=14（偶）・m*∈{0,−3,+3} ＝ γ/Z型・W−・W+
  総物質分率 f ≈ Σδ² = 8e−4（量制御窓 [1e−3,1e−1] の下端付近・氾濫しない）
  k_seed=4 は相棒が零モード k*=0 になるため禁止・k_seed=10 は相棒が自分自身。
既知の難点（走行前に明示）: ①交差項 s_i+s_j−s_k で帳簿が密になり、単一シードで測った
排他比 6.1e289 は成り立たない（V7 で定量する）②u型・d型（3∤m）は載っても単独で
電荷が読めない（閉じ込め）③γ と Z は (k,η) では縮退（質量で分かれる）④G・H・g8重項・
世代軸はどんなシードでも作れない。

追加記録: 128セル帳簿 P[k,η] のスナップショット（50步毎）・狙った相棒セルごとのパワー・
狙っていないセルの総和・排他比（V7）。

**長時間走行（T 掃引）**: 第4引数で T を上書きできる。窓は後半 [T//2, T]（宣言値）。
決定論なので T=42000 の1本の軌道が短い T をすべて前半として含む——T 掃引は独立走行を
並べるのではなく、1本の長走行から窓を切り出して読む（同一軌道なので厳密に比較できる）。
T≠4000 のときは出力名に `_T{T}` が付き、既存の T=4000 成果物を上書きしない。

**シード強度 δ の掃引（第5引数）**: セル一覧の δ を一律で上書きする。
「どこまで弱くしても同じ発展が起こるか」を測る。出力名に `_d{δ}` が付く。
走行前に登録した予言（δ=1e−3・T=42000）:
  ① τ=1 の r_nopump = 5δ²/(5δ²+3δ²) = **5/8 = 0.625000 で δ に依存しない**
     （セル数の比だけで決まる）
  ② 頂点の駆動 R = sin²θ ≈ P奇/P偶 ≈ 5δ² なので δ→1/10 で **R→1/100**。
     T=42000 で r_nopump が動く量は 1/100（距離の縮みが 8% → 0.08%）。
     「同じ発展」を見るには **T を 100 倍**必要という帰結になる
  ③ τ_space は真空側へ戻る（既知: δ=0 で 1624・δ=1e−2 1セルで 224・
     δ=1e−2 8セルで 74 → δ=1e−3 8セルは 74 と 1624 の間）
  ④ f ≈ Σδ² = 8×1e−6 = 8e−6（δ² 則の検証点）

使い方: python3 run_nsweep_three_series_v2.py <mode> [Nmin Nmax] [T] [δ] [output_suffix]
        mode ∈ {neutral, electron, vacuum, mixed, fermion_family, boson_family}
        （省略時 1 20・T=4000・δ=1e−2）
        output_suffix は安全な ASCII 識別子のみ。`_rep-<suffix>` を全成果物に
        付け、既存の正本を上書きせず完全同条件の複製走行を可能にする。
"""
from __future__ import annotations

import importlib.util
import json
import math
import os                      # ダンプ版 v2 の追加（環境変数で二段サンプリング指定）
import re
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
UF = HERE.parent / "統一万能関数_v1"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# 走行ロジックの提供元（無改変コピー）
ns = _load("ns_sweep", HERE / "run_tb_nsweep_1to20_v1.py")
assert ns.HERE == HERE, "コピーの HERE が本フォルダでない（出力先が原本になる）"

# ★力学は F v1（恒久ルール）
F1 = _load("f1_three", UF / "unified_interaction_v1.py")

# --- 宣言値 ----------------------------------------------------------------
D0 = 1e-2
MODES = {
    "neutral":  {"delta": D0, "m_pump": 0, "m_seed": 0, "label": "中性フェルミオン",
                 "cells": [(1, 0, D0)]},
    "electron": {"delta": D0, "m_pump": 0, "m_seed": 3, "label": "電子型フェルミオン",
                 "cells": [(1, 3, D0)]},
    "vacuum":   {"delta": 0.0, "m_pump": 0, "m_seed": 3, "label": "ボゾン真空（シードなし）",
                 "cells": []},
    # 混合シード: 実際の時空に近い多種共存。周期表の生成可能な住所をほぼ全て狙う。
    #   奇数帯 k=1（点火＋フェルミオン積荷）→ 相棒帯 k*=4−1=3
    #   偶数帯 k=6（ボゾン積荷）           → 相棒帯 k*=4−6=−2≡14
    #   m* = −m_seed (mod 8)
    "mixed":    {"delta": D0, "m_pump": 0, "m_seed": 0, "label": "混合シード（多種共存）",
                 "cells": [(1, 0, D0), (1, 1, D0), (1, 3, D0), (1, 5, D0), (1, 6, D0),
                           (6, 0, D0), (6, 3, D0), (6, 5, D0)]},
    # mixed から偶数帯積荷を除いた family ablation。primary は先頭セル。
    "fermion_family": {
        "delta": D0, "m_pump": 0, "m_seed": 0, "k_seed": 1,
        "label": "フェルミオン積荷のみ（5セル）",
        "cells": [(1, 0, D0), (1, 1, D0), (1, 3, D0), (1, 5, D0), (1, 6, D0)],
    },
    # mixed から奇数帯積荷を除いた family ablation。P_odd=0 が保たれるかを測る。
    "boson_family": {
        "delta": D0, "m_pump": 0, "m_seed": 0, "k_seed": 6,
        "label": "ボゾン積荷のみ（3セル）",
        "cells": [(6, 0, D0), (6, 3, D0), (6, 5, D0)],
    },
}
K_PUMP = 2
NETA, NN = ns.NETA, ns.NN
ODD_K = [k for k in range(NN) if k % 2 == 1]
EVEN_K = [k for k in range(NN) if k % 2 == 0 and k != 0]
THETA_ALPHA = 23 * math.pi / 124
R_ALPHA = math.cos(THETA_ALPHA) ** 2          # 0.697177928
T_ALPHA = math.sin(THETA_ALPHA) ** 2          # 0.302822072
ALPHA_INV = 4 * math.pi / T_ALPHA ** 2        # 137.036043

SPECIES = {0: "ν型/γ・Z型", -1: "d型", -3: "e型/W−", 2: "u型", 3: "陽電子型/W+"}

def _signed(idx: int) -> int:
    h = NETA // 2
    return int(((idx + h) % NETA) - h)


MODE = sys.argv[1] if len(sys.argv) > 1 else "neutral"
assert MODE in MODES, f"mode は {list(MODES)} のいずれか"
CFG = MODES[MODE]
DELTA, M_PUMP, M_SEED = CFG["delta"], CFG["m_pump"], CFG["m_seed"]
K_SEED = int(CFG.get("k_seed", 1))
PARTNER_K = (2 * K_PUMP - K_SEED) % NN
M_STAR = 2 * M_PUMP - M_SEED
M_STAR_IDX = M_STAR % NETA
ns.DELTA = DELTA                   # 掃引スクリプトの主条件（0 か否かの点火フラグとして働く）
# 長時間走行の宣言: 第4引数で T を上書きする。窓は後半 [T//2, T]（宣言値）。
if len(sys.argv) > 4:
    ns.T = int(sys.argv[4])
    ns.WIN = (ns.T // 2, ns.T)
T_TAG = "" if ns.T == 4000 else f"_T{ns.T}"
# シード強度 δ の上書き（第5引数）。セル一覧の δ を一律で置き換える。
# 「どこまで弱くしても同じ発展が起こるか」の測定用。出力名に _d{δ} が付く。
D_TAG = ""
if len(sys.argv) > 5:
    _dnew = float(sys.argv[5])
    CFG["cells"] = [(k, e, _dnew) for (k, e, _d) in CFG["cells"]]
    if CFG["delta"] > 0:
        CFG["delta"] = _dnew
        DELTA = _dnew
        ns.DELTA = _dnew
    D_TAG = f"_d{_dnew:g}"
# 第6引数は物理条件ではなく、完全同条件の複製走行用識別子。
# パス区切り・空白・シェル記号を入れず、既存正本の名前と分離する。
OUTPUT_SUFFIX = sys.argv[6] if len(sys.argv) > 6 else ""
if OUTPUT_SUFFIX:
    assert re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", OUTPUT_SUFFIX), (
        "output_suffix は ASCII 英数字で始まる 1〜64 文字"
        "（英数字・_ ・- のみ）"
    )
OUTPUT_TAG = f"_rep-{OUTPUT_SUFFIX}" if OUTPUT_SUFFIX else ""
TAG = T_TAG + D_TAG + OUTPUT_TAG
CELLS = CFG["cells"]               # シードセル [(k, η, δ), ...]
if CELLS:
    assert CELLS[0][0] == K_SEED, "primary k_seed と先頭シードセルが不一致"
    assert CELLS[0][1] % NETA == M_SEED % NETA, (
        "primary m_seed と先頭シードセルが不一致"
    )
# 狙った相棒セル: k* = 2·K_PUMP − k_seed, m* = 2·M_PUMP − m_seed（重複は除く）
TARGETS = []
for (_k, _e, _d) in CELLS:
    t = ((2 * K_PUMP - _k) % NN, (2 * M_PUMP - _signed(_e)) % NETA)
    if t not in TARGETS:
        TARGETS.append(t)
CELL_RECORDS = [
    {"k": int(k), "eta_index": int(e % NETA), "m_signed": _signed(e),
     "delta": float(d)}
    for (k, e, d) in CELLS
]
TARGET_RECORDS = [
    {"k": int(k), "eta_index": int(e % NETA), "m_signed": _signed(e),
     "species": SPECIES.get(_signed(e), "?")}
    for (k, e) in TARGETS
]
N_DELTA = float(sum(d for (_k, _e, d) in CELLS))
N_DELTA2 = float(sum(d * d for (_k, _e, d) in CELLS))
NF = sum(1 for (k, _e, _d) in CELLS if k % 2 == 1)
NB = sum(1 for (k, _e, _d) in CELLS if k % 2 == 0)
PF_DECLARED = float(sum(d * d for (k, _e, d) in CELLS if k % 2 == 1))
PB_DECLARED = float(sum(d * d for (k, _e, d) in CELLS if k % 2 == 0))
PSEED_DECLARED = PF_DECLARED + PB_DECLARED
R_NP_CELL_RATIO = (float(NF / (NF + NB)) if (NF + NB) else None)
R_NP_POWER_RATIO = (float(PF_DECLARED / PSEED_DECLARED)
                    if PSEED_DECLARED > 0 else None)
LEDGER_EVERY = 50                  # 128セル帳簿のスナップショット間隔（步）

_ENGINES: list = []


REC_KEYS = ("r_mean", "r_med", "r_min", "r_max", "r_raw",
            "dist_alpha", "absdist_alpha",
            "r_nopump", "dist_alpha_nopump", "even_power_nopump",
            "conc_all", "conc_k3", "dom_m", "q_hat", "readable", "k3_power",
            "odd_power", "odd_amp_max", "even_power", "total_power",
            "target_power", "nontarget_power", "excl_ratio",
            # 既存23配列の後ろにだけ追加し、common arrays の順序と値を保つ。
            "seed_power", "pump_power", "primary_seed_power")


# ===== ダンプ版 v2 の追加（0/3）===============================================
# 二段サンプリングの設定。CLI を壊さないよう環境変数で与える。
#   DUMP_TAUC   : ここまでは毎步保存（既定 4000）
#   DUMP_STRIDE : それ以降の間引き幅（既定 31）
# 実測根拠: 転移域は τ=2300〜2600、振幅の減衰時定数は約 1190 步。
# τ_c=4000 は転移後 1400 步（時定数の1.2倍）まで密に覆う。
# 間引き 31 は減衰時定数あたり 38 点で、包絡線を追うのに十分。
DUMP_TAUC = int(os.environ.get("DUMP_TAUC", "4000"))
DUMP_STRIDE = int(os.environ.get("DUMP_STRIDE", "31"))


def _dump_index(t: int):
    """步 t の保存先フレーム番号。保存しない步は None。"""
    if t < DUMP_TAUC:
        return t
    d = t - DUMP_TAUC
    return DUMP_TAUC + d // DUMP_STRIDE if d % DUMP_STRIDE == 0 else None


def _dump_frames(T: int) -> int:
    """T 步走ったときに保存されるフレーム総数。"""
    if T <= DUMP_TAUC:
        return T
    return DUMP_TAUC + (T - 1 - DUMP_TAUC) // DUMP_STRIDE + 1


def _dump_taus(T: int) -> np.ndarray:
    """各フレームが対応する τ の一覧。"""
    return np.array([t for t in range(T) if _dump_index(t) is not None],
                    dtype=np.int64)
# ==============================================================================


class RecordingEngine(F1.UnifiedEngine):
    """力学は F v1 の UnifiedEngine そのまま。step 後に記録を 1 回だけ追記する。

    _readout() は純関数（self.C を読むだけ・副作用なし）なので、毎步もう一度
    呼んで node ごとの R = scale·sin²θ を取得する。scale=1 なので R = sin²θ = r。
    """

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._rec: list = []
        self._bands: list = []
        self._ledger: list = []        # (k,η) 128セル帳簿のスナップショット
        self._ledger_t: list = []
        self._tgt: list = []           # 狙ったセルごとのパワー（毎步）
        self._t = 0
        self._dump = None              # ダンプ版 v2 の追加：C2 memmap

    def step(self, *a, **kw):
        out = super().step(*a, **kw)
        C2 = self.C2()
        # ===== ダンプ版 v2 の追加（1/3）=======================================
        # 二段サンプリングで C2（M×Nn×Nη 複素）を memmap へ書き出す。
        #   τ <  DUMP_TAUC : 毎步（転移を含む区間を密に取る）
        #   τ >= DUMP_TAUC : DUMP_STRIDE 步ごと（包絡線の減衰を追う）
        # 既存の記録・判定・図には一切触れない（読むだけ）。
        if self._dump is not None:
            k = _dump_index(self._t)
            if k is not None and k < self._dump.shape[0]:
                self._dump[k] = C2
        # =====================================================================
        A = np.abs(C2)
        P = A ** 2
        tot = float(P.sum())
        Pk = P.sum(axis=(0, 2))                       # 帯ごとのパワー（長さ Nn）
        odd_p = float(Pk[ODD_K].sum())
        even_p = float(Pk[EVEN_K].sum())
        r_raw = odd_p / (odd_p + even_p) if (odd_p + even_p) > 0 else float("nan")
        # 物質側の r: ポンプ帯（凝縮体）を分母から除く＝H1c が種に対して測った量と同型
        even_np = float(Pk[[k for k in EVEN_K if k != K_PUMP]].sum())
        r_nopump = (odd_p / (odd_p + even_np)
                    if (odd_p + even_np) > 0 else float("nan"))
        Rn = np.asarray(self._readout(), dtype=float)  # = scale·sin²θ（node ごと）
        r_mean = float(Rn.mean()) if Rn.size else float("nan")
        Pe = P.sum(axis=(0, 1))
        conc_all = float(Pe[M_STAR_IDX] / tot) if tot > 0 else float("nan")
        # `conc_k3` / `k3_power` は既存 NPZ との互換キー名。値は常に
        # 設定された primary 相棒帯を読む（既存4モードで k=3、boson_family で k=14）。
        Ppartner = P[:, PARTNER_K, :].sum(axis=0)
        s3 = float(Ppartner.sum())
        conc_k3 = (float(Ppartner[M_STAR_IDX] / s3)
                   if s3 > 0 else float("nan"))
        dom = _signed(int(np.argmax(Ppartner))) if s3 > 0 else 0
        Pkh = P.sum(axis=0)                       # (Nn, Nη) 帳簿
        tgt_vals = [float(Pkh[k, e]) for (k, e) in TARGETS]
        tgt_p = float(sum(tgt_vals))
        seed_p = float(sum(Pkh[k, e % NETA] for (k, e, _d) in CELLS))
        pump_p = float(Pkh[K_PUMP, M_PUMP % NETA])
        primary_seed_p = float(Pkh[K_SEED, M_SEED % NETA])
        nontgt_p = max(tot - pump_p - seed_p - tgt_p, 0.0)
        excl = tgt_p / nontgt_p if nontgt_p > 0 else float("inf")
        self._rec.append((
            r_mean, float(np.median(Rn)), float(Rn.min()), float(Rn.max()), r_raw,
            r_mean - R_ALPHA, abs(r_mean - R_ALPHA),
            r_nopump, r_nopump - R_ALPHA, even_np,
            conc_all, conc_k3, float(dom), dom / 3.0,
            1.0 if (s3 > 0 and dom % 3 == 0) else 0.0, s3,
            odd_p, float(A[:, ODD_K, :].max()), even_p, tot,
            tgt_p, nontgt_p, excl,
            seed_p, pump_p, primary_seed_p))
        self._bands.append(Pk.copy())
        self._tgt.append(tgt_vals)
        self._t += 1
        if (self._t % LEDGER_EVERY) == 0 or self._t == 1:
            self._ledger.append(Pkh.copy())
            self._ledger_t.append(self._t)
        return out


def build_universe(n, delta, Nn=5, Neta=8, seed=2):
    """F v1 の build_standard_universe と同一手順。帯・巻きの置き場所のみ宣言化。"""
    m = n * (n - 1) // 2
    _, v, _, _, _, _, _, Z0c, wp0 = F1.abl.build_init(n, False)
    r2 = F1.gen3.make_parent(n, seed=seed)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C2_0 = np.zeros((m, Nn, Neta), complex)
    ip = M_PUMP % Neta
    C2_0[:, K_PUMP, ip] = Z0c
    # delta==0 は真空対照の呼び出し（掃引スクリプトが各 N で必ず 1 回行う）。
    # その場合シードを一切置かない。delta>0 なら宣言したセルを全て置く。
    if delta > 0:
        for (kk, ee, dd) in CELLS:
            C2_0[:, kk, ee % Neta] += dd * seed_state
    p2 = C2_0[:, K_PUMP, ip].real / np.linalg.norm(C2_0[:, K_PUMP, ip].real)
    q2 = C2_0[:, K_PUMP, ip].imag - (C2_0[:, K_PUMP, ip].imag @ p2) * p2
    with np.errstate(divide="ignore", invalid="ignore"):
        q2 = q2 / np.linalg.norm(q2)
    eng = RecordingEngine(n, C2_0, wp0)
    # ===== ダンプ版 v2 の追加（2/3）===========================================
    # 走行ごとに C2 の memmap を用意し、初期条件と τ 対応表を併記する。
    # delta>0 が物質側（m）、delta==0 が真空対照（v）。
    _side = "m" if delta > 0 else "v"
    _stem = f"{MODE}{TAG}_N{n}_{_side}"
    _nf = _dump_frames(int(ns.T))
    _taus = _dump_taus(int(ns.T))
    assert len(_taus) == _nf, f"フレーム数不一致 {len(_taus)} vs {_nf}"
    eng._dump = np.lib.format.open_memmap(
        HERE / f"dump_C2_{_stem}_v1.npy", mode="w+",
        dtype=np.complex128, shape=(_nf, m, Nn, Neta))
    np.savez(HERE / f"dump_meta_{_stem}_v1.npz",
             p2=p2, q2=q2, C2_0=C2_0, Z0c=Z0c, wp0=wp0,
             n=n, m=m, Nn=Nn, Neta=Neta, T=int(ns.T), delta=float(delta),
             seed=seed, k_pump=K_PUMP, m_pump=M_PUMP,
             partner_k=PARTNER_K, m_star_idx=M_STAR_IDX,
             odd_k=np.array(ODD_K), even_k=np.array(EVEN_K),
             cells=np.array(CELLS, dtype=float).reshape(-1, 3),
             dump_tauc=DUMP_TAUC, dump_stride=DUMP_STRIDE,
             dump_frames=_nf, dump_taus=_taus)
    print(f"    [dump] {_stem}: {_nf} フレーム "
          f"(τ<{DUMP_TAUC} 毎步 / 以降 {DUMP_STRIDE} 步ごと) "
          f"= {_nf * m * Nn * Neta * 16 / 1e6:.1f} MB")
    # =========================================================================
    _ENGINES.append(eng)
    return eng, p2, q2


ns.F.build_standard_universe = build_universe        # メモリ上のみ差し替え


def _rename(src: str, dst: str) -> None:
    """無改変 import の fig_one/fig_summary/fig_matrix が固定名で出す図を、
    モード名付きへ改名する（図の中身は同一関数の出力そのまま）。"""
    a, b = HERE / src, HERE / dst
    if a.exists():
        a.replace(b)


def arrays(eng):
    a = np.array(eng._rec, dtype=float)
    d = {k: a[:, i] for i, k in enumerate(REC_KEYS)}
    d["bands"] = np.array(eng._bands, dtype=float)   # (steps, Nn)
    d["ledger"] = np.array(eng._ledger, dtype=float)  # (snapshots, Nn, Nη)
    d["ledger_t"] = np.array(eng._ledger_t, dtype=float)
    d["targets"] = np.array(eng._tgt, dtype=float)    # (steps, #TARGETS)
    return d


def med_win(x):
    w = x[slice(*ns.WIN)]
    w = w[np.isfinite(w)]
    return float(np.median(w)) if len(w) else float("nan")


def first_over(x, thr):
    idx = np.flatnonzero(np.nan_to_num(x, nan=-np.inf) > thr)
    return int(idx[0]) + 1 if len(idx) else None


def bit_identical(A, B):
    worst = 0.0
    for k in ns.KEYS:
        x, y = A[k], B[k]
        if not np.all(np.isnan(x) == np.isnan(y)):
            return float("inf")
        d = np.abs(np.nan_to_num(x, nan=0.0) - np.nan_to_num(y, nan=0.0))
        worst = max(worst, float(d.max()) if d.size else 0.0)
    return worst


def fig_mix(n, hm, hv, rec):
    ts = np.arange(1, len(hm["r_mean"]) + 1)
    fig, ax = plt.subplots(4, 1, figsize=(9, 11), sharex=True)
    a = ax[0]
    a.plot(ts, hm["r_mean"], lw=0.9, color="tab:blue", label="r = sin²θ（node平均・物質）")
    a.fill_between(ts, hm["r_min"], hm["r_max"], color="tab:blue", alpha=0.15,
                   label="node の min–max")
    a.plot(ts, hm["r_nopump"], lw=1.0, color="tab:green",
           label="r_nopump（ポンプ帯を除く＝物質側）")
    a.plot(ts, hv["r_mean"], "k--", lw=0.8, label="真空対照")
    a.axhline(R_ALPHA, color="tab:red", lw=1.0, ls=":",
              label=f"α根 R_α=cos²(23π/124)={R_ALPHA:.9f}")
    a.set_ylabel("混合率 r（反射率）")
    a.set_title(f"N={n}  {CFG['label']}  r 窓中央値={rec['r_med_win']:.6f}  "
                f"α根との距離={rec['dist_alpha_win']:+.6f}")
    a.legend(fontsize=7)
    a = ax[1]
    a.semilogy(ts, np.maximum(np.abs(hm["dist_alpha"]), 1e-18), lw=0.9,
               color="tab:purple", label="|r − R_α|（系全体）")
    a.semilogy(ts, np.maximum(np.abs(hm["dist_alpha_nopump"]), 1e-18), lw=0.9,
               color="tab:green", label="|r_nopump − R_α|（物質側）")
    a.set_ylabel("α根からの距離"); a.legend(fontsize=7)
    a = ax[2]
    a.semilogy(ts, np.maximum(hm["odd_power"], 1e-40), lw=0.9, color="tab:red",
               label="全奇数帯（フェルミオン帯）")
    a.semilogy(ts, np.maximum(hm["even_power"], 1e-40), lw=0.9, color="tab:blue",
               label="全偶数帯（ボゾン帯・k≠0）")
    a.semilogy(ts, np.maximum(hv["odd_power"], 1e-40), "k--", lw=0.8,
               label="真空 奇数帯")
    a.set_ylabel("帯パワー"); a.legend(fontsize=7)
    a = ax[3]
    B = hm["bands"]
    for k in dict.fromkeys((PARTNER_K, 4, 5, 6, K_SEED, K_PUMP)):
        if k < B.shape[1]:
            if k == PARTNER_K:
                kind = "（相棒・奇）" if k % 2 else "（相棒・偶）"
            else:
                kind = "（偶）" if k % 2 == 0 else "（奇）"
            a.semilogy(ts, np.maximum(B[:, k], 1e-40), lw=0.8,
                       label=f"k={k}" + kind)
    a.set_ylabel("帯別パワー"); a.set_xlabel("τ（step）"); a.legend(fontsize=7, ncol=3)
    for b in ax:
        b.axvspan(ns.WIN[0], ns.WIN[1], color="green", alpha=0.06)
        if rec["tau_space"]:
            b.axvline(rec["tau_space"], color="tab:blue", lw=0.8, ls="-.")
    fig.tight_layout()
    fig.savefig(HERE / f"fig_{MODE}{TAG}_mix_N{n}_v2.png", dpi=110)
    plt.close(fig)


def fig_ledger(n, hm, rec):
    L = hm["ledger"]
    if L.size == 0:
        return
    Lm = np.median(L, axis=0)
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
    a = ax[0]
    im = a.imshow(np.log10(np.maximum(Lm.T, 1e-40)), aspect="auto", origin="lower",
                  cmap="viridis")
    a.set_xlabel("帯 k（偶=ボゾン・奇=フェルミオン）")
    a.set_ylabel("毛 η（巻き）")
    a.set_yticks(range(NETA))
    a.set_yticklabels([f"{e}({_signed(e):+d})" for e in range(NETA)], fontsize=7)
    a.set_xticks(range(NN)); a.set_xticklabels(range(NN), fontsize=7)
    for (k, e) in TARGETS:
        a.add_patch(plt.Rectangle((k - .5, e - .5), 1, 1, fill=False,
                                  edgecolor="red", lw=1.4))
    for (k, e, _d) in CELLS:
        a.add_patch(plt.Rectangle((k - .5, (e % NETA) - .5), 1, 1, fill=False,
                                  edgecolor="white", lw=1.0, ls="--"))
    a.add_patch(plt.Rectangle((K_PUMP - .5, (M_PUMP % NETA) - .5), 1, 1, fill=False,
                              edgecolor="orange", lw=1.6))
    fig.colorbar(im, ax=a, label="log₁₀ パワー（窓中央値）")
    a.set_title(f"N={n} 128セル帳簿（赤=狙った相棒・白破線=シード・橙=ポンプ）",
                fontsize=10)
    a = ax[1]
    ts = np.arange(1, hm["targets"].shape[0] + 1)
    for i, (k, e) in enumerate(TARGETS):
        a.semilogy(ts, np.maximum(hm["targets"][:, i], 1e-40), lw=0.9,
                   label=f"k={k} m={_signed(e):+d} {SPECIES.get(_signed(e),'?')}")
    a.semilogy(ts, np.maximum(hm["nontarget_power"], 1e-40), "k--", lw=1.0,
               label="狙っていないセルの総和")
    a.axvspan(ns.WIN[0], ns.WIN[1], color="green", alpha=0.06)
    a.set_xlabel("τ（step）"); a.set_ylabel("パワー")
    a.set_title(f"狙った相棒セルの成長と非狙いセル（排他比中央値 "
                f"{rec['excl_ratio_med']:.2e}）", fontsize=10)
    a.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_{MODE}{TAG}_ledger_N{n}_v2.png", dpi=120)
    plt.close(fig)


def _assert_replicate_paths_unused(nmin: int, nmax: int) -> None:
    """suffix 付き複製走行で、既存の同名成果を上書きしない。

    import 元の図関数は一時的に固定名 `fig_nsweep_*_v1.png` を使うため、
    その残骸がある場合も開始しない。複数プロセスの並列走行は不可。
    """
    if not OUTPUT_SUFFIX:
        return
    finals = [
        HERE / f"result_nsweep_{MODE}{TAG}_v2.json",
        HERE / f"fig_{MODE}{TAG}_summary_v2.png",
        HERE / f"fig_{MODE}{TAG}_birth_matrix_v2.png",
    ]
    fixed = [HERE / "fig_nsweep_summary_v1.png",
             HERE / "fig_nsweep_birth_matrix_v1.png"]
    for n in range(nmin, nmax + 1):
        finals.extend([
            HERE / f"nsweep_{MODE}{TAG}_N{n}_v2.npz",
            HERE / f"fig_{MODE}{TAG}_4panel_N{n}_v2.png",
            HERE / f"fig_{MODE}{TAG}_mix_N{n}_v2.png",
            HERE / f"fig_{MODE}{TAG}_ledger_N{n}_v2.png",
        ])
        fixed.append(HERE / f"fig_nsweep_N{n}_v1.png")
    occupied = [p.name for p in finals + fixed if p.exists()]
    if occupied:
        raise FileExistsError(
            "複製走行の出力先または固定名一時図が既に存在: "
            + ", ".join(occupied)
        )


def main():
    t0 = time.time()
    nmin = int(sys.argv[2]) if len(sys.argv) > 3 else 1
    nmax = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    _assert_replicate_paths_unused(nmin, nmax)
    print(f"=== N 掃引 {nmin}→{nmax}・mode={MODE}（{CFG['label']}）"
          f"・F=v1・Nn={NN}・Nη={NETA}・δ={DELTA}・m*={M_STAR}"
          f"（毛idx {M_STAR_IDX}）・相棒帯 k={PARTNER_K}・T={ns.T} ===")
    if OUTPUT_SUFFIX:
        print(f"    複製走行 output_suffix={OUTPUT_SUFFIX!r} / tag={OUTPUT_TAG}")
    print(f"    α根 R_α = cos²(23π/124) = {R_ALPHA:.9f} / 1−R_α = {T_ALPHA:.9f} / "
          f"4π/(1−R_α)² = {ALPHA_INV:.6f}（α⁻¹実測 137.035999206）")
    print(f"    奇数帯 k={ODD_K} / 偶数帯 k={EVEN_K}")
    recs, fails = [], []
    out = {"env": {"mode": MODE, "label": CFG["label"], "Nn": NN, "Neta": NETA,
                   "T": ns.T, "delta": DELTA, "seed": ns.SEED,
                   "cell": list(ns.CELL), "order": ns.ORDER,
                   "window": list(ns.WIN),
                   "output_suffix": OUTPUT_SUFFIX or None,
                   "seed_cells": CELL_RECORDS,
                   "target_cells": TARGET_RECORDS,
                   "seed_strength": {
                       "n_cells": len(CELLS),
                       "n_delta": N_DELTA,
                       "n_delta2": N_DELTA2,
                       "nF": NF,
                       "nB": NB,
                       "PF": PF_DECLARED,
                       "PB": PB_DECLARED,
                       "Pseed": PSEED_DECLARED,
                       "Acoh": N_DELTA,
                       "r_np_cell_ratio": R_NP_CELL_RATIO,
                       "r_np_power_ratio": R_NP_POWER_RATIO,
                       "definitions": {
                           "n_delta": "Σ_j δ_j（一様強度で nδ）",
                           "n_delta2": "Σ_j δ_j^2（一様強度で nδ^2）",
                           "PF_PB": "奇数帯/偶数帯シードの Σ_j δ_j^2（宣言値）",
                           "Pseed": "PF+PB（宣言シードパワー）",
                           "Acoh": "Σ_j δ_j（一様強度で nδ）",
                           "r_np_cell_ratio": "nF/(nF+nB)",
                           "r_np_power_ratio": "PF/(PF+PB)",
                       },
                   },
                   "primary": {
                       "pump": {"k": K_PUMP, "m_signed": M_PUMP,
                                "eta_index": M_PUMP % NETA},
                       "seed": {"k": K_SEED, "m_signed": M_SEED,
                                "eta_index": M_SEED % NETA,
                                "declared_only": not bool(CELLS)},
                       "partner": {"k": PARTNER_K, "m_signed": M_STAR,
                                   "eta_index": M_STAR_IDX},
                   },
                   "legacy_diagnostics": {
                       "seed_closure_cell": [1, 0],
                       "seed_closure_note": (
                           "import元 run_one の固定診断。boson_family の実シード"
                           " k=6 の閉塞値ではなく、同 mode では適用外"
                       ),
                       "conc_k3_and_k3_power_note": (
                           "互換キー名。読出し帯は primary.partner.k"
                       ),
                       "matter_born_note": (
                           "import元 summarize の legacy matter_born は "
                           "g_matter_fraction=f_seed（奇数帯パワー分率）>1e-30 "
                           "の判定。ボゾン積荷の存在判定ではない"
                       ),
                   },
                   "functions": ["unified_interaction_v1（★v1・恒久ルール）",
                                 "unified_dimension_v1", "unified_readout_v3",
                                 "selection_v1"],
                   "recipe": {"k_pump": K_PUMP, "k_seed": K_SEED,
                              "partner_k": PARTNER_K,
                              "m_pump": M_PUMP, "m_seed": M_SEED,
                              "m_star": M_STAR, "hair_index": M_STAR_IDX,
                              "sum_rule": "m*=2m_pump−m_seed / k*=2k_pump−k_seed"},
                   "alpha": {"theta": "23π/124", "R_alpha": R_ALPHA,
                             "T_alpha": T_ALPHA, "alpha_inv_from_R": ALPHA_INV,
                             "alpha_inv_measured": 137.035999206},
                   "odd_k": ODD_K, "even_k": EVEN_K,
                   "base_script_md5": "bfa5d854b637fe97c33a7148be9c7f86"},
           "N": {}, "failed": {}}
    for n in range(nmin, nmax + 1):
        t1 = time.time()
        _ENGINES.clear()
        try:
            Hm, Rm, Am, Ccm, Csm = ns.run_one(n, DELTA)
            # ダンプ版 v1 の追加（3/3）：全τの瞬時フレーム束を書き出す
            np.savez(HERE / f"dump_frame_{MODE}{TAG}_N{n}_m_v1.npz",
                     **ns.FRAME_DUMP[-1])
            eng_m = _ENGINES[-1]
            Hv, Rv, Av, Ccv, Csv = ns.run_one(n, 0.0)
            np.savez(HERE / f"dump_frame_{MODE}{TAG}_N{n}_v_v1.npz",
                     **ns.FRAME_DUMP[-1])
            eng_v = _ENGINES[-1]
        except Exception as ex:
            msg = f"{type(ex).__name__}: {ex}"
            fails.append(n); out["failed"][n] = msg
            out["N"][n] = {"N": n, "M": n * (n - 1) // 2, "built": False,
                           "error": msg}
            print(f"N={n:3d} M={n*(n-1)//2:4d}: **構成不能** {msg[:70]}")
            (HERE / f"result_nsweep_{MODE}{TAG}_v2.json").write_text(
                json.dumps(out, indent=1, ensure_ascii=False, default=float))
            continue
        hm, hv = arrays(eng_m), arrays(eng_v)
        assert len(hm["r_mean"]) == ns.T, f"記録数 {len(hm['r_mean'])} ≠ T"
        rec = ns.summarize(n, Hm, Rm, Am, Ccm, Csm, Hv, Av)
        B = hm["bands"]
        even_rise = {str(k): first_over(B[:, k], 0.0) for k in EVEN_K if k != K_PUMP}
        odd_rise = {str(k): first_over(B[:, k], 0.0) for k in ODD_K if k != K_SEED}
        rec.update({
            "r_med_win": med_win(hm["r_mean"]),
            "r_nopump_med_win": med_win(hm["r_nopump"]),
            "dist_alpha_nopump_win": med_win(hm["dist_alpha_nopump"]),
            "r_nopump_final": float(hm["r_nopump"][-1]),
            "even_power_nopump_med": med_win(hm["even_power_nopump"]),
            "r_raw_med_win": med_win(hm["r_raw"]),
            "r_min_win": med_win(hm["r_min"]), "r_max_win": med_win(hm["r_max"]),
            "dist_alpha_win": med_win(hm["dist_alpha"]),
            "absdist_alpha_win": med_win(hm["absdist_alpha"]),
            "r_final": float(hm["r_mean"][-1]),
            "dist_alpha_final": float(hm["dist_alpha"][-1]),
            "r_med_win_vacuum": med_win(hv["r_mean"]),
            "conc_k3_med": med_win(hm["conc_k3"]),
            "dom_m_med": med_win(hm["dom_m"]), "q_hat_med": med_win(hm["q_hat"]),
            "readable_rate": med_win(hm["readable"]),
            "k3_power_med": med_win(hm["k3_power"]),
            # primary 相棒帯の汎用名。上の既存キーと同じ配列の alias。
            "conc_partner_med": med_win(hm["conc_k3"]),
            "partner_dom_m_med": med_win(hm["dom_m"]),
            "partner_q_hat_med": med_win(hm["q_hat"]),
            "partner_readable_rate": med_win(hm["readable"]),
            "partner_power_med": med_win(hm["k3_power"]),
            "odd_power_max": float(np.nanmax(hm["odd_power"])),
            "odd_amp_max": float(np.nanmax(hm["odd_amp_max"])),
            "odd_exact_zero_steps": int((hm["odd_power"] == 0.0).sum()),
            "even_power_med": med_win(hm["even_power"]),
            "odd_power_med": med_win(hm["odd_power"]),
            "band_rise_even": even_rise, "band_rise_odd": odd_rise,
            "band_power_med": {str(k): med_win(B[:, k]) for k in range(NN)},
            "target_power_med": med_win(hm["target_power"]),
            "nontarget_power_med": med_win(hm["nontarget_power"]),
            "excl_ratio_med": med_win(hm["excl_ratio"]),
            "payload_seed_power_med": med_win(hm["seed_power"]),
            "payload_seed_power_max": float(np.nanmax(hm["seed_power"])),
            "primary_seed_power_med": med_win(hm["primary_seed_power"]),
            "primary_seed_power_max": float(np.nanmax(hm["primary_seed_power"])),
            "pump_power_med": med_win(hm["pump_power"]),
            "pump_power_max": float(np.nanmax(hm["pump_power"])),
            "target_exact_zero_steps": int((hm["target_power"] == 0.0).sum()),
            # R は非負。node 最大値が厳密 0 なら全 node で R≡0。
            "readout_R_exact_zero_steps": int((hm["r_max"] == 0.0).sum()),
            "targets": [[int(k), int(_signed(e)), SPECIES.get(_signed(e), "?")]
                        for (k, e) in TARGETS],
            "target_power_each_med": {
                f"k{k}_m{_signed(e):+d}": med_win(hm["targets"][:, i])
                for i, (k, e) in enumerate(TARGETS)},
            "ledger_med": [[float(np.median(hm["ledger"][:, k, e]))
                            for e in range(NETA)] for k in range(NN)],
            "determinism_max_abs": bit_identical(Hm, Hv) if DELTA == 0 else None,
        })
        recs.append(rec); out["N"][n] = rec
        ns.fig_one(n, Hm, Hv, Am, Ccm, Csm, rec)
        _rename(f"fig_nsweep_N{n}_v1.png", f"fig_{MODE}{TAG}_4panel_N{n}_v2.png")
        fig_mix(n, hm, hv, rec)
        fig_ledger(n, hm, rec)
        np.savez_compressed(
            HERE / f"nsweep_{MODE}{TAG}_N{n}_v2.npz",
            **{f"m_{k}": Hm[k] for k in ns.KEYS},
            **{f"v_{k}": Hv[k] for k in ns.KEYS},
            m_resid=Rm, m_acq=Am, v_acq=Av,
            m_cond_closure=Ccm, m_seed_closure=Csm,
            **{f"rec_m_{k}": v for k, v in hm.items()},
            **{f"rec_v_{k}": v for k, v in hv.items()},
            # 互換キー rec_*_conc_k3 / rec_*_k3_power を残したままの alias。
            rec_m_conc_partner=hm["conc_k3"],
            rec_m_partner_dom_m=hm["dom_m"],
            rec_m_partner_q_hat=hm["q_hat"],
            rec_m_partner_readable=hm["readable"],
            rec_m_partner_power=hm["k3_power"],
            rec_v_conc_partner=hv["conc_k3"],
            rec_v_partner_dom_m=hv["dom_m"],
            rec_v_partner_q_hat=hv["q_hat"],
            rec_v_partner_readable=hv["readable"],
            rec_v_partner_power=hv["k3_power"],
            seed_cells_index=np.array(
                [(k, e % NETA) for (k, e, _d) in CELLS], dtype=int).reshape(-1, 2),
            seed_cells_delta=np.array([d for (_k, _e, d) in CELLS], dtype=float),
            targets_index_2d=np.array(TARGETS, dtype=int).reshape(-1, 2),
            # legacy キーは vacuum 時の shape=(0,) を含め、既存配列と同一に保つ。
            targets_index=np.array(TARGETS, dtype=float))
        print(f"N={n:3d} M={rec['M']:4d}: 空間τ={str(rec['tau_space']):>5} "
              f"物質={str(rec['matter_born']):>5} 時間τ={str(rec['tau_time']):>5} "
              f"| r={rec['r_med_win']:.6f} r_np={rec['r_nopump_med_win']:.6f} "
              f"Δα_np={rec['dist_alpha_nopump_win']:+.3e} "
              f"奇={rec['odd_power_med']:.2e} 偶={rec['even_power_med']:.2e} "
              f"集中度={rec['conc_k3_med']:.4f} m̂={rec['dom_m_med']:+.1f} "
              f"狙い={rec['target_power_med']:.2e} 排他比={rec['excl_ratio_med']:.2e} "
              f"[{time.time()-t1:.0f}s]")
        (HERE / f"result_nsweep_{MODE}{TAG}_v2.json").write_text(
            json.dumps(out, indent=1, ensure_ascii=False, default=float))
    if recs:
        ns.fig_summary(recs, fails)
        _rename("fig_nsweep_summary_v1.png", f"fig_{MODE}{TAG}_summary_v2.png")
    ns.fig_matrix(recs, fails, nmin, nmax)
    _rename("fig_nsweep_birth_matrix_v1.png", f"fig_{MODE}{TAG}_birth_matrix_v2.png")
    out["failed_N"] = fails
    # V4 は neutral/electron の単一種巻き集中度試験。多種・無シード・
    # family に適用すると、NaN 空集合が all(empty)=True になる場合がある。
    v4_applicable = MODE in {"neutral", "electron"}
    if v4_applicable:
        v4_conc_ge_09 = bool(recs) and all(
            (r["conc_k3_med"] >= 0.9)
            for r in recs if np.isfinite(r["conc_k3_med"])
        )
    else:
        v4_conc_ge_09 = None
    b3_applicable = MODE == "boson_family"
    out["judgments"] = {
        "V1_failed_N": fails,
        "V2_space_born_all_N": bool(recs) and all(r["space_born"] for r in recs),
        "V3_matter_born_all_N": bool(recs) and all(r["matter_born"] for r in recs),
        "V3_odd_exact_zero_all_N": bool(recs) and all(
            r["odd_exact_zero_steps"] == ns.T for r in recs),
        "V4_conc_ge_09_all_N": v4_conc_ge_09,
        "V4_applicable": v4_applicable,
        "V4_note": (
            "neutral/electron の単一種 primary 相棒巻き集中度のみに適用。"
            "vacuum/mixed/family ablation は N/A"
        ),
        "B3_applicable": b3_applicable,
        "B3_note": (
            "boson_family で payload は存在するが、奇数帯 R と"
            " primary 相棒生成が厳密零かを分離する判定。他 mode は N/A"
        ),
        "B3_payload_seed_present_all_N": (
            bool(recs) and all(r["payload_seed_power_max"] > 0.0 for r in recs)
            if b3_applicable else None),
        "B3_target_exact_zero_all_steps_all_N": (
            bool(recs) and all(r["target_exact_zero_steps"] == ns.T for r in recs)
            if b3_applicable else None),
        "B3_readout_R_exact_zero_all_steps_all_N": (
            bool(recs) and all(r["readout_R_exact_zero_steps"] == ns.T for r in recs)
            if b3_applicable else None),
        "V5_r_med_range": [min((r["r_med_win"] for r in recs), default=None),
                           max((r["r_med_win"] for r in recs), default=None)],
        "V5_absdist_alpha_min": min((r["absdist_alpha_win"] for r in recs),
                                    default=None),
        "V5_r_nopump_range": [min((r["r_nopump_med_win"] for r in recs), default=None),
                              max((r["r_nopump_med_win"] for r in recs), default=None)],
        "V5_absdist_alpha_nopump_min": min(
            (abs(r["dist_alpha_nopump_win"]) for r in recs
             if r["dist_alpha_nopump_win"] == r["dist_alpha_nopump_win"]),
            default=None),
        "V7_excl_ratio_med_range": [
            min((r["excl_ratio_med"] for r in recs), default=None),
            max((r["excl_ratio_med"] for r in recs), default=None)],
        "V7_targets": recs[0]["targets"] if recs else None,
        "V6_even_bands_risen": {str(k): sum(
            1 for r in recs if r["band_rise_even"].get(str(k)) is not None)
            for k in EVEN_K if k != K_PUMP},
    }
    out["runtime_sec"] = time.time() - t0
    (HERE / f"result_nsweep_{MODE}{TAG}_v2.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"\n構成できなかった N: {fails if fails else 'なし'}")
    print(f"判定: {json.dumps(out['judgments'], ensure_ascii=False, default=float)}")
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
