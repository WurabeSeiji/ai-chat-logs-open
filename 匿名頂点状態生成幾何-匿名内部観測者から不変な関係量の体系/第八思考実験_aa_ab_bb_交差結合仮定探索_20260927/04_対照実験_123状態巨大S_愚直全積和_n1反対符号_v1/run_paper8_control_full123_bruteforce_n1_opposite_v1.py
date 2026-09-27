#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper 8 control experiment -- strict full-123 brute-force implementation, v1.

PURPOSE
=======
Run exactly one condition:
    n = 1, opposite-sign pair
    aa: (C,D) = (0.91, 0.00)
    ab: (C,D) = (1.09, 0.36)
    bb: (C,D) = (0.91, 0.00)

This program intentionally prioritizes auditability over speed.

NON-NEGOTIABLE IMPLEMENTATION RULES
====================================
1. One persistent state vector only:
       psi.shape == (123,)
   It is allocated once and used for the whole run. No aa/ab/bb state slices,
   no 41-state copies, no concatenate/reshape-based local transitions.

2. One persistent full interaction matrix only:
       S.shape == (123,123)
   It is allocated once and reused for the whole run. No 41x41 local matrices,
   no block copies, no sparse matrix, no list of local actions.

3. Every microstep reconstructs ALL 15,129 entries of the same S array.
   Structural zeros are explicitly assigned as zeros. Entries are never skipped
   because a state component or q component happens to be zero.

4. All 11 one-hot candidate phases are represented simultaneously as q-column
   coefficients. There is NO "active phase = argmax(q)" branch in state
   generation. A q value of zero makes its candidate contribution zero only in
   the final product-sum; the candidate coefficient still exists in S.

5. The state update is one full 123x123 by 123 product-sum. The strict engine
   explicitly visits every (i,j) pair and performs
       acc += S[i,j] * psi[j]
   with no zero-value skip and no sparse optimization.

6. One macrostep = 11 strict microsteps. Microstep states are NOT persisted.
   The initial full 123-state and EVERY completed macrostep full 123-state are
   persisted with no thinning and no omission.

7. The generator never reads Paper 6 / Paper 7 reference trajectories. Reference
   comparison must be a separate post-run program.

8. Readout / progress / storage never feed back into psi or S.

WARNING
=======
This implementation is intentionally very slow. That is expected. Do not add
state-value pruning, q-phase pruning, sparse products, 41-state sub-transitions,
or reference-driven shortcuts without changing the experiment identity.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

# -----------------------------------------------------------------------------
# Fixed experiment constants inherited from Paper 6 / Paper 7 state machine.
# -----------------------------------------------------------------------------
STEPS_PER_ORBIT = 4000
HSTEP = 2.0 * math.pi / STEPS_PER_ORBIT
PH_RE = math.cos(HSTEP)
PH_IM = math.sin(HSTEP)
TAU0 = 1.0e4
P0 = 50.0
TARGET_P = 20.0
N0 = 0.25

# Local 41-state layout, used only as numeric index constants into global psi.
UR, UI, P, E, HR, HI, Q, N, C, D = range(10)
K1 = 10
K2 = 13
K3 = 16
K4 = 19
PS = 22
ES = 23
DV = 24
PM = 27
EM = 28
DPHI = 29
QB = 30
NST = 41
NCH = 3
NTOT = 123
NQ = 11
MICROSTEPS_PER_MACRO = 11
PRODUCTS_PER_MICRO = NTOT * NTOT
PRODUCTS_PER_MACRO = MICROSTEPS_PER_MACRO * PRODUCTS_PER_MICRO
CHANNEL_NAMES = ("aa", "ab", "bb")
CHANNEL_BASES = (0, NST, 2 * NST)
CHANNEL_CD = ((0.91, 0.0), (1.09, 0.36), (0.91, 0.0))

# Raw persistence policy.
ROWS_PER_PART = 100_000
HDF5_CHUNK_ROWS = 1024
PROGRESS_EVERY_MACROS = 100
PROGRESS_EVERY_SECONDS = 30.0


# -----------------------------------------------------------------------------
# Small scalar helpers. They never accept a 41-state subvector.
# -----------------------------------------------------------------------------
def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _require_rate_domain(p: float, c: float, label: str) -> None:
    if not (math.isfinite(p) and math.isfinite(c) and p > 0.0 and c > 0.0):
        raise ValueError(f"rates domain failure at {label}: p={p!r}, c={c!r}")


def rate_dp(p: float, n: float, c: float, d: float) -> float:
    _require_rate_domain(p, c, "rate_dp")
    rp = math.sqrt(p)
    rc = math.sqrt(c)
    return -(4.0 / 3.0) * n * d * rc / rp - (64.0 / 5.0) * n * c * rc / (p * rp)


def rate_de(_p: float, _n: float, _c: float, _d: float) -> float:
    return 0.0


def rate_dt(p: float, _n: float, c: float, _d: float) -> float:
    _require_rate_domain(p, c, "rate_dt")
    return p * math.sqrt(p) / math.sqrt(c)


def gidx(base: int, local_index: int) -> int:
    """Global index helper; returns an integer only, never a state slice."""
    return base + local_index


# -----------------------------------------------------------------------------
# Initialization: one 123-state vector, allocated once.
# -----------------------------------------------------------------------------
def initialize_psi(psi: np.ndarray) -> None:
    if psi.shape != (NTOT,):
        raise ValueError("psi must be shape (123,)")
    # Explicitly touch every state component once at initialization.
    for i in range(NTOT):
        psi[i] = 0.0

    for ch in range(NCH):
        base = CHANNEL_BASES[ch]
        c0, d0 = CHANNEL_CD[ch]
        psi[gidx(base, UR)] = 1.0
        psi[gidx(base, UI)] = 0.0
        psi[gidx(base, P)] = P0
        psi[gidx(base, E)] = 0.0
        psi[gidx(base, HR)] = 1.0
        psi[gidx(base, HI)] = 0.0
        psi[gidx(base, Q)] = 1.0
        psi[gidx(base, N)] = N0
        psi[gidx(base, C)] = c0
        psi[gidx(base, D)] = d0
        psi[gidx(base, PS)] = P0
        psi[gidx(base, PM)] = P0
        psi[gidx(base, QB)] = 1.0  # q0=1, q1..q10=0


# -----------------------------------------------------------------------------
# Candidate value F_{row,phase}(psi), expressed only with absolute global reads.
# This function is called for every local row 0..29 and every phase 0..10 when
# the corresponding q-column coefficient of S is constructed. It never checks
# whether q_phase is zero.
# -----------------------------------------------------------------------------
def candidate_value(psi: np.ndarray, base: int, row: int, phase: int) -> float:
    # Default candidate: OLD value survives this phase.
    value = float(psi[gidx(base, row)])

    p = float(psi[gidx(base, P)])
    e = float(psi[gidx(base, E)])
    n = float(psi[gidx(base, N)])
    c = float(psi[gidx(base, C)])
    d = float(psi[gidx(base, D)])
    ps = float(psi[gidx(base, PS)])
    dv0 = float(psi[gidx(base, DV)])
    dv1 = float(psi[gidx(base, DV + 1)])
    dv2 = float(psi[gidx(base, DV + 2)])
    dphi = float(psi[gidx(base, DPHI)])

    # phase 0: k1 = rates(P,N,C,D)
    if phase == 0:
        if row == K1:
            value = rate_dp(p, n, c, d)
        elif row == K1 + 1:
            value = rate_de(p, n, c, d)
        elif row == K1 + 2:
            value = rate_dt(p, n, c, d)

    # phase 1: midpoint from k1
    elif phase == 1:
        if row == PS:
            value = p + 0.5 * HSTEP * float(psi[gidx(base, K1)])
        elif row == ES:
            value = e + 0.5 * HSTEP * float(psi[gidx(base, K1 + 1)])

    # phase 2: k2 = rates(PS,N,C,D)
    elif phase == 2:
        if row == K2:
            value = rate_dp(ps, n, c, d)
        elif row == K2 + 1:
            value = rate_de(ps, n, c, d)
        elif row == K2 + 2:
            value = rate_dt(ps, n, c, d)

    # phase 3: midpoint from k2
    elif phase == 3:
        if row == PS:
            value = p + 0.5 * HSTEP * float(psi[gidx(base, K2)])
        elif row == ES:
            value = e + 0.5 * HSTEP * float(psi[gidx(base, K2 + 1)])

    # phase 4: k3 = rates(PS,N,C,D)
    elif phase == 4:
        if row == K3:
            value = rate_dp(ps, n, c, d)
        elif row == K3 + 1:
            value = rate_de(ps, n, c, d)
        elif row == K3 + 2:
            value = rate_dt(ps, n, c, d)

    # phase 5: endpoint predictor from k3
    elif phase == 5:
        if row == PS:
            value = p + HSTEP * float(psi[gidx(base, K3)])
        elif row == ES:
            value = e + HSTEP * float(psi[gidx(base, K3 + 1)])

    # phase 6: k4 = rates(PS,N,C,D)
    elif phase == 6:
        if row == K4:
            value = rate_dp(ps, n, c, d)
        elif row == K4 + 1:
            value = rate_de(ps, n, c, d)
        elif row == K4 + 2:
            value = rate_dt(ps, n, c, d)

    # phase 7: RK4 combined increment
    elif phase == 7:
        if DV <= row <= DV + 2:
            a = row - DV
            k1a = float(psi[gidx(base, K1 + a)])
            k2a = float(psi[gidx(base, K2 + a)])
            k3a = float(psi[gidx(base, K3 + a)])
            k4a = float(psi[gidx(base, K4 + a)])
            value = (HSTEP / 6.0) * (k1a + 2.0 * k2a + 2.0 * k3a + k4a)

    # phase 8: midpoint readout state
    elif phase == 8:
        if row == PM:
            value = p + 0.5 * dv0
        elif row == EM:
            value = e + 0.5 * dv1

    # phase 9: phase increment
    elif phase == 9:
        if row == DPHI:
            value = HSTEP

    # phase 10: macro-boundary update and work-register reset
    elif phase == 10:
        if K1 <= row <= K1 + 2:
            value = 0.0
        elif K2 <= row <= K2 + 2:
            value = 0.0
        elif K3 <= row <= K3 + 2:
            value = 0.0
        elif K4 <= row <= K4 + 2:
            value = 0.0
        elif DV <= row <= DV + 2:
            value = 0.0
        elif row == PM:
            value = p + dv0
        elif row == EM:
            value = e + dv1
        elif row == DPHI:
            value = 0.0
        elif row == UR:
            ur = float(psi[gidx(base, UR)])
            ui = float(psi[gidx(base, UI)])
            value = ur * PH_RE - ui * PH_IM
        elif row == UI:
            ur = float(psi[gidx(base, UR)])
            ui = float(psi[gidx(base, UI)])
            value = ur * PH_IM + ui * PH_RE
        elif row == P:
            value = p + dv0
        elif row == E:
            value = e + dv1
        elif row == HR:
            hr = float(psi[gidx(base, HR)])
            hi = float(psi[gidx(base, HI)])
            cd = math.cos(dphi)
            sd = math.sin(dphi)
            value = hr * cd - hi * sd
        elif row == HI:
            hr = float(psi[gidx(base, HR)])
            hi = float(psi[gidx(base, HI)])
            cd = math.cos(dphi)
            sd = math.sin(dphi)
            value = hr * sd + hi * cd
        elif row == Q:
            q_old = float(psi[gidx(base, Q)])
            value = q_old * math.exp(dv2 / TAU0)
        elif row == PS:
            value = p + dv0
        elif row == ES:
            value = e + dv1

    else:
        raise ValueError(f"invalid phase {phase}")

    return float(value)


# -----------------------------------------------------------------------------
# Every S[i,j] is explicitly assigned on every microstep.
# Structural zero is determined by index geometry only, never by psi[j]==0 or
# q_phase==0. No state-value pruning exists here.
# -----------------------------------------------------------------------------
def interaction_entry(psi: np.ndarray, out_i: int, in_j: int) -> float:
    out_ch = out_i // NST
    in_ch = in_j // NST

    # Control experiment: all cross-channel blocks are structurally zero.
    if out_ch != in_ch:
        return 0.0

    base = CHANNEL_BASES[out_ch]
    li = out_i - base
    lj = in_j - base

    # Logical/work rows 0..29 are the full 11-candidate one-hot product-sum.
    # Every q-column exists in S regardless of whether that q value is 0 or 1.
    if li < QB:
        if QB <= lj < QB + NQ:
            phase = lj - QB
            return candidate_value(psi, base, li, phase)
        return 0.0

    # q cyclic shift rows 30..40.
    if li == QB:
        return 1.0 if lj == QB + 10 else 0.0
    return 1.0 if lj == li - 1 else 0.0


def rebuild_full_S(psi: np.ndarray, S: np.ndarray) -> None:
    """Rebuild the same persistent 123x123 S object, touching all 15,129 cells."""
    if S.shape != (NTOT, NTOT):
        raise ValueError("S must be shape (123,123)")
    for i in range(NTOT):
        for j in range(NTOT):
            S[i, j] = interaction_entry(psi, i, j)


# -----------------------------------------------------------------------------
# Strict full product-sum. No if on S[i,j], psi[j], or q. Every pair is visited.
# -----------------------------------------------------------------------------
def full_product_sum(S: np.ndarray, psi: np.ndarray, psi_next: np.ndarray) -> None:
    if S.shape != (NTOT, NTOT) or psi.shape != (NTOT,) or psi_next.shape != (NTOT,):
        raise ValueError("shape mismatch in full_product_sum")
    for i in range(NTOT):
        acc = 0.0
        for j in range(NTOT):
            acc += float(S[i, j]) * float(psi[j])
        psi_next[i] = acc


# -----------------------------------------------------------------------------
# State validation and macro stop condition. These are readout/control-flow only;
# they never modify psi or S.
# -----------------------------------------------------------------------------
def q_sum(psi: np.ndarray, base: int) -> float:
    s = 0.0
    for k in range(NQ):
        s += float(psi[gidx(base, QB + k)])
    return s


def q_nonzero_count(psi: np.ndarray, base: int) -> int:
    count = 0
    for k in range(NQ):
        if float(psi[gidx(base, QB + k)]) != 0.0:
            count += 1
    return count


def validate_macro_boundary(psi: np.ndarray) -> None:
    for i in range(NTOT):
        if not math.isfinite(float(psi[i])):
            raise FloatingPointError(f"nonfinite psi[{i}]={psi[i]!r}")
    for ch in range(NCH):
        base = CHANNEL_BASES[ch]
        p = float(psi[gidx(base, P)])
        c = float(psi[gidx(base, C)])
        ps = float(psi[gidx(base, PS)])
        if p <= 0.0 or c <= 0.0 or ps <= 0.0:
            raise ValueError(f"domain failure channel={CHANNEL_NAMES[ch]} P={p} C={c} PS={ps}")
        qs = q_sum(psi, base)
        qnz = q_nonzero_count(psi, base)
        if qs != 1.0 or qnz != 1:
            raise ValueError(f"q one-hot failure channel={CHANNEL_NAMES[ch]} qsum={qs} nonzero={qnz}")


def all_channels_reached_target(psi: np.ndarray) -> bool:
    for ch in range(NCH):
        base = CHANNEL_BASES[ch]
        if float(psi[gidx(base, P)]) > TARGET_P:
            return False
    return True


# -----------------------------------------------------------------------------
# Full macrostate persistence. Initial row + EVERY completed macrostep row.
# -----------------------------------------------------------------------------
class MacroStateWriter:
    def __init__(self, out_dir: Path, rows_per_part: int = ROWS_PER_PART):
        self.out_dir = out_dir
        self.rows_per_part = rows_per_part
        self.part_index = -1
        self.file = None
        self.ds_step = None
        self.ds_psi = None
        self.rows_in_part = 0
        self.total_rows = 0
        self.parts = []

    def _open_next_part(self, start_step: int) -> None:
        self.close_part()
        self.part_index += 1
        path = self.out_dir / f"raw_macro_part{self.part_index:03d}.h5"
        f = h5py.File(path, "w")
        f.attrs["format"] = "paper8_full123_macrostate_v1"
        f.attrs["state_dimension"] = NTOT
        f.attrs["condition"] = "n1_opposite_sign"
        f.attrs["channels"] = json.dumps(CHANNEL_NAMES)
        f.attrs["channel_CD"] = json.dumps(CHANNEL_CD)
        f.attrs["start_macro_step"] = int(start_step)
        f.attrs["created_utc"] = utc_now()
        self.ds_step = f.create_dataset(
            "macro_step", shape=(0,), maxshape=(None,), dtype="<i8",
            chunks=(HDF5_CHUNK_ROWS,)
        )
        self.ds_psi = f.create_dataset(
            "psi123", shape=(0, NTOT), maxshape=(None, NTOT), dtype="<f8",
            chunks=(HDF5_CHUNK_ROWS, NTOT)
        )
        self.file = f
        self.rows_in_part = 0
        self.parts.append({
            "part_index": self.part_index,
            "file": path.name,
            "start_macro_step": int(start_step),
            "end_macro_step": None,
            "rows": 0,
        })

    def append(self, macro_step: int, psi: np.ndarray) -> None:
        if self.file is None or self.rows_in_part >= self.rows_per_part:
            self._open_next_part(macro_step)
        r = self.rows_in_part
        self.ds_step.resize((r + 1,))
        self.ds_psi.resize((r + 1, NTOT))
        self.ds_step[r] = np.int64(macro_step)
        # Persist the complete 123-state row, with no thinning.
        self.ds_psi[r, :] = psi
        self.rows_in_part += 1
        self.total_rows += 1
        p = self.parts[-1]
        p["end_macro_step"] = int(macro_step)
        p["rows"] = int(self.rows_in_part)

    def flush(self) -> None:
        if self.file is not None:
            self.file.flush()

    def close_part(self) -> None:
        if self.file is not None:
            self.file.flush()
            self.file.close()
            self.file = None
            self.ds_step = None
            self.ds_psi = None

    def close(self) -> None:
        self.close_part()


# -----------------------------------------------------------------------------
# Progress reporting: stdout + atomic progress.json. No state feedback.
# -----------------------------------------------------------------------------
def physical_progress_fraction(psi: np.ndarray) -> float:
    slowest_p = -math.inf
    for ch in range(NCH):
        p = float(psi[gidx(CHANNEL_BASES[ch], P)])
        if p > slowest_p:
            slowest_p = p
    frac = (P0 - slowest_p) / (P0 - TARGET_P)
    return max(0.0, min(1.0, frac))


def progress_payload(
    *, status: str, macro_completed: int, micro_in_macro: int,
    psi: np.ndarray, writer: MacroStateWriter, started_monotonic: float,
    last_message: str = ""
) -> dict:
    elapsed = max(0.0, time.monotonic() - started_monotonic)
    rate = (macro_completed / elapsed) if elapsed > 0.0 else 0.0
    frac = physical_progress_fraction(psi)
    eta_seconds = (elapsed * (1.0 - frac) / frac) if frac > 0.0 and frac < 1.0 else None
    ps = {}
    cds = {}
    qinfo = {}
    for ch in range(NCH):
        name = CHANNEL_NAMES[ch]
        base = CHANNEL_BASES[ch]
        ps[name] = float(psi[gidx(base, P)])
        cds[name] = {
            "C": float(psi[gidx(base, C)]),
            "D": float(psi[gidx(base, D)]),
        }
        qinfo[name] = {
            "sum": q_sum(psi, base),
            "nonzero_count": q_nonzero_count(psi, base),
        }
    completed_microsteps = macro_completed * MICROSTEPS_PER_MACRO + micro_in_macro
    return {
        "schema": "paper8_full123_progress_v1",
        "status": status,
        "pid": os.getpid(),
        "updated_utc": utc_now(),
        "macro_completed": int(macro_completed),
        "micro_in_current_macro": int(micro_in_macro),
        "completed_microsteps": int(completed_microsteps),
        "products_per_microstep": PRODUCTS_PER_MICRO,
        "logical_product_terms_visited": int(completed_microsteps * PRODUCTS_PER_MICRO),
        "raw_rows_saved": int(writer.total_rows),
        "current_raw_part": int(writer.part_index),
        "P": ps,
        "CD": cds,
        "q": qinfo,
        "target_P": TARGET_P,
        "physical_progress_fraction": frac,
        "elapsed_seconds": elapsed,
        "macrosteps_per_second": rate,
        "eta_seconds_from_state_fraction": eta_seconds,
        "message": last_message,
    }


def write_json_atomic(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def emit_progress(progress_path: Path, payload: dict) -> None:
    write_json_atomic(progress_path, payload)
    p = payload["P"]
    frac = 100.0 * float(payload["physical_progress_fraction"])
    eta = payload["eta_seconds_from_state_fraction"]
    eta_text = "n/a" if eta is None else f"{eta:.1f}s"
    print(
        f"[{payload['updated_utc']}] status={payload['status']} "
        f"macro={payload['macro_completed']} micro={payload['micro_in_current_macro']}/11 "
        f"P(aa,ab,bb)=({p['aa']:.12g},{p['ab']:.12g},{p['bb']:.12g}) "
        f"progress={frac:.6f}% raw_rows={payload['raw_rows_saved']} "
        f"part={payload['current_raw_part']:03d} rate={payload['macrosteps_per_second']:.6g} macro/s "
        f"eta_state={eta_text} product_terms={payload['logical_product_terms_visited']}",
        flush=True,
    )


# -----------------------------------------------------------------------------
# Main experiment. No reference files are opened anywhere in this program.
# -----------------------------------------------------------------------------
def run(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    progress_path = out_dir / "progress.json"
    manifest_path = out_dir / "raw_manifest.json"
    run_info_path = out_dir / "run_info.json"

    # Allocate these exactly once for the uninterrupted run.
    psi = np.empty(NTOT, dtype=np.float64)
    psi_next = np.empty(NTOT, dtype=np.float64)
    S = np.empty((NTOT, NTOT), dtype=np.float64)
    initialize_psi(psi)

    writer = MacroStateWriter(out_dir)
    started_monotonic = time.monotonic()
    started_utc = utc_now()

    run_info = {
        "schema": "paper8_full123_run_info_v1",
        "started_utc": started_utc,
        "condition": "n=1 opposite-sign",
        "channel_CD": {"aa": [0.91, 0.0], "ab": [1.09, 0.36], "bb": [0.91, 0.0]},
        "state_dimension": NTOT,
        "S_dimension": [NTOT, NTOT],
        "microsteps_per_macro": MICROSTEPS_PER_MACRO,
        "products_per_microstep": PRODUCTS_PER_MICRO,
        "products_per_macro": PRODUCTS_PER_MACRO,
        "target_P": TARGET_P,
        "raw_policy": "initial state plus every completed macrostep full psi[123], no thinning",
        "reference_data_used_by_generator": False,
        "state_slicing_used": False,
        "local_41x41_matrices_used": False,
        "zero_state_pruning_used": False,
        "active_q_phase_pruning_used": False,
    }
    write_json_atomic(run_info_path, run_info)

    # Row 0: complete initial macrostate.
    writer.append(0, psi)
    writer.flush()
    emit_progress(
        progress_path,
        progress_payload(
            status="running", macro_completed=0, micro_in_macro=0,
            psi=psi, writer=writer, started_monotonic=started_monotonic,
            last_message="initial full 123-state saved"
        ),
    )

    macro_completed = 0
    last_progress_time = time.monotonic()

    try:
        while not all_channels_reached_target(psi):
            # Exactly 11 microsteps per macrostep.
            for micro in range(MICROSTEPS_PER_MACRO):
                rebuild_full_S(psi, S)            # touches all 15,129 S cells
                full_product_sum(S, psi, psi_next) # visits all 15,129 product terms
                # Synchronous state replacement; psi object itself is retained.
                for i in range(NTOT):
                    psi[i] = psi_next[i]

                now = time.monotonic()
                if now - last_progress_time >= PROGRESS_EVERY_SECONDS:
                    emit_progress(
                        progress_path,
                        progress_payload(
                            status="running", macro_completed=macro_completed,
                            micro_in_macro=micro + 1, psi=psi, writer=writer,
                            started_monotonic=started_monotonic,
                            last_message="inside macrostep; microstate not persisted"
                        ),
                    )
                    last_progress_time = now

            macro_completed += 1
            validate_macro_boundary(psi)

            # Persist EVERY macrostate, with no thinning.
            writer.append(macro_completed, psi)

            # Flush raw data periodically so an external observer sees progress.
            if macro_completed % PROGRESS_EVERY_MACROS == 0:
                writer.flush()
                emit_progress(
                    progress_path,
                    progress_payload(
                        status="running", macro_completed=macro_completed,
                        micro_in_macro=0, psi=psi, writer=writer,
                        started_monotonic=started_monotonic,
                        last_message="macro boundary; full psi[123] persisted"
                    ),
                )
                last_progress_time = time.monotonic()

        writer.flush()
        final_payload = progress_payload(
            status="completed", macro_completed=macro_completed,
            micro_in_macro=0, psi=psi, writer=writer,
            started_monotonic=started_monotonic,
            last_message="all three channels reached P<=20"
        )
        emit_progress(progress_path, final_payload)

    except BaseException as exc:
        writer.flush()
        fail_payload = progress_payload(
            status="failed", macro_completed=macro_completed,
            micro_in_macro=0, psi=psi, writer=writer,
            started_monotonic=started_monotonic,
            last_message=f"{type(exc).__name__}: {exc}"
        )
        emit_progress(progress_path, fail_payload)
        raise

    finally:
        writer.close()
        manifest = {
            "schema": "paper8_full123_raw_manifest_v1",
            "condition": "n=1 opposite-sign",
            "initial_row_included": True,
            "raw_rows_total": int(writer.total_rows),
            "state_dimension": NTOT,
            "parts": writer.parts,
            "completed_utc": utc_now(),
        }
        write_json_atomic(manifest_path, manifest)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Strict full-123 brute-force Paper 8 control run")
    p.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "results_full123_bruteforce_n1_opposite_v1",
        help="output directory for HDF5 raw parts, progress.json, and manifests",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    print("STRICT FULL-123 CONTROL RUN", flush=True)
    print("condition: n=1 opposite-sign; aa=(0.91,0), ab=(1.09,0.36), bb=(0.91,0)", flush=True)
    print("generator reference input: NONE", flush=True)
    print("state: one persistent psi[123]; interaction: one persistent S[123,123]", flush=True)
    print(f"every microstep visits {PRODUCTS_PER_MICRO} product terms; 11 microsteps/macro", flush=True)
    print("raw: initial state + every macrostate full psi[123], no thinning", flush=True)
    run(args.out)


if __name__ == "__main__":
    main()
