#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read-only progress viewer for the strict full-123 Paper 8 control run."""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path


def render(p: dict) -> str:
    P = p.get("P", {})
    frac = 100.0 * float(p.get("physical_progress_fraction", 0.0))
    eta = p.get("eta_seconds_from_state_fraction")
    eta_text = "n/a" if eta is None else f"{float(eta):.1f}s"
    return (
        f"updated={p.get('updated_utc')} status={p.get('status')} "
        f"macro={p.get('macro_completed')} micro={p.get('micro_in_current_macro')}/11 "
        f"P(aa,ab,bb)=({P.get('aa')},{P.get('ab')},{P.get('bb')}) "
        f"progress={frac:.6f}% rows={p.get('raw_rows_saved')} "
        f"part={p.get('current_raw_part')} rate={p.get('macrosteps_per_second')} macro/s "
        f"eta_state={eta_text} product_terms={p.get('logical_product_terms_visited')} "
        f"message={p.get('message','')}"
    )


def read_once(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def main() -> None:
    ap = argparse.ArgumentParser(description="Show progress.json from full123 run")
    ap.add_argument("progress", type=Path, help="path to progress.json")
    ap.add_argument("--watch", action="store_true", help="keep displaying changes")
    ap.add_argument("--interval", type=float, default=5.0, help="watch interval seconds")
    args = ap.parse_args()

    last = None
    while True:
        try:
            p = read_once(args.progress)
            marker = (p.get("updated_utc"), p.get("macro_completed"), p.get("micro_in_current_macro"))
            if marker != last:
                print(render(p), flush=True)
                last = marker
        except FileNotFoundError:
            print(f"progress file not found: {args.progress}", flush=True)
        if not args.watch:
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
