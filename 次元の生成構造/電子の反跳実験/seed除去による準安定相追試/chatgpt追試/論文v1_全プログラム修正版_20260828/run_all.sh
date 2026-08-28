#!/bin/bash
set -e; cd "$(dirname "$0")"; export LC_ALL=en_US.UTF-8 MPLBACKEND=Agg
O=../論文v1_全再現テスト_20260828/original; [ -d "$O" ] || { echo "原本 $O が必要（論文v1_全再現テスト_20260828/run_all.sh で展開）"; exit 1; }
rm -rf fixed fixed_baseline; mkdir fixed; cp -R "$O"/* fixed/; find fixed -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
( cd fixed && patch -p1 -N -s < ../../論文v1_全再現テスト_20260828/results/path_patches.diff )
python3 apply_fixes.py; python3 verify_fixes_all.py | tee results/verify_fixes_all.log
cp -R fixed fixed_baseline
./run_all_fixed.sh fixed amplitude; ./run_all_fixed.sh fixed_baseline phase
python3 compare_three_way.py
NF="$(cd ../../../.. && pwd)/note_figs_self_consistent_inflation"
for R in fixed fixed_baseline; do OUTD=note_figs_fixed/$R; mkdir -p "$OUTD/zip"; ( cd $R/complex_simplex_decompactification_N5_N16_20260826 && zip -q -r "../../$OUTD/zip/complex_simplex_decompactification_N5_N16_20260826.zip" results/N5_geometry_summary.csv results/N16_geometry_summary.csv ); python3 note_figs_fixed/make_note_figs_fixed.py "$OUTD" "$R/N5_dynamics_followup_theorems_and_stability_20260826/followup_dynamics_20260826"; python3 note_figs_fixed/make_note_fig0_fixed.py "$OUTD" "$OUTD/zip"; done
