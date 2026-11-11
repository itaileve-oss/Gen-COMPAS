#!/usr/bin/env bash
# v2 endpoints, FROM SCRATCH (2026-11-07). Replaces build_endpoints.sh, whose B endpoint
# collapsed toward symmetric during an unrestrained 0.5 ns hot-start equilibration.
#   * staged restrained equilibration (~9 ns, src/equil_protocol.py) before production
#   * bioRxiv 10.1101/2025.11.28.690881 settings: ff14SB/TIP3P/0.15 M KCl, 310 K,
#     HMR 4 fs, vdW cutoff 12 A switched from 10 A (prep_system.py defaults; the old
#     script overrode them with 10/9)
# Stops after the 5 ns unbiased runs + initial_noh.dcd: the B endpoint must be checked
# (does it stay asymmetric?) BEFORE any Gen-COMPAS training.
# Run ON the instance from /workspace/gc. Log: /workspace/gc/build_endpoints_v2.log
set -euo pipefail
GC=/workspace/gc; PY=/opt/micromamba/envs/tnfmd/bin/python
cd "$GC"
export GC_CPU_BLOCKS="0-9;10-19"

# 1) archive everything derived from the v1 endpoints (nothing is deleted)
ARCH=archive_v1_collapsedB
mkdir -p $ARCH
for x in systems Iterations data/initial_noh.dcd analysis iter0_diffusion.log build_endpoints.log; do
  [ -e "$x" ] && mv "$x" "$ARCH/" || true
done
mkdir -p systems

# 2) A fix stage first (B copies its HIS tautomers), then full A and B in parallel
$PY src/prep_system.py data/1TNF.pdb systems/A --mutate LEU-143-ASP --renumber-start 6 --fix-only
cp systems/A/fixed.pdb systems/A_fixed_ref.pdb
COMMON="--renumber-start 6 --platform CUDA"
($PY -u src/prep_system.py data/1TNF.pdb systems/A --mutate LEU-143-ASP $COMMON > systems/A.prep.log 2>&1) &
($PY -u src/prep_system.py data/7JRA.pdb systems/B --chain-map C:A,B:B,A:C --res-range 82-233 \
    --graft-loops-from C --variants-from systems/A_fixed_ref.pdb --pocket-ref data/7JRA.pdb $COMMON > systems/B.prep.log 2>&1) &
wait
echo "=== equilibration done ==="

# 3) 5 ns unbiased production from each (unchanged settings: 10 ps frames, C->A CVs)
for P in A B; do
  cat > systems/$P/unbiased.conf <<CONF
systemXML    $GC/systems/$P/system.xml
proteinPDB   $GC/systems/$P/equilibrated_protein.pdb
inputState   $GC/systems/$P/equilibrated.xml
outputName   unbiased
temperature  310
timestep     4.0
run          1250000
dcdfreq      2500
cvPair       C,A
CONF
done
(cd systems/A && $PY $GC/src/gc_md.py --device 0 unbiased.conf > unbiased.out 2>&1) &
(cd systems/B && $PY $GC/src/gc_md.py --device 0 unbiased.conf > unbiased.out 2>&1) &
wait
for P in A B; do cp systems/$P/unbiased.state.xml systems/$P/start.xml; done
mkdir -p data
/opt/micromamba/envs/gen-compas/bin/python gencompas/combine_initial.py
echo "=== endpoints v2 done - CHECK B BEFORE TRAINING ==="
