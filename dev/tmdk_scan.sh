#!/usr/bin/env bash
# Scan the TMD force constant on a few REAL iteration-0 targets before the namd stage.
# v1 test: TMDk 10000 pulled the C->A CVs only ~13% of the way to a far target in 100 ps.
# Usage (on the GC instance): bash gencompas/tmdk_scan.sh TARGET1.pdb [TARGET2.pdb ...]
set -euo pipefail
GC=/workspace/gc; PY=/opt/micromamba/envs/tnfmd/bin/python
export GC_CPU_BLOCKS="0-9;10-19"
S=$GC/tmdk_scan; mkdir -p $S
for tgt in "$@"; do
  tn=$(basename "$tgt" .pdb)
  for K in 10000 30000 100000; do
    for P in A B; do
      j=$S/${tn}_k${K}_$P; mkdir -p $j; cp "$tgt" $j/target.pdb
      sed -e "s|^TMDk .*|TMDk          $K|" -e "s|^TMDFile .*|TMDFile       target.pdb|" \
          -e "s|^outputName .*|outputName    TMD.$P|" $GC/gencompas/NAMD_inputs/TMD.$P.conf > $j/TMD.conf
      (cd $j && $PY $GC/src/gc_md.py --device 0 TMD.conf > tmd.out 2>&1) &
    done
    wait
  done
done
echo "target k protocol  rmsd_start rmsd_end  S1 S3 theta (end)"
for j in $S/*_k*_[AB]; do
  f=$j/TMD.*.colvars.traj
  awk -v n=$(basename $j) '!/^#/{if(!s){s=$6}; e=$6; a=$2; c=$4; t=$5} END{printf "%-40s %6.2f %6.2f  %6.2f %6.2f %5.2f\n", n, s, e, a, c, t}' $f
done
