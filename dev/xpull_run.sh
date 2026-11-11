#!/usr/bin/env bash
# Cross-pull test (user idea 2026-11-08): pull each origin to the OTHER origin's heavy-atom TMD endpoint,
# then 3 x 2 ns unbiased with fresh velocities. Waits for tmdh_seeds.sh to finish.
set -uo pipefail
C=/workspace/gc_control; T=$C/tmd_heavy_test; MPY=/opt/micromamba/envs/tnfmd/bin/python; export GC_CPU_BLOCKS="0-9;10-19"
until grep -q "SEEDS DONE" $T/seeds.log; do sleep 60; done
$MPY $C/xpull_targets.py 2>&1 | grep -v dcdplugin || exit 1
for f in frame_001 frame_002; do
  for P in A B; do
    j=$T/xpull/${f}__$P; mkdir -p $j; cp $T/xpull_${f}_for$P.pdb $j/output.pdb
    cp $T/jobs/${f}__$P/TMD.$P.conf $j/
    ( cd $j && $MPY /workspace/gc/src/gc_md.py --device 0 TMD.$P.conf > tmd.log 2>&1; echo "$f $P xTMD exit $? $(date -u +%T)" ) &
  done
  wait
done
for s in 1 2 3; do
  for f in frame_001 frame_002; do
    for P in A B; do
      d=$T/xpull/${f}__$P/seed$s; mkdir -p $d
      sed -e "s#^inputState .*#inputState    ../TMD.$P.state.xml#" -e "s/^outputName .*/outputName    unbias.$P/" $T/jobs/${f}__$P/Unbiased.$P.conf > $d/Unbiased.$P.conf
      echo "resetVelocities on" >> $d/Unbiased.$P.conf
      ( cd $d && $MPY /workspace/gc/src/gc_md.py --device 0 Unbiased.$P.conf > unbiased.log 2>&1; echo "$f $P seed$s exit $? $(date -u +%T)" ) &
    done
    wait
  done
done
echo "XPULL DONE $(date -u +%T)"
