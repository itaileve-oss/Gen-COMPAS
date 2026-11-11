#!/usr/bin/env bash
# Memory test: 3 extra 2-ns unbiased runs with NEW random velocities from each heavy-atom TMD endpoint.
# If A-origin endpoints still always go to A and B-origin to B -> origin memory; if mixed -> intermediates are real.
set -uo pipefail
T=/workspace/gc_control/tmd_heavy_test/jobs; MPY=/opt/micromamba/envs/tnfmd/bin/python; export GC_CPU_BLOCKS="0-9;10-19"
for s in 1 2 3; do
  for f in frame_001 frame_002; do
    for P in A B; do
      d=$T/${f}__$P/seed$s; mkdir -p $d
      sed -e "s#^inputState .*#inputState    ../TMD.$P.state.xml#" -e "s/^outputName .*/outputName    unbias.$P/" $T/${f}__$P/Unbiased.$P.conf > $d/Unbiased.$P.conf
      echo "resetVelocities on" >> $d/Unbiased.$P.conf
      ( cd $d && $MPY /workspace/gc/src/gc_md.py --device 0 Unbiased.$P.conf > unbiased.log 2>&1; echo "$f $P seed$s exit $? $(date -u +%T)" ) &
    done
    wait
  done
done
echo "SEEDS DONE $(date -u +%T)"
