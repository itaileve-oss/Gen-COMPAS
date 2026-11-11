#!/usr/bin/env bash
# Heavy-atom TMD test: 2 targets x {A,B} origin, TMD 500 ps on heavy atoms, k scaled x7.83
# (NAMD k/N convention: same per-atom stiffness as the CA-only k=1e5), then 2 ns unbiased.
set -uo pipefail
T=/workspace/gc_control/tmd_heavy_test; MPY=/opt/micromamba/envs/tnfmd/bin/python; GPY=/opt/micromamba/envs/gen-compas/bin/python
GC=/workspace/gc; TPL=$GC/Iterations/0th.Iteration/namd/frame_001__A
$GPY $T/../tmdh_setup.py > $T/setup.log 2>&1 || { echo setup failed; exit 1; }
cat $T/setup.log | grep -E "^frame|Selection"
export GC_CPU_BLOCKS="0-9;10-19"
for f in $T/targets/frame_*.pdb; do
  n=$(basename $f .pdb)
  for P in A B; do
    j=$T/jobs/${n}__$P; mkdir -p $j; cp $f $j/output.pdb
    sed -e "s/^run .*/run           125000         # 500 ps heavy-atom TMD test/" -e "s/^TMDLastStep .*/TMDLastStep   125000/" \
        -e "s/^TMDk .*/TMDk          783000/" -e "s#systems/A/#systems/$P/#g" -e "s/TMD\.A/TMD.$P/" $TPL/TMD.A.conf > $j/TMD.$P.conf
    sed -e "s#systems/A/#systems/$P/#g" -e "s/TMD\.A/TMD.$P/" -e "s/unbias\.A/unbias.$P/" $TPL/Unbiased.A.conf > $j/Unbiased.$P.conf
    ( cd $j && $MPY $GC/src/gc_md.py --device 0 TMD.$P.conf > tmd.log 2>&1 && $MPY $GC/src/gc_md.py --device 0 Unbiased.$P.conf > unbiased.log 2>&1; echo "$n $P exit $? $(date -u +%T)" ) &
  done
  wait
done
echo "TEST DONE $(date -u +%T)"
