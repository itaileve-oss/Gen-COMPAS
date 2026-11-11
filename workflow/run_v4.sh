#!/usr/bin/env bash
# After the cross-pull test: summarise it, then the fresh v4 Gen-COMPAS run, iterations 0-3.
# On a failed step, retry with --resume (max 3 retries). Log: /workspace/gc/v4_run.log
set -uo pipefail
C=/workspace/gc_control; T=$C/tmd_heavy_test
until grep -q "XPULL DONE" $T/xpull.log; do sleep 120; done
{ echo "cross-pull summary $(date -u)"; cd $T/xpull
  for d in */seed*; do f=$(ls $d/unbias.*.colvars.traj); awk -v j=$d '!/^#/{n++; a[n]=$2; b[n]=$4} END{s=0;t=0;for(i=n-19;i<=n;i++){s+=a[i];t+=b[i]}; printf "%s last-400ps S1 %.2f S3 %.2f\n", j, s/20, t/20}' $f; done
  for d in */; do tail -n 1 $d/TMD.*.colvars.traj | awk -v j=$d '{printf "%s TMD end S1 %.2f S3 %.2f rmsd %.2f\n", j, $2, $4, $6}'; done; } > $T/xpull_summary.txt 2>&1
cd /workspace/gc/gencompas; export GC_CPU_BLOCKS="0-9;10-19"
GC=/opt/micromamba/envs/gen-compas/bin/gen-compas
echo "=== v4 start $(date -u)"
$GC --config tnf.v4.workflow.yaml --iteration 0 1 2 3; rc=$?
for try in 1 2 3; do
  [ $rc -eq 0 ] && break
  echo "=== v4 exited $rc, resume attempt $try $(date -u)"
  sleep 60
  $GC --config tnf.v4.workflow.yaml --iteration 0 1 2 3 --resume; rc=$?
done
echo "=== v4 finished rc=$rc $(date -u)"
