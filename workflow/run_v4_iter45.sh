#!/usr/bin/env bash
# Resume v4 after train_committor failed (torch.compile needs g++, absent). Disable dynamo -> eager mode.
set -uo pipefail
export TORCHDYNAMO_DISABLE=1 GC_CPU_BLOCKS="0-9;10-19"
cd /workspace/gc/gencompas
GC=/opt/micromamba/envs/gen-compas/bin/gen-compas
rc=1
for try in 1 2 3; do
  echo "=== v4 iter4-5 (no-compile) attempt $try $(date -u)"
  $GC --config tnf.v4.workflow.yaml --iteration 4 5 --resume; rc=$?
  [ $rc -eq 0 ] && break
  sleep 60
done
echo "=== v4 iter4-5 finished rc=$rc $(date -u)"
