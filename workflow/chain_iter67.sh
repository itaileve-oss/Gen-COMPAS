#!/usr/bin/env bash
# wait for iter4-5 launcher to exit; start iter6-7 only if it ended rc=0
PID=${1:?usage: chain_iter67.sh PID_OF_RUNNING_LAUNCHER}
while kill -0 $PID 2>/dev/null; do sleep 300; done
if grep -q "iter4-5 finished rc=0" /workspace/gc/v4_iter45_run.log; then
  /workspace/gc_control/run_v4_iter67.sh > /workspace/gc/v4_iter67_run.log 2>&1
else
  echo "iter4-5 did not finish rc=0; not starting 6-7 $(date -u)" > /workspace/gc/v4_iter67_run.log
fi
