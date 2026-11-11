#!/usr/bin/env bash
# Waits for the real iteration-0 MD (namd stage) to finish, then runs diffusion experiments
# in /workspace/gc_control only. Never touches /workspace/gc/Iterations.
set -uo pipefail
C=/workspace/gc_control; GPY=/opt/micromamba/envs/gen-compas/bin/python; MPY=/opt/micromamba/envs/tnfmd/bin/python
until python3 -c "import json,sys; m=json.load(open('/workspace/gc/Iterations/0th.Iteration/workflow_manifest.json')); sys.exit(0 if m['steps'].get('namd',{}).get('status')=='completed' else 1)" 2>/dev/null; do sleep 120; done
echo "namd stage completed $(date -u +%T)"
# A: same 800 frames, 300 epochs (only change: epochs)
$GPY $C/exp_train.py $C/expA_300ep 300 /workspace/gc/data/initial_noh.dcd 10,2 > $C/expA.log 2>&1
for nz in 10 2; do $MPY $C/gate_any.py $C/expA_300ep/candidates_noise$nz.dcd > $C/expA_gate_noise$nz.txt 2>&1; done
echo "expA done $(date -u +%T)"
# B: initial + all iteration-0 unbiased MD, 60 epochs
$MPY $C/make_iter0_noh.py > $C/expB_data.log 2>&1
$GPY $C/exp_train.py $C/expB_iter0data 60 $C/iter0_plus_initial_noh.dcd 10,2 > $C/expB.log 2>&1
for nz in 10 2; do $MPY $C/gate_any.py $C/expB_iter0data/candidates_noise$nz.dcd > $C/expB_gate_noise$nz.txt 2>&1; done
echo "expB done $(date -u +%T)"
