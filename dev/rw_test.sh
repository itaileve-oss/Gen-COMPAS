#!/usr/bin/env bash
# RiteWeight robustness test on iteration-0 data, control folder only.
C=/workspace/gc_control/rw_test; mkdir -p $C
GPY=/opt/micromamba/envs/gen-compas/bin/python
for K in 50 30; do
  for SEED in 2026 7; do
    d=$C/k${K}_s${SEED}; mkdir -p $d
    $GPY - <<PY
import yaml
c=yaml.safe_load(open('/workspace/gc/Iterations/0th.Iteration/effective_config.yaml'))
rw=c['RiteWeight']; rw['io']['out']='$d'; rw['riteweight']['n_clusters']=$K; rw['riteweight']['seed']=$SEED
def redirect(o):
    if isinstance(o,dict):
        return {k:(v.replace('/workspace/gc/Iterations/0th.Iteration/riteweight','$d') if isinstance(v,str) else redirect(v)) for k,v in o.items()}
    return o
rw['outputs']=redirect(rw.get('outputs',{}))
yaml.safe_dump(c,open('$d/config.yaml','w'))
print('outputs ->',rw['outputs'])
PY
    (cd /workspace/gc/external/Gen-COMPAS && OMP_NUM_THREADS=4 taskset -c 16-19 $GPY -m common.runner --step riteweight --config $d/config.yaml > $d/rw.log 2>&1); echo "K=$K seed=$SEED exit $?" | tee -a $C/summary.txt
    grep -a -i "converged\|delta\|ESS\|effective\|error" $d/rw.log | tail -3 >> $C/summary.txt
  done
done
