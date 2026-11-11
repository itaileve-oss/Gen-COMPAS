import sys; sys.path.insert(0,'/workspace/gc/external/Gen-COMPAS')
from utils.mdtraj_io import load_dcd
import numpy as np
t=load_dcd('/workspace/gc/data/initial_noh.dcd',top='/workspace/gc/gencompas/TOP_files/protein_noh.psf')
t.superpose(t[0])
X=t.xyz*10; A=X[:400]; B=X[400:]
mA,mB=A.mean(0),B.mean(0)
d=np.linalg.norm(mB-mA,axis=1)                     # mean A -> mean B displacement per atom
fl=np.sqrt(((X-X.mean(0))**2).sum(2).mean(0))      # per-atom RMS fluctuation over all 800 frames
print('A-mean vs B-mean per-atom displacement (A): median %.2f  p75 %.2f  p90 %.2f  p99 %.2f  max %.2f'%tuple(np.percentile(d,[50,75,90,99,100])))
print('atoms moving >3 A: %d  >5 A: %d  >10 A: %d  (of %d)'%((d>3).sum(),(d>5).sum(),(d>10).sum(),len(d)))
print('per-atom RMS fluctuation over all frames: median %.2f p90 %.2f max %.2f'%tuple(np.percentile(fl,[50,90,100])))
top=t.topology
big=np.argsort(d)[::-1][:12]
print('largest movers:',[(str(top.atom(i).residue),top.atom(i).residue.chain.index,top.atom(i).name,round(float(d[i]),1)) for i in big])
# by chain
for c in top.chains:
    idx=[a.index for a in c.atoms]; print('chain',c.index,'median disp %.2f p90 %.2f'%(np.median(d[idx]),np.percentile(d[idx],90)))
# normalisation scale used by the loader
print('loader global std per axis (A):',(X.std(axis=(0,1))).round(2))
