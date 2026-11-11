"""Heavy-atom TMD test (gc_control only). Pick 2 'between' candidates from expA noise2
(lambda closest to 0.5 on the A->B axis of gate_any.py), hydrogenate + flag heavy atoms."""
import sys, os, numpy as np, mdtraj as md
sys.path.insert(0,'/workspace/gc/src'); sys.path.insert(0,'/workspace/gc/external/Gen-COMPAS')
from gc_cv_check import frame_cvs, ca_map
from tools.occupancy import hydrogenate_and_set_occupancy
T='/workspace/gc_control/tmd_heavy_test'; os.makedirs(f'{T}/noh',exist_ok=True)
top='/workspace/gc/gencompas/TOP_files/protein_noh.pdb'
g=md.load('/workspace/gc_control/expA_300ep/candidates_noise2.dcd',top=top)
t=md.load('/workspace/gc/data/initial_noh.dcd',top=top); idx=ca_map(g.topology)
cv=lambda tr: np.array([[(r:=frame_cvs(x*10,idx))['S1_CA'],r['S2_CA'],r['S3_CA'],r['theta_spread']] for x in tr.xyz])
G=cv(g); A=cv(t[:400]); B=cv(t[400:]); mA,mB=A.mean(0),B.mean(0); sc=np.sqrt((A.var(0)+B.var(0))/2); w=(mB-mA)/sc
L=((G-mA)/sc)@w/(w@w)
ca=g.topology.select('name CA'); ch=np.array([g.topology.atom(i).residue.chain.index for i in ca]); same=ch[1:]==ch[:-1]
bad=(np.linalg.norm(np.diff(g.xyz[:,ca]*10,axis=1),axis=2)[:,same]>4.2).mean(1)
cand=np.where((L>0.3)&(L<0.7))[0]; cand=cand[np.argsort(bad[cand])]      # cleanest geometry first
pick=[cand[0]]
for c in cand[1:]:
    if abs(L[c]-L[pick[0]])>0.1 or len(cand)<3: pick.append(c); break
for k,i in enumerate(pick,1):
    g[i].save_pdb(f'{T}/noh/frame_{k:03d}.pdb')
    print(f'frame_{k:03d} = cand {i}: lambda {L[i]:.2f} S1 {G[i,0]:.2f} S2 {G[i,1]:.2f} S3 {G[i,2]:.2f} theta {G[i,3]:.2f} CA-CA>4.2 {bad[i]:.3f}')
hydrogenate_and_set_occupancy(f'{T}/noh','/workspace/gc/gencompas/TOP_files/protein.psf',
    '/workspace/gc/gencompas/TOP_files/protein.pdb',f'{T}/targets','element != H')
