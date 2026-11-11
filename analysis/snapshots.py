"""Remote: CA RMSD of every iteration-3 frame to 1TNF and 7JRA; write real (not averaged)
snapshots: highest-weight frame in the +/-0.3 A window of the A and B FES minima, the frame
closest to 7JRA, and the frame closest to 1TNF."""
import numpy as np, pandas as pd, mdtraj as md, sys
import os
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","src")); from gc_cv_check import ca_map, frame_cvs
RW=os.environ.get("RW_DIR","/workspace/gc/Iterations_v4/3rd.Iteration/riteweight")
top=md.load_psf(f"{RW}/diffusion_training.psf")
for c,n in zip(top.chains,"ABC"): c.chain_id=n
T=md.load(f"{RW}/diffusion_training.dcd",top=top); d=pd.read_csv(f"{RW}/frame_weights.csv")
def cakeys(t): return {(a.residue.chain.chain_id,a.residue.resSeq):a.index for a in t.topology.atoms if a.name=="CA"}
ks=cakeys(T); res={}
for name,f in [("1TNF","1TNF_crystal_symmetric.pdb"),("7JRA","7JRA_crystal_asymmetric.pdb")]:
    X=md.load(f); kx=cakeys(X); keys=sorted(set(ks)&set(kx))
    it=[ks[k] for k in keys]; ix=[kx[k] for k in keys]
    ref=md.Trajectory(X.xyz[:,ix],None); sub=md.Trajectory(T.xyz[:,it],None)
    res[name]=10*md.rmsd(sub,ref,0); print(name,"common CA",len(keys))
    r=frame_cvs(X.xyz[0]*10,ca_map(X.topology)); print(f"  crystal CVs S1 {r['S1_CA']:.2f} S2 {r['S2_CA']:.2f} S3 {r['S3_CA']:.2f}")
d["rmsd_1TNF"]=res["1TNF"]; d["rmsd_7JRA"]=res["7JRA"]; d[["S1","S3","fel_weight","rmsd_1TNF","rmsd_7JRA"]].to_csv("frame_rmsd_crystals.csv",index=False)
print("RMSD to 1TNF: min %.2f median %.2f | RMSD to 7JRA: min %.2f median %.2f"%(d.rmsd_1TNF.min(),d.rmsd_1TNF.median(),d.rmsd_7JRA.min(),d.rmsd_7JRA.median()))
print("frames closer to 7JRA than to 1TNF: %d / %d (weight %.3f)"%((d.rmsd_7JRA<d.rmsd_1TNF).sum(),len(d),d.fel_weight[d.rmsd_7JRA<d.rmsd_1TNF].sum()))
picks={}
for nm,(c1,c3) in {"A_min":(20.0,9.8),"B_min":(22.7,12.0)}.items():
    m=(abs(d.S1-c1)<=0.3)&(abs(d.S3-c3)<=0.3); picks[nm]=int(d.fel_weight[m].idxmax())
picks["closest_to_7JRA"]=int(d.rmsd_7JRA.idxmin()); picks["closest_to_1TNF"]=int(d.rmsd_1TNF.idxmin())
for nm,i in picks.items():
    T[i].save_pdb(f"sim_{nm}_frame{i}.pdb")
    print(f"{nm:16s} frame {i:5d}: S1 {d.S1[i]:.2f} S3 {d.S3[i]:.2f} | RMSD to 1TNF {d.rmsd_1TNF[i]:.2f}  to 7JRA {d.rmsd_7JRA[i]:.2f} | weight {d.fel_weight[i]:.2e}")
for lo,hi in [(0,22),(22,23.5),(23.5,30)]:
    m=(d.S1>=lo)&(d.S1<hi); print(f"S1 {lo}-{hi}: n={m.sum()} median RMSD 1TNF {d.rmsd_1TNF[m].median():.2f}, 7JRA {d.rmsd_7JRA[m].median():.2f}")
