"""Memorisation check: CA RMSD of each generated frame to its NEAREST training frame,
vs the same statistic for held-in training frames to each other (excluding self + +/-5 neighbours)."""
import sys, mdtraj as md, numpy as np
top='/workspace/gc/gencompas/TOP_files/protein_noh.pdb'
g=md.load(sys.argv[1],top=top); t=md.load(sys.argv[2],top=top)
ca=g.topology.select('name CA'); g=g.atom_slice(ca)[::5]; t=t.atom_slice(ca)
ts=t[::max(1,len(t)//1600)]
near=np.array([md.rmsd(ts,g,i).min() for i in range(len(g))])*10
ref=t[::max(1,len(t)//200)]; nn=[]
for i in range(len(ref)):
    r=md.rmsd(ts,ref,i)*10; r=r[r>1e-3]; nn.append(np.sort(r)[0])
nn=np.array(nn)
print(f"gen->nearest-train CA RMSD median {np.median(near):.2f} p5 {np.percentile(near,5):.2f} p95 {np.percentile(near,95):.2f} frac<0.5 {(near<0.5).mean():.3f}")
print(f"train->nearest-other-train median {np.median(nn):.2f}")
print(f"COPY_MEDIAN {np.median(near):.3f}")
