import mdtraj as md, numpy as np
D='/workspace/gc_control/trpcage'
top=D+'/Dataset/trp_noh.psf'
t=md.load(D+'/Dataset/0th.Iteration.noh.dcd',top=top)
g=md.load(D+'/workflow_output/0th.Iteration/generated/candidates.dcd',top=top)
fold=md.load(D+'/NAMD_inputs/folded.pdb'); unf=md.load(D+'/NAMD_inputs/unfolded.pdb')
ca=t.topology.select('name CA')
print('training frames',t.n_frames,'atoms',t.n_atoms,'| generated',g.n_frames)
def ref_rmsd(tr,ref):
    rc=ref.topology.select('name CA'); return md.rmsd(tr,ref,atom_indices=ca,ref_atom_indices=rc)*10
for name,tr in [('train',t),('gen',g)]:
    x=tr.xyz[:,ca]*10; d=np.linalg.norm(x[:,1:]-x[:,:-1],axis=2)
    s=tr[::max(1,tr.n_frames//200)]; r=np.array([md.rmsd(s,s,k,atom_indices=ca) for k in range(s.n_frames)])*10
    rf=ref_rmsd(tr,fold); ru=ref_rmsd(tr,unf)
    print(f"{name:5s} CA-CA {d.mean():.2f}+/-{d.std():.2f} frac>4.2 {(d>4.2).mean():.3f} | pairwise CA RMSD median {np.median(r[np.triu_indices(s.n_frames,1)]):.2f} A"
          f" | RMSD to folded pct5/50/95 {np.percentile(rf,[5,50,95]).round(2)} | to unfolded {np.percentile(ru,[5,50,95]).round(2)}")
    print('      hist RMSD-to-folded bins 0..12 step 1:',np.histogram(rf,bins=np.arange(0,13,1))[0])
