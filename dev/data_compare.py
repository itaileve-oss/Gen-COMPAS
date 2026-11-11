import sys; sys.path.insert(0,'/workspace/gc/external/Gen-COMPAS')
from utils.mdtraj_io import load_dcd
import mdtraj as md, numpy as np
sets={'TNF':('/workspace/gc/gencompas/TOP_files/protein_noh.psf','/workspace/gc/data/initial_noh.dcd'),
      'TRP':('/workspace/gc_control/trpcage/Dataset/trp_noh.psf','/workspace/gc_control/trpcage/Dataset/0th.Iteration.noh.dcd')}
for name,(psf,dcd) in sets.items():
    t=load_dcd(dcd,top=psf); top=t.topology
    print(f"\n=== {name}: {t.n_frames} frames, {t.n_atoms} atoms, {top.n_residues} residues, {top.n_chains} chains, {top.n_bonds} bonds")
    # bond sanity in frame 0 and across frames
    b=np.array([[x.index,y.index] for x,y in top.bonds])
    L=np.linalg.norm(t.xyz[:,b[:,0]]-t.xyz[:,b[:,1]],axis=2)*10
    print(f"  bond lengths (A) over all frames: median {np.median(L):.2f}, max {L.max():.2f}, frac>2.0 {(L>2.0).mean():.4f}")
    nbonded=np.setdiff1d(np.arange(t.n_atoms),b.ravel()); print('  atoms with no bond:',len(nbonded))
    # inter-chain bonds?
    cb=sum(1 for x,y in top.bonds if x.residue.chain.index!=y.residue.chain.index); print('  bonds crossing chains:',cb)
    # peptide bonds C-N between consecutive residues
    pep=sum(1 for x,y in top.bonds if {x.name,y.name}=={'C','N'} and x.residue.index!=y.residue.index); print('  peptide C-N bonds:',pep)
    ca=top.select('name CA')
    x=t.xyz[:,ca]*10; same=np.array([top.atom(i).residue.chain.index for i in ca]); same=same[1:]==same[:-1]
    d=np.linalg.norm(x[:,1:]-x[:,:-1],axis=2)[:,same]; print(f"  CA-CA: {d.mean():.2f}+/-{d.std():.2f}, max {d.max():.2f}")
    if top.n_chains>1:
        com=np.stack([t.xyz[:,[a.index for a in c.atoms]].mean(1) for c in top.chains],1)*10
        dd=np.linalg.norm(com[:,[0,1,2]]-com[:,[1,2,0]],axis=2); print(f"  inter-chain COM distances: min {dd.min():.1f} max {dd.max():.1f}")
    # coordinates: centred? spread
    cen=t.xyz.mean(1)*10; print(f"  raw frame-centre spread (A): std {cen.std(0).round(2)}  range {np.ptp(cen,0).round(1)}")
    raw=md.rmsd(t,t,0,precentered=False)*10
    t2=t[:]; t2.superpose(t2[0]); al=np.sqrt((((t2.xyz-t2.xyz[0])*10)**2).sum(2).mean(1))
    print(f"  RMSD to frame 0: optimal-fit median {np.median(raw):.2f} | after loader-style superpose (no-refit) median {np.median(al):.2f}")
    # per-atom std after alignment (how much each atom moves) and overall coordinate std used for normalisation
    print(f"  per-axis coord std after superpose (nm): {t2.xyz.std(axis=(0,1)).round(3)} ; per-atom fluctuation median {np.median(t2.xyz.std(0))*10:.2f} A")
    # psf residue/segment info
    print('  first chain ids/segments:', [c.chain_id if hasattr(c,'chain_id') else c.index for c in top.chains][:4], ' first residues:', [str(r) for r in list(top.residues)[:2]], ' last:', str(list(top.residues)[-1]))
