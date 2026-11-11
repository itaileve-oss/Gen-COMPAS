"""Which residues does the cross-pull fail to move? Per-residue heavy-atom deviation (after global
heavy-atom superposition) between each cross-pull TMD endpoint and its target."""
import mdtraj as md, numpy as np
T='/workspace/gc_control/tmd_heavy_test'
for f in ('frame_001','frame_002'):
    for P in ('A','B'):
        tgt=md.load(f'{T}/xpull/{f}__{P}/output.pdb'); end=md.load(f'{T}/xpull/{f}__{P}/TMD.{P}.dcd',top=f'{T}/xpull/{f}__{P}/output.pdb')[-1]
        h=tgt.topology.select('element != H'); end.superpose(tgt,atom_indices=h)
        d=np.linalg.norm(end.xyz[0,h]-tgt.xyz[0,h],axis=1)*10
        res={}
        for i,a in zip(h,d):
            r=tgt.topology.atom(i).residue; res.setdefault((r.chain.index,r.resSeq,r.name),[]).append(a)
        rr=sorted(((np.sqrt(np.mean(np.square(v))),k) for k,v in res.items()),reverse=True)
        tot=np.sqrt(np.mean(d**2)); top=rr[:12]
        share=sum(len(res[k])*v**2 for v,k in rr[:20])/np.sum(d**2)
        print(f'{f} {P}-origin: heavy RMSD {tot:.2f} A; top-20 residues carry {100*share:.0f}% of the squared deviation')
        print('   '+', '.join(f'{"ABC"[k[0]]}:{k[2]}{k[1]} {v:.1f}' for v,k in top))
