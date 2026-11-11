"""Cross-pull targets: A-origin is pulled to where the B-origin heavy-atom TMD ended, and vice versa.
Coordinates = last frame of the other origin's TMD dcd; atom order/occupancy (heavy=1) from targets/frame_NNN.pdb."""
import sys, numpy as np, mdtraj as md
sys.path.insert(0,'/workspace/gc/src'); from gc_cv_check import frame_cvs, ca_map
T='/workspace/gc_control/tmd_heavy_test'
for f in ('frame_001','frame_002'):
    tpl=f'{T}/targets/{f}.pdb'; lines=open(tpl).read().splitlines()
    for P,other in (('A','B'),('B','A')):
        tr=md.load(f'{T}/jobs/{f}__{other}/TMD.{other}.dcd',top=tpl)[-1]
        xyz=tr.xyz[0]*10; k=0; out=[]
        for l in lines:
            if l.startswith(('ATOM','HETATM')):
                x,y,z=xyz[k]; k+=1; l=f'{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}'
            out.append(l)
        assert k==tr.n_atoms
        open(f'{T}/xpull_{f}_for{P}.pdb','w').write('\n'.join(out)+'\n')
        noh=tr.atom_slice(tr.topology.select('element != H')); r=frame_cvs(noh.xyz[0]*10,ca_map(noh.topology))
        print(f'{f} target for {P}-origin = {other}-origin TMD end: S1 {r["S1_CA"]:.2f} S3 {r["S3_CA"]:.2f} theta {r["theta_spread"]:.2f}')
