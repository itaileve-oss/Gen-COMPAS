import mdtraj as md, numpy as np
pos={'ARG':['NE','NH1','NH2'],'LYS':['NZ'],'HIS':['ND1','NE2']}; neg={'GLU':['OE1','OE2'],'ASP':['OD1','OD2']}
def sb(t):
    top=t.topology; P=[(a.index,a.residue) for a in top.atoms if a.residue.name in pos and a.name in pos[a.residue.name]]
    N=[(a.index,a.residue) for a in top.atoms if a.residue.name in neg and a.name in neg[a.residue.name]]
    pairs=set()
    for i,ri in P:
        for j,rj in N:
            if np.linalg.norm(t.xyz[0,i]-t.xyz[0,j])*10<4.0: pairs.add((f'{"ABC"[ri.chain.index]}:{ri.name}{ri.resSeq}',f'{"ABC"[rj.chain.index]}:{rj.name}{rj.resSeq}'))
    return pairs
A=sb(md.load('/workspace/gc/systems/A/equilibrated_protein.pdb')); B=sb(md.load('/workspace/gc/systems/B/equilibrated_protein.pdb'))
keys=('ARG103','GLU104','ARG6','GLU146','HIS73')
f=lambda S: sorted(p for p in S if any(k in p[0] or k in p[1] for k in keys))
print('only in A:', f(A-B)); print('only in B:', f(B-A)); print('both:', f(A&B))
T='/workspace/gc_control/tmd_heavy_test/xpull'
for d,P in (('frame_001__A','A'),('frame_001__B','B')):
    e=md.load(f'{T}/{d}/TMD.{P}.dcd',top=f'{T}/{d}/output.pdb')[-1]; E=sb(e)
    print(d,'pulled end keeps A-only bridges:',len(f(E&(A-B))),'/',len(f(A-B)),' B-only:',len(f(E&(B-A))),'/',len(f(B-A)))
