import sys; sys.path.insert(0,'/workspace/gc/src')
import mdtraj as md, numpy as np
from gc_cv_check import frame_cvs, ca_map
top='/workspace/gc/gencompas/TOP_files/protein_noh.pdb'
g=md.load(''+sys.argv[1]+'',top=top)
t=md.load('/workspace/gc/data/initial_noh.dcd',top=top)
idx=ca_map(g.topology)
def cvs(tr): return np.array([[ (r:=frame_cvs(x*10,idx))['S1_CA'], r['S2_CA'], r['S3_CA'], r['theta_spread']] for x in tr.xyz])
G=cvs(g); A=cvs(t[:400]); B=cvs(t[400:])
pass
names=['S1','S2','S3','theta']
for n,X in [('A',A),('B',B),('GEN',G)]:
    print(n.ljust(4),'  '.join(f"{c} {X[:,i].mean():.2f}+/-{X[:,i].std():.2f} [{np.percentile(X[:,i],5):.2f},{np.percentile(X[:,i],95):.2f}]" for i,c in enumerate(names)))
# progress along A->B axis in standardized S1,S3
mA,mB=A.mean(0),B.mean(0); d=mB-mA; sc=np.sqrt((A.var(0)+B.var(0))/2)
def lam(X): v=(X-mA)/sc; w=d/sc; return v@w/(w@w)
L=lam(G)
print('gen lambda (0=A,1=B) pct 5/25/50/75/95:',np.percentile(L,[5,25,50,75,95]).round(2))
print('histogram lambda bins [-0.5..1.5 step .25]:',np.histogram(L,bins=np.arange(-0.5,1.76,0.25))[0])
print('between (0.2<lam<0.8):',((L>0.2)&(L<0.8)).sum(),' A-like (<0.2):',(L<=0.2).sum(),' B-like (>=0.8):',(L>=0.8).sum())
ca=g.topology.select('name CA'); ch=np.array([g.topology.atom(i).residue.chain.index for i in ca]); same=ch[1:]==ch[:-1]
for n,tr in [('gen',g),('train',t)]:
    x=tr.xyz[:,ca]*10; dd=np.linalg.norm(x[:,1:]-x[:,:-1],axis=2)[:,same]
    s=tr.atom_slice(ca)[::max(1,len(tr)//100)]; r=np.array([md.rmsd(s,s,k) for k in range(len(s))])*10
    print(f"{n}: CA-CA {dd.mean():.2f}+/-{dd.std():.2f} frac>4.2 {(dd>4.2).mean():.3f} | pairwise CA RMSD median {np.median(r[np.triu_indices(len(s),1)]):.2f}")
