import numpy as np, glob, mdtraj as md, pandas as pd
D='/workspace/gc/Iterations/0th.Iteration/namd'
top='/workspace/gc/gencompas/TOP_files/protein.psf'
rows=[]
for f in sorted(glob.glob(D+'/*/unbias.[AB].colvars.traj')):
    c=pd.read_csv(f,sep=r'\s+',comment='#',header=None,names=['step','S1','S2','S3','theta','r'])
    tag=f.split('/')[-2]; rows.append((tag, c.S1.iloc[-50:].mean(), c.S2.iloc[-50:].mean(), c.S3.iloc[-50:].mean(), c.theta.iloc[-50:].mean()))
df=pd.DataFrame(rows,columns=['run','S1','S2','S3','theta']); df['prot']=df.run.str[-1]
print(df.groupby('prot')[['S1','S2','S3','theta']].describe().loc[:,(slice(None),['mean','min','max'])].round(2).T.to_string())
print(df.to_string(index=False, float_format=lambda x:f"{x:.2f}"))
# structural difference A-runs vs B-runs (last frame of each), per residue CA
ref=md.load(glob.glob(D+'/frame_001__A/unbias.A.dcd')[0],top=top,frame=0)
ca=ref.topology.select('name CA')
last={p:[] for p in 'AB'}
for f in sorted(glob.glob(D+'/*/unbias.[AB].dcd')):
    t=md.load(f,top=top)[-1]; last[f[-5]].append(t)
A=md.join(last['A']); B=md.join(last['B']); allt=md.join([A,B]); allt.superpose(allt[0],atom_indices=ca)
XA=allt.xyz[:20][:,ca]*10; XB=allt.xyz[20:][:,ca]*10
d=np.linalg.norm(XA.mean(0)-XB.mean(0),axis=1)
res=[ref.topology.atom(i).residue for i in ca]
o=np.argsort(d)[::-1][:20]
print('largest A-run vs B-run mean CA differences:',[(f"{res[i].chain.index}:{res[i]}",round(float(d[i]),1)) for i in o])
for ch in range(3):
    m=[k for k,r in enumerate(res) if r.chain.index==ch]; print('chain',ch,'median %.2f p90 %.2f max %.2f'%(np.median(d[m]),np.percentile(d[m],90),d[m].max()))
