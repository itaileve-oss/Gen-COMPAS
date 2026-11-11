"""Remote: choose ONE real frame per node minimising the summed CA RMSD between consecutive
nodes (chain-B fitted; Viterbi over the candidate pools), interpolate 5 steps between them,
superpose all on chain B, write PDB/DCD/CVs."""
import json, sys, numpy as np, mdtraj as md
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","src")); from gc_cv_check import ca_map, frame_cvs
import os
RW=os.environ.get("RW_DIR","/workspace/gc/Iterations_v4/3rd.Iteration/riteweight")   # riteweight/ dir with diffusion_training.dcd/.psf + frame_weights.csv
top=md.load_psf(f"{RW}/diffusion_training.psf")
for c,n in zip(top.chains,"ABC"): c.chain_id=n
J=json.load(open("nodes_v3.json"))["nodes"]
T=md.load(f"{RW}/diffusion_training.dcd",top=top)
cb=top.select("chainid 1 and name CA"); ca=top.select("name CA")
T.superpose(T,0,atom_indices=cb)
Xca=T.xyz[:,ca]*10
import pandas as pd
dfw=pd.read_csv(f"{RW}/frame_weights.csv"); CV=dfw[["S1","S3"]].values
pools=[n["pool"] for n in J]; LAM=1.5
dev=lambda k: LAM*np.hypot(CV[pools[k],0]-J[k]["S1"],CV[pools[k],1]-J[k]["S3"])
cost=[dev(0)]; back=[]
# start: prefer frames closest to the node CVs -> cost 0; transitions: CA RMSD (after chain-B fit, no refit)
for k in range(1,len(pools)):
    a=Xca[pools[k-1]]; b=Xca[pools[k]]
    D=np.sqrt(((a[:,None]-b[None])**2).sum(-1).mean(-1))   # (na,nb)
    D[np.array(pools[k-1])[:,None]==np.array(pools[k])[None]]=np.inf   # no frame reuse
    tot=cost[-1][:,None]+D+dev(k)[None]; back.append(tot.argmin(0)); cost.append(tot.min(0))
sel=[int(cost[-1].argmin())]
for k in range(len(pools)-1,0,-1): sel.append(int(back[k-1][sel[-1]]))
sel=sel[::-1]; frames=[pools[k][i] for k,i in enumerate(sel)]
steps=[np.sqrt(((Xca[frames[k]]-Xca[frames[k+1]])**2).sum(-1).mean()) for k in range(len(frames)-1)]
print("chosen frames",frames); print("consecutive CA RMSD: mean %.2f max %.2f A"%(np.mean(steps),np.max(steps)))
X=T.xyz[frames]; sub=5; out=[]
for k in range(len(X)-1):
    for a in range(sub): out.append(X[k]+(X[k+1]-X[k])*a/sub)
out.append(X[-1]); M=md.Trajectory(np.array(out,dtype=np.float32),top)
M[0].save_pdb("v3_top.pdb"); M.save_dcd("v3_path.dcd")
idx=ca_map(top); rows=[(i,)+tuple(frame_cvs(x*10,idx)[k] for k in ("S1_CA","S2_CA","S3_CA")) for i,x in enumerate(M.xyz)]
np.savetxt("v3_cvs.csv",np.array(rows),delimiter=",",header="frame,S1,S2,S3",comments="",fmt="%.3f")
json.dump(dict(frames=frames,steps=[float(s) for s in steps]),open("v3_frames.json","w"))
print("frames out",len(out),"start CVs",rows[0][1:],"end",rows[-1][1:])
