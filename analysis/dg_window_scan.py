"""dG(B-A) from +/-w windows (in S1 and S3) around the A and B minima (minima located on the
0.15-A-smoothed weighted histogram, 310 K), for w = 0.1-0.5 A; 200 bootstraps over trajectories."""
import sys, numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
kT=0.0019872*310; W=[0.1,0.2,0.3,0.4,0.5]
e0=np.linspace(17,26,201); e1=np.linspace(7,17,201); c0=.5*(e0[1:]+e0[:-1]); c1=.5*(e1[1:]+e1[:-1]); X,Y=np.meshgrid(c0,c1,indexing='ij')
for f in sys.argv[1:]:
    d=pd.read_csv(f); w0=d.fel_weight.values/d.fel_weight.sum(); S1,S3=d.S1.values,d.S3.values; tid=d.trajectory_id.values; T=np.unique(tid)
    ix=np.clip(np.digitize(S1,e0)-1,0,199); iy=np.clip(np.digitize(S3,e1)-1,0,199)
    def run(w):
        H=np.zeros((200,200)); np.add.at(H,(ix,iy),w); F=-kT*np.log(np.maximum(gaussian_filter(H,(0.15/0.045,0.15/0.05)),1e-15))
        a=np.unravel_index(np.argmin(np.where(X<21,F,np.inf)),F.shape); b=np.unravel_index(np.argmin(np.where(X>21.8,F,np.inf)),F.shape)
        out=[]
        for h in W:
            pa=w[(abs(S1-c0[a[0]])<=h)&(abs(S3-c1[a[1]])<=h)].sum(); pb=w[(abs(S1-c0[b[0]])<=h)&(abs(S3-c1[b[1]])<=h)].sum()
            out.append(-kT*np.log(pb/pa) if pa>0 and pb>0 else np.nan)
        return np.array(out),(c0[a[0]],c1[a[1]]),(c0[b[0]],c1[b[1]])
    base,A,B=run(w0); rng=np.random.default_rng(0); bs=[]
    for _ in range(200):
        cnt=np.bincount(rng.choice(len(T),len(T)),minlength=len(T)); wb=w0*cnt[np.searchsorted(T,tid)]; bs.append(run(wb/wb.sum())[0])
    bs=np.array(bs); nA=[((abs(S1-A[0])<=h)&(abs(S3-A[1])<=h)).sum() for h in W]; nB=[((abs(S1-B[0])<=h)&(abs(S3-B[1])<=h)).sum() for h in W]
    print(f"{f.split('/')[-2]}: A min ({A[0]:.2f}, {A[1]:.2f})  B min ({B[0]:.2f}, {B[1]:.2f})")
    for k,h in enumerate(W): print(f"   +/-{h:.1f} A: dG {base[k]:+.2f} +/- {np.nanstd(bs[:,k]):.2f} kcal/mol   (frames A {nA[k]}, B {nB[k]})")
