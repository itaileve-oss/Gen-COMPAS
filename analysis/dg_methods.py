import sys, pandas as pd, numpy as np
from scipy.ndimage import gaussian_filter, label as cclabel

def watershed(F,mk,mask):
    # steepest descent on 8-neighbour grid; each local minimum -> A (1) if its S1<21.4 else B (2)
    n0,n1=F.shape; G=np.where(mask,F,np.inf); tgt=np.full(F.shape,-1); 
    pad=np.pad(G,1,constant_values=np.inf); best=G.copy(); bi=np.zeros(F.shape,int); bj=np.zeros(F.shape,int)
    I,J=np.indices(F.shape)
    for di in (-1,0,1):
        for dj in (-1,0,1):
            v=pad[1+di:1+di+n0,1+dj:1+dj+n1]; m=v<best; best=np.where(m,v,best); bi=np.where(m,di,bi); bj=np.where(m,dj,bj)
    ni=I+bi; nj=J+bj
    for _ in range(2000):
        ni2=ni[ni,nj]; nj2=nj[ni,nj]
        if np.array_equal(ni2,ni) and np.array_equal(nj2,nj): break
        ni,nj=ni2,nj2
    s1=c0[ni]; out=np.where(s1<21.4,1,2); return np.where(mask,out,0)
kT=0.0019872*310
d=pd.read_csv(sys.argv[1]); w=d.fel_weight.values; w=w/w.sum()
S1,S2,S3=d.S1.values,d.S2.values,d.S3.values
print(f'frames {len(d)}, trajectories {d.trajectory_id.nunique()}, ESS (1/sum w^2) = {1/np.sum(w**2):.0f}, max frame weight {w.max():.4f}')
print('RiteWeight cores: A at S1 %.1f-%.1f, B at S1 %.1f-%.1f, S3 %.1f-%.1f'%(S1[d.state=="A"].min(),S1[d.state=="A"].max(),S1[d.state=="B"].min(),S1[d.state=="B"].max(),S3[d.state=="B"].min(),S3[d.state=="B"].max()))
e0=np.linspace(17,26,201); e1=np.linspace(7,17,201); c0=.5*(e0[1:]+e0[:-1]); c1=.5*(e1[1:]+e1[:-1]); X,Y=np.meshgrid(c0,c1,indexing='ij')
ix=np.clip(np.digitize(S1,e0)-1,0,199); iy=np.clip(np.digitize(S3,e1)-1,0,199)
boxes={'Sym':(19.8,20.2,7.8,8.2,9.1,9.5),'Asym1':(21.8,22.2,9.6,10.0,10.8,11.2),'Asym2':(23.7,24.1,11.1,11.5,11.8,12.2),'Asym3':(24.8,25.2,10.6,11.0,10.8,11.2)}
def methods(w):
    r={}; H=np.zeros((200,200)); np.add.at(H,(ix,iy),w)
    Hs=gaussian_filter(H,(0.15/0.045,0.15/0.05)); F=-kT*np.log(np.maximum(Hs,1e-15)); F-=F.min()
    # basins: watershed on smoothed F, seeded at A and B minima
    mA=np.unravel_index(np.argmin(np.where(X<21,F,np.inf)),F.shape); mB=np.unravel_index(np.argmin(np.where(X>21.8,F,np.inf)),F.shape)
    mk=np.zeros(F.shape,int); mk[mA]=1; mk[mB]=2; ws=watershed(F,mk,mask=Hs>1e-9)
    lab=ws[ix,iy]; r['1 whole-basin (watershed)']=-kT*np.log(w[lab==2].sum()/w[lab==1].sum())
    r['2 split S1<21 | S1>21.8']=-kT*np.log(w[S1>21.8].sum()/w[S1<21].sum())
    r['3 min-to-min (smoothed F)']=F[mB]-F[mA]
    win=lambda m:(np.abs(S1-c0[m[0]])<=.3)&(np.abs(S3-c1[m[1]])<=.3)
    r['4 +/-0.3 A window at minima']=-kT*np.log(w[win(mB)].sum()/w[win(mA)].sum())
    for cut in (1.0,2.0):
        reg=(F<=cut); lab2=ws*reg; L=lab2[ix,iy]
        r[f'5 basin region F<={cut:.0f} kcal']=-kT*np.log(w[L==2].sum()/w[L==1].sum())
    P={k:w[(S1>=a-.3)&(S1<=b+.3)&(S2>=c-.3)&(S2<=e+.3)&(S3>=f-.3)&(S3<=g+.3)].sum() for k,(a,b,c,e,f,g) in boxes.items()}
    for k in ('Asym1','Asym2','Asym3'): r[f'6 paper 3D box (pad 0.3) {k}']=-kT*np.log(P[k]/P['Sym']) if P[k]>0 else np.nan
    r['7 RiteWeight cores B vs A']=-kT*np.log(w[d.state.values=="B"].sum()/w[d.state.values=="A"].sum())
    r['_areaA_F<1']=int(((ws==1)&(F<=1)).sum()); r['_areaB_F<1']=int(((ws==2)&(F<=1)).sum())
    return r
base=methods(w)
rng=np.random.default_rng(0); tid=d.trajectory_id.values; T=np.unique(tid); boots=[]
for b in range(200):
    cnt=np.bincount(rng.choice(len(T),len(T)),minlength=len(T)); wb=w*cnt[np.searchsorted(T,tid)]
    boots.append(methods(wb/wb.sum()))
for k,v in base.items():
    if k.startswith('_'): print(f'   {k}: {v} bins'); continue
    bs=np.array([b[k] for b in boots],float); bs=bs[np.isfinite(bs)]
    print(f'   {k:38s} {v:+6.2f} kcal/mol  (bootstrap over trajectories: {np.mean(bs):+.2f} +/- {np.std(bs):.2f}, n={len(bs)})')
