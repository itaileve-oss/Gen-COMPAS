import numpy as np, json
z=np.load('/workspace/gc/Iterations/0th.Iteration/riteweight/features_internal_zmat.npz')
X=z['X'].astype(np.float64); meta=str(z['meta'])
nT=40; nF=100; traj=np.repeat(np.arange(nT),nF); isA=(traj%2==0)
# feature-wise separation: |meanA-meanB| / pooled std (circular-unsafe but indicative)
mA,mB=X[isA].mean(0),X[~isA].mean(0); sd=np.sqrt((X[isA].var(0)+X[~isA].var(0))/2)+1e-9
sep=np.abs(mA-mB)/sd
o=np.argsort(sep)[::-1]
print('features with separation >5 pooled sd:',(sep>5).sum(),'  >10:',(sep>10).sum(),' of',X.shape[1])
print('top 15 feature idx/sep/meanA/meanB:',[(int(i),round(float(sep[i]),1),round(float(mA[i]),2),round(float(mB[i]),2)) for i in o[:15]])
# nearest neighbour across groups
sub=np.arange(0,4000,8); Xs=X[sub]; a=isA[sub]
D=np.sqrt(((Xs[:,None]-Xs[None])**2).sum(2)); np.fill_diagonal(D,np.inf)
print('A frames: NN in A %.2f, NN in B %.2f'%(np.median(D[a][:,a].min(1)),np.median(D[a][:,~a].min(1))))
print('meta (first 600 chars):',meta[:600])
