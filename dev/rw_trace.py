import numpy as np, sys
sys.path.insert(0,'/workspace/gc/external/Gen-COMPAS')
from tools.riteweight import assign_clusters_random_centers, build_transition_matrix
from scipy.linalg import eig
X=np.load('/workspace/gc/Iterations/0th.Iteration/riteweight/features_internal_zmat.npz')['X'].astype(np.float64)
L=25; nT=40; nF=100
starts=np.concatenate([np.arange(0,nF-L)+i*nF for i in range(nT)]); ends=starts+L
segtraj=np.repeat(np.arange(nT),nF-L)
for K in (100,50):
    rng=np.random.default_rng(2026); w=np.full(starts.size,1/starts.size)
    print(f"--- K={K}")
    for it in range(1,13):
        lab=assign_clusters_random_centers(X,K,rng)
        T,wI=build_transition_matrix(lab[starts],lab[ends],w,K)
        vals,vecs=eig(T.T); n1=(np.abs(vals-1)<1e-8).sum(); idx=np.argmin(np.abs(vals-1))
        v=np.abs(np.real(vecs[:,idx])); pi=v/v.sum()
        empty=(wI<=1e-15).sum(); support=(pi>1e-12).sum()
        scale=pi/(wI+1e-15); w_new=w*scale[lab[starts]]; s=w_new.sum()
        trajw=np.bincount(segtraj,weights=w,minlength=nT)
        print(f"it {it:2d}: empty clusters {empty:3d}, eig==1 count {n1}, pi support {support:3d}/{K}, w min {w.min():.1e} max {w.max():.1e}, "
              f"trajectories with weight>1e-6: {(trajw>1e-6).sum():2d}/40, sum w_new {s:.2e}")
        if not np.isfinite(s) or s<=0: print('   -> collapse'); break
        w=np.clip(w_new,0,np.inf)/s
