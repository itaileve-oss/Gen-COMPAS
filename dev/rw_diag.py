import numpy as np, sys
sys.path.insert(0,'/workspace/gc/external/Gen-COMPAS')
from tools.riteweight import assign_clusters_random_centers, build_transition_matrix, stationary_distribution
from scipy.sparse.csgraph import connected_components
X=np.load('/workspace/gc/Iterations/0th.Iteration/riteweight/features_internal_zmat.npz')['X'].astype(np.float64)
print('X finite:',np.isfinite(X).all(),' NaN count',np.isnan(X).sum(),' feature ranges min/max',X.min().round(2),X.max().round(2))
L=25; nT=40; nF=100
starts=np.concatenate([np.arange(0,nF-L)+i*nF for i in range(nT)]); ends=starts+L
traj=np.repeat(np.arange(nT),nF)
for K in (100,50,30,10):
    rng=np.random.default_rng(2026)
    lab=assign_clusters_random_centers(X,K,rng)
    T,wI=build_transition_matrix(lab[starts],lab[ends],np.full(starts.size,1/starts.size),K)
    occ=wI>0
    n,cc=connected_components((T>0).astype(int),directed=True,connection='strong')
    # closed classes = SCCs with no outgoing edges
    closed=[c for c in range(n) if not ((T[cc==c][:, cc!=c])>0).any()]
    trajs_per_cluster=[len(set(traj[lab==k])) for k in range(K)]
    vals=np.linalg.eigvals(T.T); n1=(np.abs(vals-1)<1e-8).sum()
    print(f"K={K}: occupied-start clusters {occ.sum()}/{K}, strong components {n}, CLOSED classes {len(closed)}, eigenvalues==1: {n1}, "
          f"median trajectories per cluster {np.median(trajs_per_cluster):.0f}, clusters with 1 trajectory {sum(t==1 for t in trajs_per_cluster)}")
# how separated are trajectories? nearest-neighbour from another trajectory vs own
from numpy.linalg import norm
sub=np.arange(0,4000,10); Xs=X[sub]; ts=traj[sub]
D=np.sqrt(((Xs[:,None]-Xs[None])**2).sum(2)); np.fill_diagonal(D,np.inf)
own=np.array([D[i][ts==ts[i]].min() for i in range(len(sub))]); other=np.array([D[i][ts!=ts[i]].min() for i in range(len(sub))])
print('NN distance within own trajectory median %.2f | to nearest OTHER trajectory median %.2f'%(np.median(own),np.median(other)))
