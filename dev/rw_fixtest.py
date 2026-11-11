import numpy as np, sys, importlib.util
sys.path.insert(0,'/workspace/gc/external/Gen-COMPAS')
spec=importlib.util.spec_from_file_location('rwp','/workspace/gc_control/riteweight_patched.py'); rwp=importlib.util.module_from_spec(spec); sys.modules["rwp"]=rwp; spec.loader.exec_module(rwp)
X=np.load('/workspace/gc/Iterations/0th.Iteration/riteweight/features_internal_zmat.npz')['X'].astype(np.float64)
L=25; nT=40; nF=100
starts=np.concatenate([np.arange(0,nF-L)+i*nF for i in range(nT)]); ends=starts+L; segtraj=np.repeat(np.arange(nT),nF-L)
for K,seed in ((100,2026),(100,7),(50,2026)):
    r=rwp.riteweight(X,starts,ends,K,50,1e-6,5,30,seed)
    w=r.w_segment; tw=np.bincount(segtraj,weights=w,minlength=nT)
    ess=1/np.sum(w**2)
    print(f"K={K} seed={seed}: finite {np.isfinite(w).all()}, passes {len(r.delta_history)}, last delta {r.delta_history[-1]:.2e}, "
          f"segment ESS {ess:.0f}/{w.size}, trajectories >1% weight {(tw>0.01).sum()}/40, max traj weight {tw.max():.3f}")
    print('   per-trajectory weight (A runs even idx, B odd):', np.round(tw,3).tolist())
