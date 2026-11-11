"""Patch Gen-COMPAS tools/riteweight.py: stationary distribution over occupied, strongly
connected clusters only (ergodic trimming). Original kept as riteweight.py.orig."""
import shutil, sys
p = sys.argv[1]
src = open(p).read()
old = """        T, wI = build_transition_matrix(start_lab, end_lab, w, n_clusters)
        pi = stationary_distribution(T)
"""
new = """        T, wI = build_transition_matrix(start_lab, end_lab, w, n_clusters)
        pi = stationary_distribution_trimmed(T, wI)
"""
assert src.count(old) == 1, "anchor not found"
helper = '''
def stationary_distribution_trimmed(T, wI, eps=1e-15):
    """[TNF PATCH 2026-11-07] Stationary distribution on the occupied, strongly connected
    part of T only (standard MSM ergodic trimming).

    build_transition_matrix() gives clusters with no segment-start weight a self-loop,
    which makes them ABSORBING states with eigenvalue exactly 1; with several eigenvalues
    at 1, stationary_distribution() may return one supported on an empty cluster, and all
    segment weights collapse to zero (seen on TNF iteration 0: NaN at pass 5 with 100
    clusters, pass 7 with 50). Here empty clusters get pi = 0, transitions into them are
    dropped and rows renormalised, and pi is solved on the largest strongly connected set.
    """
    from scipy.sparse.csgraph import connected_components
    occ = np.where(wI > eps)[0]
    pi = np.zeros(T.shape[0], dtype=np.float64)
    if occ.size == 0:
        raise RuntimeError("No occupied clusters.")
    Ts = T[np.ix_(occ, occ)].copy()
    n, cc = connected_components((Ts > 0).astype(np.int8), directed=True, connection="strong")
    if n > 1:
        sizes = np.bincount(cc, weights=wI[occ])
        keep = np.where(cc == int(np.argmax(sizes)))[0]
        occ, Ts = occ[keep], Ts[np.ix_(keep, keep)]
    rs = Ts.sum(axis=1, keepdims=True)
    Ts = Ts / np.where(rs > 0, rs, 1.0)
    pi[occ] = stationary_distribution(Ts)
    return pi

'''
anchor = "@dataclass\nclass RiteWeightResult"
assert src.count(anchor) == 1
src = src.replace(old, new).replace(anchor, helper + anchor)
shutil.copy(p, p + ".orig")
open(p, "w").write(src)
print("patched", p)
