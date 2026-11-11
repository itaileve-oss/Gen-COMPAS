"""Build Gen-COMPAS TMD targets from eBDIMS2 path frames (documented deviation).

Why: with clean v2 A/B endpoints the iteration-0 diffusion model still collapsed to the
training mean (all 1000 candidates at lambda ~0.47, pairwise CA RMSD 0.53 A vs 2.8 A in
the data, broken CA-CA geometry). Agreed fallback (user, 2026-11-06): use evenly spaced
eBDIMS2 frames (1TNF -> 7JRA, CA-only, CG elastic network + DIMS) as TMD TARGETS ONLY.
They are never training data; RiteWeight only reweights the post-TMD unbiased MD.

Each target = the occupancy-stage template PDB (same atoms/order as the MD protein) with
the 456 CA moved onto the eBDIMS2 frame and occupancy 1 on CA only (gc_md.py pulls only
occupancy-1 atoms; the RMSD force is superposition-invariant, so frame orientation does
not matter). Frames are chosen evenly in cumulative CA-RMSD arc length along the path.

Usage: python make_ebdims_targets.py TEMPLATE.pdb OUT_DIR N FRAME.pdb [FRAME.pdb ...]
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from gc_cv_check import frame_cvs  # noqa: E402


def read_ca(path):
    out = {}
    for l in Path(path).read_text().splitlines():
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            out[(l[21], int(l[22:26]))] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out


def kabsch_rmsd(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0)
    u, s, vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    R = vt.T @ np.diag([1, 1, d]) @ u.T
    return float(np.sqrt(((P @ R.T - Q) ** 2).sum(1).mean()))


def main():
    template, out_dir, n = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    frames = sys.argv[4:]
    tlines = template.read_text().splitlines()
    keys = [(l[21], int(l[22:26])) for l in tlines if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    cas = [read_ca(f) for f in frames]
    for f, c in zip(frames, cas):
        missing = [k for k in keys if k not in c]
        if missing:
            raise SystemExit(f"ABORT: {f} lacks {len(missing)} template CA, e.g. {missing[:3]}")
    X = [np.array([c[k] for k in keys]) for c in cas]
    step = [0.0] + [kabsch_rmsd(X[i], X[i + 1]) for i in range(len(X) - 1)]
    s = np.cumsum(step); s /= s[-1]
    picks = []
    for target in np.linspace(0, 1, n + 1)[1:]:   # n points, skip the A end itself
        i = int(np.argmin(np.abs(s - target)))
        if i not in picks:
            picks.append(i)
    out_dir.mkdir(parents=True, exist_ok=True)
    idx = {k: i for i, k in enumerate(keys)}
    ca_idx_map = {}  # for CVs: positions in a CA-only array keyed (chain, resid)
    log = ["target\tframe\tarc_s\tS1_CA\tS2_CA\tS3_CA\ttheta_spread"]
    for j, i in enumerate(picks):
        new, a = [], 0
        for l in tlines:
            if l.startswith(("ATOM", "HETATM")):
                is_ca = l[12:16].strip() == "CA"
                if is_ca:
                    x = X[i][idx[(l[21], int(l[22:26]))]]
                    l = l[:30] + f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}" + l[54:]
                l = l[:54] + ("  1.00" if is_ca else "  0.00") + l[60:]
            new.append(l)
        name = f"ebdims_target_{j:02d}.pdb"
        (out_dir / name).write_text("\n".join(new) + "\n")
        ca_idx_map = {k: m for m, k in enumerate(keys)}
        cv = frame_cvs(X[i], ca_idx_map)
        log.append(f"{name}\t{Path(frames[i]).name}\t{s[i]:.3f}\t{cv['S1_CA']:.2f}\t{cv['S2_CA']:.2f}"
                   f"\t{cv['S3_CA']:.2f}\t{cv['theta_spread']:.2f}")
    (out_dir / "ebdims_targets.tsv").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
