"""Project structures (generated candidates, training frames, eBDIMS2 path) onto the
paper CVs S1/S2/S3 for all three cyclic chain pairs, plus our theta-spread.

Works on anything with CA atoms in ECD numbering (6-157): the matched WT topology
(protein_noh) or CA-only eBDIMS2 frames. Writes a CSV, one row per frame.

    python gc_cv_check.py OUT.csv --top TOP.pdb --traj X.dcd --label gen
    python gc_cv_check.py OUT.csv --pdbs DIMS_MD*.pdb --label ebdims
"""
from __future__ import annotations

import argparse
import csv

import mdtraj as md
import numpy as np

CVS = {  # full-length precursor numbering; ECD = paper - 76
    "S1": ([130, 131, 132, 200, 201, 202], [130, 131, 132, 200, 201, 202]),
    "S2": ([169, 170, 171, 172, 173, 174, 193, 194, 195, 196, 197, 198],
           [136, 137, 138, 139, 140, 192, 193, 194, 195, 226, 227, 228]),
    "S3": ([151, 152, 153, 173, 174, 175], [140, 141, 142, 190, 191, 192]),
}
OFF = 76
PAIRS = [("A", "B"), ("B", "C"), ("C", "A")]


def ca_map(top):
    return {(a.residue.chain.chain_id, a.residue.resSeq): a.index
            for a in top.atoms if a.name == "CA"}


def frame_cvs(x, idx):
    row = {}
    for p, q in PAIRS:
        for cv, (ga, gc) in CVS.items():
            a = [idx[(p, r - OFF)] for r in ga]
            c = [idx[(q, r - OFF)] for r in gc]
            row[f"{cv}_{p}{q}"] = float(np.linalg.norm(x[a].mean(0) - x[c].mean(0)))
    common = sorted(set.intersection(*[{r for (ch, r) in idx if ch == c} for c in "ABC"]))
    sel = {c: np.array([idx[(c, r)] for r in common]) for c in "ABC"}
    angs = []
    for p, q in PAIRS:
        P, Q = x[sel[p]], x[sel[q]]
        u, _, vt = np.linalg.svd((P - P.mean(0)).T @ (Q - Q.mean(0)))
        d = np.sign(np.linalg.det(vt.T @ u.T))
        R = vt.T @ np.diag([1, 1, d]) @ u.T
        angs.append(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
    row["theta_spread"] = float(np.std(angs))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--top")
    ap.add_argument("--traj")
    ap.add_argument("--pdbs", nargs="*")
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    if a.traj:
        t = md.load(a.traj, top=a.top)
        frames = [(i, t.xyz[i] * 10) for i in range(len(t))]
        idx = ca_map(t.topology)
    else:
        frames = []
        for i, f in enumerate(a.pdbs):
            t = md.load(f)
            frames.append((i, t.xyz[0] * 10))
        idx = ca_map(t.topology)
    rows = []
    for i, x in frames:
        r = {"label": a.label, "frame": i}
        r.update(frame_cvs(x, idx))
        rows.append(r)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{a.label}: {len(rows)} frames -> {a.out}")


if __name__ == "__main__":
    main()
