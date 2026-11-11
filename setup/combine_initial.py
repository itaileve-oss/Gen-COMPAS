"""Combine the A and B 5 ns unbiased runs into Gen-COMPAS's initial diffusion DCD
(heavy atoms, first 1 ns of each discarded as re-equilibration) and print basin
statistics on S1/S3 for RiteWeight.committor_labels."""
from pathlib import Path
import mdtraj as md, numpy as np, pandas as pd
GC = Path("/workspace/gc"); TOP = GC / "gencompas/TOP_files"
skip = 100   # frames of 10 ps = 1 ns
parts = []
for P in "AB":
    t = md.load(str(GC / f"systems/{P}/unbiased.dcd"), top=str(TOP / "protein.psf"))[skip:]
    parts.append(t.atom_slice(t.top.select("not element H")))
    cv = pd.read_csv(GC / f"systems/{P}/unbiased.colvars.traj", sep=r"\s+", comment="#", header=None,
                     names=["step", "S1", "S2", "S3", "theta_spread", "tmd_rmsd"]).iloc[skip:]
    print(P, "frames", t.n_frames, " ".join(f"{c} {cv[c].mean():.2f}+/-{cv[c].std():.2f}"
                                           for c in ["S1", "S2", "S3", "theta_spread"]))
joined = md.join(parts)
(GC / "data").mkdir(exist_ok=True)
joined.save_dcd(str(GC / "data/initial_noh.dcd"))
chk = md.load(str(GC / "data/initial_noh.dcd"), top=str(TOP / "protein_noh.psf"))
print("initial_noh.dcd", chk.n_frames, "frames x", chk.n_atoms, "atoms (expect 3570)")
