"""Write the PSF/PDB topologies Gen-COMPAS needs, from the matched wild-type protein.

  TOP_files/protein.psf/.pdb      protein WITH H (7134 atoms) = RiteWeight io.topology,
                                  Occupancy reference, and the shim's TMD target format
  TOP_files/protein_noh.psf/.pdb  heavy atoms only = diffusion training topology
Both A and B systems share this exact atom list (verified in prep_system.py).
Run in the gen-compas env:  python gencompas/make_topologies.py runs/A_wt_1TNF/fixed.pdb
"""
import sys
from pathlib import Path
import mdtraj as md
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "external" / "Gen-COMPAS"))
from utils.mdtraj_io import save_psf

src = Path(sys.argv[1]); out = Path(__file__).resolve().parent / "TOP_files"
t = md.load(str(src))
noh = t.atom_slice(t.top.select("not element H"))
for tr, name in ((t, "protein"), (noh, "protein_noh")):
    save_psf(tr.top, out / f"{name}.psf"); tr.save_pdb(str(out / f"{name}.pdb"))
    back = md.load(str(out / f"{name}.pdb"), top=str(out / f"{name}.psf"))
    print(name, tr.n_atoms, "atoms; reload via PSF ok:", back.n_atoms == tr.n_atoms,
          "| chains", [c.chain_id for c in tr.top.chains], "| first/last res", tr.top.residue(0), tr.top.residue(tr.n_residues - 1),
          "| bonds", tr.top.n_bonds)
