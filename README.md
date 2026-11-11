# Gen-COMPAS workflow with OpenMM (TNF-α example)

A practical workflow for running **[Gen-COMPAS](https://github.com/Tangcyu/Gen-COMPAS)**
(generative, committor-guided path sampling; Nature 2026,
[s41586-026-11025-1](https://www.nature.com/articles/s41586-026-11025-1)) on a large
protein system, with **OpenMM instead of NAMD** as the MD engine. It also includes
analysis tools for free energies between states, with bootstrap errors, and for
building transition-path movies and structure comparisons.

> This repository does **not** contain Gen-COMPAS itself. It contains our configs,
> MD driver, a small patch, and analysis/visualization scripts. Install upstream
> Gen-COMPAS (MIT) as described below.

The worked example is the **TNF-α trimer**: symmetric apo state (1TNF) → asymmetric
state stabilized by small-molecule inhibitors (7JRA, ligand removed). See
[`examples/tnf/`](examples/tnf/).

## How Gen-COMPAS works (one paragraph)

Each iteration: (1) a diffusion model is trained on all MD frames so far and generates
new structures; (2) a variational committor network (VCN) is trained with endpoint
labels (state A core → q = 0, state B core → q = 1); (3) generated structures with
q ≈ 0.5 are picked as targets; (4) short targeted-MD (TMD) runs pull from A and from B
to each target, followed by short unbiased MD; (5) **RiteWeight** builds a Markov model
over all unbiased segments and gives equilibrium frame weights, from which the free-energy
landscape is computed. No single trajectory has to cross the barrier. The landscape is
stitched together from many overlapping short trajectories.

## Install

```bash
git clone https://github.com/Tangcyu/Gen-COMPAS external/Gen-COMPAS
cd external/Gen-COMPAS && git checkout 5004c0d   # version this workflow was run with
git apply ../../patches/riteweight_ergodic_trimming.patch   # see "Patch" below
pip install -e .
```

Two conda/micromamba environments, deliberately separate (torch and OpenMM CUDA pins
conflict otherwise); `setup/setup_gc_instance.sh` builds both on a bare GPU box:

- `tnfmd`: OpenMM ≥ 8.1, PDBFixer, mdtraj (runs `src/gc_md.py`)
- `gen-compas`: PyTorch (CUDA) + Gen-COMPAS (diffusion, VCN, RiteWeight)

## Layout

| Path | Contents |
|---|---|
| `src/gc_md.py` | **OpenMM stand-in for NAMD.** Reads Gen-COMPAS's NAMD-style `.conf` files (TMD + unbiased) and writes the DCD/colvars files the workflow expects. Point `NAMD.execution.*.command` at it |
| `src/gc_cv_check.py` | Project any structures/trajectory onto the TNF CVs (S1/S2/S3 for all chain pairs, θ-spread) |
| `src/prep_system.py`, `src/equil_protocol.py`, `src/make_plumed.py` | System prep + staged restrained equilibration of both end states (shared with the metaD repo) |
| `setup/` | Instance setup, endpoint building (`build_endpoints_v2.sh`), Gen-COMPAS topologies (`make_topologies.py`), initial diffusion data (`combine_initial.py`), optional eBDIMS2 TMD targets |
| `workflow/tnf.workflow.yaml` | Final TNF config (v4 run) with the reasoning for every non-default value as comments |
| `workflow/NAMD_inputs/` | TMD/unbiased templates (heavy-atom TMD, 500 ps; 2 ns unbiased) |
| `workflow/*.sh` | Launch / resume / chain-iterations scripts |
| `patches/` | RiteWeight ergodic-trimming patch (+ the script that applies it) |
| `analysis/dg_methods.py` | ΔG(B−A) seven ways (whole basin, split, min-to-min, ±0.3 Å window, F-threshold region, paper 3D state boxes, RiteWeight cores) with 200-fold bootstrap over trajectories |
| `analysis/dg_window_scan.py` | ΔG vs window size (±0.1–0.5 Å) |
| `analysis/snapshots.py` | RMSD of every frame to reference crystals; extracts representative real frames |
| `visualization/` | Minimum free-energy path on the FES → one real frame per node chosen for continuity (Viterbi over candidate pools) → PyMOL render with crystal ghosts → composite movie |
| `dev/` | Exploratory/diagnostic scripts from the TNF project (TMD force-constant scans, cross-pulls, RiteWeight diagnostics, salt-bridge analysis). Not needed to run the workflow |

Paths in the shell scripts assume an instance layout with `WORK=/workspace/gc`. Adjust as needed.

## Patch: RiteWeight ergodic trimming

`build_transition_matrix()` gives clusters with no segment-start weight a self-loop,
which makes them absorbing states with eigenvalue 1. The stationary distribution can then
land on an empty cluster, and all weights collapse to NaN; on TNF this happened in
iteration 0. The patch solves for the stationary distribution only on the occupied,
strongly connected set of clusters (standard MSM ergodic trimming).

## Lessons learned on a large system

1. **Train the iteration-0 diffusion model long enough.** At 50 epochs it collapsed to
   the mean structure; 300 epochs worked.
2. **TMD on all heavy atoms** with the force constant scaled per atom (k_total = k × N_heavy / N_CA).
   Cα-only TMD did not reach the targets.
3. **Put the committor cores on the free-energy minima, not on your prepared endpoint
   structures.** Our prepared B structure relaxed away from its starting point. Once
   the cores were moved to the FES minima (after iteration 3), the B core held 11× more
   weight, and iteration 4's effective sample size tripled.
4. **The VCN gets sharp.** By iteration 5, 2985 of 3000 generated structures had q > 0.9 and only
   8 fell in 0.3–0.7. Increase `num_samples` and/or set `require_n_targets: false`, rather
   than letting the step fail and retry with the same samples.
5. **Disable `torch.compile`** (`TORCHDYNAMO_DISABLE=1`) on boxes without a C++ compiler.
6. **Report several ΔG definitions with bootstrap errors.** Whole-basin and
   minimum-window ΔG answer different questions. A broad, shallow basin can win on
   population while losing on depth.

## License

MIT for the code in this repository (see `LICENSE`). Gen-COMPAS is MIT-licensed by its authors.
