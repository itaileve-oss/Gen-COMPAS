"""OpenMM stand-in for NAMD inside the Gen-COMPAS workflow.

Gen-COMPAS runs MD through a command template whose only required placeholders are
{namd} and {config}; its downstream stages read FILES only. So this script takes
the place of `namd3`: point NAMD.execution.<mode>.command at

    ["{namd}", "/abs/path/src/gc_md.py", "--device", "{device}", "{config}"]

with NAMD.namd_path = the python interpreter of the MD env.

Config format: NAMD-like `key value` lines, '#' comments, keys case-insensitive.
Gen-COMPAS rewrites exactly the `TMDk` and `TMDFile` lines of the TMD template.

    systemXML      /abs/system.xml          serialized OpenMM System (from prep_system.py)
    proteinPDB     /abs/equilibrated_protein.pdb   protein atoms, numbering 6-157
    inputState     /abs/equilibrated.xml    or e.g. TMD.A.state.xml (relative = job dir)
    outputName     TMD.A
    run            25000                    steps
    dcdfreq        2500                     steps between frames (protein-only DCD)
    temperature    310
    timestep       4.0                      fs
    seed           0                        0 = random
    TMD            on
    TMDk           10000.0                  kcal/mol/A^2, NAMD convention: U = k/(2N) (RMSD - RMSD*(t))^2
    TMDFile        output.pdb               protein PDB; atoms with occupancy 1 are targeted
    TMDLastStep    25000                    RMSD* goes linearly from initial RMSD to TMDFinalRMSD
    TMDFinalRMSD   0.0                      A
    cvPair         C,A                      ordered subunit pair for the paper's S1/S2/S3

Outputs in the job dir: <outputName>.dcd (protein only), <outputName>.colvars.traj
(one row per DCD frame: step S1 S2 S3 theta_spread tmd_rmsd, distances in A),
<outputName>.state.xml (to continue from), <outputName>.log.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np
import openmm
import openmm.app as app
from openmm import unit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_plumed import CVS  # noqa: E402  paper CV residue groups, full-length numbering

ECD_OFFSET = 76          # full-length precursor numbering -> soluble-domain numbering
KCAL_A2_TO_KJ_NM2 = 4.184 * 100.0


def read_conf(path: Path) -> dict[str, str]:
    conf = {}
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip().rstrip(";").strip()
        if not line:
            continue
        key, _, val = line.partition(" ")
        conf[key.lower()] = val.strip()
    return conf


def flag(conf, key, default=False) -> bool:
    return conf.get(key, "on" if default else "off").lower() in ("on", "yes", "true", "1")


def resolve(job: Path, p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else job / q


def cv_groups(top: app.Topology, pair: tuple[str, str]):
    """Atom indices of every CV's two CA groups, from the protein topology's numbering."""
    ca = {(r.chain.id, int(r.id)): a.index for r in top.residues() for a in r.atoms() if a.name == "CA"}
    out = {}
    for name, (ga, gc) in CVS.items():
        sides = []
        for chain, groups in ((pair[0], ga), (pair[1], gc)):
            idx = [ca[(chain, n - ECD_OFFSET)] for g in groups for n in g]
            sides.append(np.array(idx))
        out[name] = sides
    chains = sorted({c for c, _ in ca})
    common = sorted(set.intersection(*[{n for c, n in ca if c == ch} for ch in chains]))
    theta_sel = [np.array([ca[(ch, n)] for n in common]) for ch in chains]
    return out, theta_sel


def make_whole(xyz: np.ndarray, box: np.ndarray, chain_slices) -> np.ndarray:
    """Shift each protein chain by box vectors so its centre is nearest chain 0's.

    Chains are separate molecules, so per-molecule wrapping (used when OpenMM writes
    states) can put subunits in different periodic images; the CVs would then read
    ~box length. Rectangular box assumed (our systems are).
    """
    xyz = xyz.copy()
    L = np.diag(box)
    ref = xyz[chain_slices[0]].mean(0)
    for sl in chain_slices[1:]:
        shift = np.round((xyz[sl].mean(0) - ref) / L) * L
        xyz[sl] -= shift
    return xyz


def theta_spread(xyz, sel) -> float:
    angs = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        P, Q = xyz[sel[a]], xyz[sel[b]]
        u, _, vt = np.linalg.svd((P - P.mean(0)).T @ (Q - Q.mean(0)))
        d = np.sign(np.linalg.det(vt.T @ u.T))
        R = vt.T @ np.diag([1, 1, d]) @ u.T
        angs.append(math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(R) - 1) / 2)))))
    return float(np.std(angs))


def kabsch_rmsd(P, Q) -> float:
    P, Q = P - P.mean(0), Q - Q.mean(0)
    u, _, vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(u @ vt))
    R = u @ np.diag([1, 1, d]) @ vt
    return float(np.sqrt(((P @ R - Q) ** 2).sum(1).mean()))


def pin_to_free_block():
    """Pin this process to the first free CPU block in $GC_CPU_BLOCKS (e.g. '0-9;10-19').

    Measured on the 4070 Ti S host (2-socket Xeon, GPU on NUMA node 0): pinning to the
    GPU's socket +35%; two runs on separate blocks 516 ns/day total vs 345 for one.
    A flock per block keeps concurrent Gen-COMPAS jobs from sharing cores; the lock
    is held for the process lifetime (fd kept open) and released on exit.
    """
    import fcntl, os
    blocks = os.environ.get("GC_CPU_BLOCKS")
    if not blocks:
        return None, None
    for i, spec in enumerate(blocks.split(";")):
        fd = open(f"/tmp/gc_cpu_block_{i}.lock", "w")
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            fd.close()
            continue
        cpus = set()
        for part in spec.split(","):
            lo, _, hi = part.partition("-")
            cpus.update(range(int(lo), int(hi or lo) + 1))
        os.sched_setaffinity(0, cpus)
        return spec, fd
    return "none free (unpinned)", None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path)
    ap.add_argument("--device", default="0", help="CUDA device index ('' or 'cpu' = CPU)")
    ap.add_argument("--platform", default=None, help="override: CUDA, OpenCL, CPU")
    args = ap.parse_args()

    job = args.config.resolve().parent
    conf = read_conf(args.config)
    out = conf.get("outputname", args.config.stem)
    log = open(job / f"{out}.log", "w")

    def say(msg):
        print(msg, flush=True)
        log.write(msg + "\n")
        log.flush()

    system = openmm.XmlSerializer.deserialize(resolve(job, conf["systemxml"]).read_text())
    prot = app.PDBFile(str(resolve(job, conf["proteinpdb"])))
    n_prot = prot.topology.getNumAtoms()
    temp = float(conf.get("temperature", 310)) * unit.kelvin
    dt = float(conf.get("timestep", 4.0)) * unit.femtosecond
    nsteps = int(float(conf["run"]))
    dcdfreq = int(float(conf.get("dcdfreq", 2500)))
    pair = tuple(c.strip() for c in conf.get("cvpair", "C,A").split(","))

    groups, theta_sel = cv_groups(prot.topology, pair)
    chain_slices = []
    for c in prot.topology.chains():
        idx = [a.index for a in c.atoms()]
        chain_slices.append(slice(idx[0], idx[-1] + 1))

    # The three chains are separate molecules, and OpenMM (CUDA) periodically re-images
    # each MOLECULE into the box. RMSDForce has no PBC, so once a chain is re-imaged to
    # the far side the TMD force sees a torn-apart trimer and pulls the chains apart
    # (found 2026-11-07: B start sat at the box corner; k=1e5 TMD tore inter-chain COM
    # distances from 21-23 to 25-31 A within 10 ps). A zero-energy bond between chains
    # makes OpenMM treat the trimer as ONE molecule for re-imaging. No energy changes.
    join = openmm.CustomBondForce("0")
    first = [s.start for s in chain_slices]
    for j in first[1:]:
        join.addBond(first[0], j, [])
    system.addForce(join)

    # ---- targeted MD -------------------------------------------------------------
    tmd = flag(conf, "tmd")
    tmd_idx = None
    if tmd:
        tgt = app.PDBFile(str(resolve(job, conf["tmdfile"])))
        if tgt.topology.getNumAtoms() != n_prot:
            raise SystemExit(f"TMDFile has {tgt.topology.getNumAtoms()} atoms, protein has {n_prot}")
        names_t = [a.name for a in tgt.topology.atoms()]
        names_p = [a.name for a in prot.topology.atoms()]
        if names_t != names_p:
            raise SystemExit("TMDFile atom names/order differ from proteinPDB")
        # occupancy column is not kept by OpenMM's reader -> parse it directly
        occ = [float(l[54:60]) for l in resolve(job, conf["tmdfile"]).read_text().splitlines()
               if l.startswith(("ATOM", "HETATM"))]
        tmd_idx = [i for i, o in enumerate(occ) if o > 0.5]
        if not tmd_idx:
            raise SystemExit("TMDFile has no atoms with occupancy 1")
        tgt_xyz = np.array(tgt.positions.value_in_unit(unit.nanometer))
        ref = [openmm.Vec3(0, 0, 0)] * system.getNumParticles()
        for i in tmd_idx:
            ref[i] = openmm.Vec3(*tgt_xyz[i])
        rmsd_force = openmm.RMSDForce(ref, tmd_idx)
        k = float(conf["tmdk"]) * KCAL_A2_TO_KJ_NM2 / len(tmd_idx)
        # NAMD applies TMD force only when RMSD lags behind the schedule (one-sided);
        # with TMDFinalRMSD > initial (steering AWAY) the sign flips.
        cv = openmm.CustomCVForce("0.5*k_tmd*max(0, dir*(rmsd-rmsd0))^2")
        cv.addGlobalParameter("dir", 1.0)
        cv.addCollectiveVariable("rmsd", rmsd_force)
        cv.addGlobalParameter("k_tmd", k)
        cv.addGlobalParameter("rmsd0", 0.0)
        system.addForce(cv)
        last = int(float(conf.get("tmdlaststep", nsteps)))
        final = float(conf.get("tmdfinalrmsd", 0.0)) / 10.0
        say(f"TMD on {len(tmd_idx)} atoms, k={conf['tmdk']} kcal/mol/A^2 (per-atom {k:.1f} kJ/mol/nm^2), "
            f"to {final*10:.2f} A by step {last}")

    integ = openmm.LangevinMiddleIntegrator(temp, 1 / unit.picosecond, dt)
    seed = int(conf.get("seed", 0))
    if seed:
        integ.setRandomNumberSeed(seed)
    plat_name = args.platform or ("CPU" if args.device in ("", "cpu") else "CUDA")
    plat = openmm.Platform.getPlatformByName(plat_name)
    # single: measured 345 vs 274 ns/day for mixed on this host; prep/metaD runs used single too
    precision = conf.get("precision", "single")
    props = {"DeviceIndex": args.device, "Precision": precision} if plat_name == "CUDA" else {}
    cpu_block, _lock = pin_to_free_block()
    ctx = openmm.Context(system, integ, plat, props)

    st_path = resolve(job, conf["inputstate"])
    state = openmm.XmlSerializer.deserialize(st_path.read_text())
    ctx.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    pos = np.array(state.getPositions().value_in_unit(unit.nanometer))
    box = np.array(state.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
    pos[:n_prot] = make_whole(pos[:n_prot], box, chain_slices)  # rigid shifts only
    # Recentre: translate EVERYTHING so the protein COM sits at the box centre (a rigid
    # translation of the whole periodic system changes nothing physical).
    pos += np.diag(box) / 2 - pos[:n_prot].mean(0)
    ctx.setPositions(pos)
    if flag(conf, "resetvelocities"):
        ctx.setVelocitiesToTemperature(temp, seed or np.random.randint(1, 2**31))
    else:
        ctx.setVelocities(state.getVelocities())
    ctx.setTime(0.0)
    say(f"precision {precision}; CPU block {cpu_block}")
    say(f"platform {plat.getName()}; {system.getNumParticles()} atoms; {nsteps} steps; "
        f"frame every {dcdfreq}; CV pair {pair[0]}->{pair[1]}; input {st_path.name}")

    rmsd_init = None
    if tmd:
        rmsd_init = kabsch_rmsd(pos[tmd_idx], tgt_xyz[tmd_idx])
        ctx.setParameter("rmsd0", rmsd_init)
        ctx.setParameter("dir", 1.0 if final <= rmsd_init else -1.0)
        say(f"initial RMSD to target {rmsd_init*10:.2f} A")

    from mdtraj.formats import DCDTrajectoryFile
    dcd = DCDTrajectoryFile(str(job / f"{out}.dcd"), "w")
    colv = open(job / f"{out}.colvars.traj", "w")
    colv.write("#            step               S1               S2               S3     theta_spread         tmd_rmsd\n")

    def frame(step):
        s = ctx.getState(getPositions=True)
        x = np.array(s.getPositions().value_in_unit(unit.nanometer))[:n_prot]
        bx = np.array(s.getPeriodicBoxVectors().value_in_unit(unit.nanometer))
        x = make_whole(x, bx, chain_slices)
        vals = [np.linalg.norm(x[a].mean(0) - x[c].mean(0)) * 10 for a, c in groups.values()]
        th = theta_spread(x, theta_sel)
        r = kabsch_rmsd(x[tmd_idx], tgt_xyz[tmd_idx]) * 10 if tmd else float("nan")
        L = np.diag(bx) * 10
        dcd.write(x * 10, cell_lengths=L, cell_angles=[90.0, 90.0, 90.0])
        colv.write(f"{step:16d} {vals[0]:16.6f} {vals[1]:16.6f} {vals[2]:16.6f} {th:16.6f} {r:16.6f}\n")
        colv.flush()
        return vals, th, r

    t0 = time.time()
    step = 0
    chunk = min(dcdfreq, 100) if tmd else dcdfreq   # TMD: update RMSD* every <=100 steps
    while step < nsteps:
        n = min(chunk, nsteps - step)
        if tmd:
            frac = min(1.0, (step + n) / last)
            ctx.setParameter("rmsd0", rmsd_init + (final - rmsd_init) * frac if step < last else final)
        integ.step(n)
        step += n
        if step % dcdfreq == 0 or step == nsteps:
            vals, th, r = frame(step)
            if step % (dcdfreq * 10) == 0 or step == nsteps:
                rate = step * dt.value_in_unit(unit.nanosecond) / ((time.time() - t0) / 86400)
                say(f"step {step:>9d}  S1 {vals[0]:6.2f} S2 {vals[1]:6.2f} S3 {vals[2]:6.2f}  "
                    f"theta-spread {th:5.2f}  rmsd->target {r:6.2f} A  {rate:6.1f} ns/day")
    dcd.close()
    colv.close()

    # Final state WITHOUT the TMD force's parameters mattering: State holds x, v, box.
    final_state = ctx.getState(getPositions=True, getVelocities=True)
    (job / f"{out}.state.xml").write_text(openmm.XmlSerializer.serialize(final_state))
    say(f"done in {time.time()-t0:.0f} s -> {out}.dcd / .colvars.traj / .state.xml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
