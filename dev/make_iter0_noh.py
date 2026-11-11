"""Heavy-atom DCD = initial 800 frames + every iteration-0 unbiased run (control use only)."""
import glob, mdtraj as md
top_h = '/workspace/gc/gencompas/TOP_files/protein.psf'; top_noh = '/workspace/gc/gencompas/TOP_files/protein_noh.psf'
init = md.load('/workspace/gc/data/initial_noh.dcd', top=top_noh)
ref = md.load(glob.glob('/workspace/gc/Iterations/0th.Iteration/namd/*/unbias.A.dcd')[0], top=top_h, frame=0)
heavy = ref.topology.select('not element H')
parts = [init]
for f in sorted(glob.glob('/workspace/gc/Iterations/0th.Iteration/namd/*/unbias.[AB].dcd')):
    t = md.load(f, top=top_h).atom_slice(heavy)
    t = md.Trajectory(t.xyz, init.topology, unitcell_lengths=t.unitcell_lengths, unitcell_angles=t.unitcell_angles)
    parts.append(t)
j = md.join(parts, check_topology=False)
assert j.n_atoms == init.n_atoms == 3570
j.save_dcd('/workspace/gc_control/iter0_plus_initial_noh.dcd')
print('frames', j.n_frames, '(initial', init.n_frames, '+ iteration-0 unbiased', j.n_frames - init.n_frames, 'from', len(parts) - 1, 'runs)')
