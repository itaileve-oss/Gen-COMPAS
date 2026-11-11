# Run from this folder:  cd examples/tnf/structures  then  pymol movie_path.pml
# The 40 real MD frames used in tnf_gencompas_iter3_v4_AtoB_movie_v3.mp4 (state 1 = symmetric side, state 40 = 7JRA-like).
# Superposed on chain B. Ghosts: 1TNF grey, 7JRA gold. Use the movie controls (play / arrow keys) to step through states.
reinitialize
load ./movie_path_40_real_frames.pdb, path
load ./1TNF_crystal_symmetric.pdb, x1TNF
load ./7JRA_crystal_asymmetric.pdb, x7JRA
remove hydro
super x1TNF and chain B and name CA, path and chain B and name CA
super x7JRA and chain B and name CA, path and chain B and name CA
dss path, state=0
dss x1TNF
dss x7JRA
hide everything
show cartoon
show sticks, x7JRA and resn VGY
color yellow, x7JRA and resn VGY and elem C
color slate, path and chain A
color grey70, path and chain B
color salmon, path and chain C
color red, path and chain C and resi 102-109
color grey50, x1TNF
color gold, x7JRA and polymer
set cartoon_transparency, 0.6, x1TNF
set cartoon_transparency, 0.5, x7JRA
bg_color white
orient path
set movie_fps, 6
