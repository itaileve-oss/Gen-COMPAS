TNF-alpha A vs B: crystal structures and Gen-COMPAS snapshots (made 2026-11-10)
====================================================================

Open compare_AB.pml in PyMOL from this folder (cd here first, then File > Run Script). Press F1-F5 for the scenes:
  F1  1TNF vs 7JRA (the two crystals)
  F2  1TNF vs simulated A minimum
  F3  7JRA vs simulated B minimum
  F4  7JRA vs closest sampled frame
  F5  simulated A minimum vs simulated B minimum
Everything is superposed on chain B CA. Orange dashes = S1, green = S3 (chain C -> chain A centroids).
All simulation files are single, real MD frames (heavy atoms only); nothing is averaged or interpolated.

File                                   S1 (A)  S3 (A)  CA RMSD to 1TNF / 7JRA   notes
1TNF_crystal_symmetric.pdb             19.90    9.43        -                    apo crystal, protein only
7JRA_crystal_asymmetric.pdb            24.16   14.01        -                    inhibitor VGY kept; chains remapped
                                                                                 C->A, B->B, A->C and residues -76
                                                                                 to match the simulation
sim_A_min_frame6496.pdb                19.93    9.80    1.74 / 2.94              top-weight frame at A minimum
sim_B_min_frame3463.pdb                22.92   12.28    2.40 / 2.25              top-weight frame at B minimum
sim_closest_to_7JRA_frame3733.pdb      23.49   14.12    2.85 / 1.38              rare (weight 1.5e-5), from a
                                                                                 trajectory started near 7JRA
sim_closest_to_1TNF_frame6811.pdb      19.80    9.60    1.57 / 2.83

Source: Gen-COMPAS v4 iteration 3, all 16,000 frames (160 unbiased 2-ns trajectories), RiteWeight weights.
