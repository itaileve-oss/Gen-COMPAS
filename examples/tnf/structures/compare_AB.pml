# Run from this folder:  cd examples/tnf/structures  then  pymol compare_AB.pml
# TNF-alpha: crystal structures vs Gen-COMPAS (v4, iteration 3) snapshots. All REAL MD frames (no averaging).
# Everything superposed on chain B CA (chain B barely moves). S1 orange / S3 green = C->A centroid distances (ECD numbering).
#   1TNF   symmetric crystal (apo)                      S1 19.90  S3  9.43
#   7JRA   asymmetric crystal (+ inhibitor VGY)          S1 24.16  S3 14.01
#   A_min  highest-weight frame at the A FES minimum     S1 19.93  S3  9.80   CA RMSD to 1TNF 1.74 / 7JRA 2.94
#   B_min  highest-weight frame at the B FES minimum     S1 22.92  S3 12.28   CA RMSD to 1TNF 2.40 / 7JRA 2.25
#   near7JRA  sampled frame closest to 7JRA              S1 23.49  S3 14.12   CA RMSD to 1TNF 2.85 / 7JRA 1.38 (rare, high F)
# Scenes: F1 crystals | F2 1TNF vs A_min | F3 7JRA vs B_min | F4 7JRA vs near7JRA | F5 A_min vs B_min
reinitialize
load ./1TNF_crystal_symmetric.pdb, x1TNF
load ./7JRA_crystal_asymmetric.pdb, x7JRA
load ./sim_A_min_frame6496.pdb, A_min
load ./sim_B_min_frame3463.pdb, B_min
load ./sim_closest_to_7JRA_frame3733.pdb, near7JRA
load ./sim_closest_to_1TNF_frame6811.pdb, near1TNF
remove hydro
remove solvent
python
from pymol import cmd
objs=["x1TNF","x7JRA","A_min","B_min","near7JRA","near1TNF"]
for o in objs[1:]:
    cmd.super(f"{o} and chain B and name CA", "x1TNF and chain B and name CA")
cmd.dss("polymer"); cmd.hide("everything"); cmd.show("cartoon","polymer"); cmd.show("sticks","x7JRA and resn VGY"); cmd.color("yellow","x7JRA and resn VGY and elem C")
pal={"x1TNF":("grey60","grey60","grey60"),"x7JRA":("tv_blue","grey80","salmon"),"A_min":("lightblue","grey90","lightpink"),
     "B_min":("slate","grey70","red"),"near7JRA":("deepblue","grey50","firebrick"),"near1TNF":("palecyan","grey95","pink")}
for o,(a,b,c) in pal.items():
    cmd.color(a,f"{o} and chain A"); cmd.color(b,f"{o} and chain B"); cmd.color(c,f"{o} and chain C")
cmd.color("grey60","x1TNF")
cmd.set("dash_width",4); cmd.set("dash_gap",0.2); cmd.set("label_size",18); cmd.set("sphere_scale",0.7)
for o in objs:
    for k,(sa,sb),col in [("S1",("chain C and resi 54-56+124-126","chain A and resi 54-56+124-126"),"orange"),
                          ("S3",("chain C and resi 75-77+97-99","chain A and resi 64-66+114-116"),"forest")]:
        cmd.pseudoatom(f"{k}a_{o}",selection=f"{o} and {sa} and name CA"); cmd.pseudoatom(f"{k}b_{o}",selection=f"{o} and {sb} and name CA")
        cmd.distance(f"d{k}_{o}",f"{k}a_{o}",f"{k}b_{o}"); cmd.color(col,f"d{k}_{o}"); cmd.show("spheres",f"{k}a_{o} or {k}b_{o}"); cmd.color(col,f"{k}a_{o} or {k}b_{o}")
cmd.group("1TNF_crystal","x1TNF *_x1TNF"); cmd.group("7JRA_crystal","x7JRA *_x7JRA")
for o in ["A_min","B_min","near7JRA","near1TNF"]: cmd.group(f"sim_{o}",f"{o} *_{o}")
members={"1TNF_crystal":"x1TNF","7JRA_crystal":"x7JRA","sim_A_min":"A_min","sim_B_min":"B_min","sim_near7JRA":"near7JRA","sim_near1TNF":"near1TNF"}
def only(*g):
    cmd.disable("all")
    for x in g:
        o=members[x]; cmd.enable(x); cmd.enable(o)
        for n in cmd.get_names("objects"):
            if n.endswith("_"+o): cmd.enable(n)
only("1TNF_crystal","7JRA_crystal"); cmd.orient("x1TNF"); cmd.scene("F1","store",message="crystals: 1TNF (grey) vs 7JRA (blue/grey/salmon, ligand yellow)")
only("1TNF_crystal","sim_A_min"); cmd.scene("F2","store",message="1TNF vs simulated A minimum")
only("7JRA_crystal","sim_B_min"); cmd.scene("F3","store",message="7JRA vs simulated B minimum")
only("7JRA_crystal","sim_near7JRA"); cmd.scene("F4","store",message="7JRA vs closest sampled frame (1.38 A)")
only("sim_A_min","sim_B_min"); cmd.scene("F5","store",message="simulated A minimum vs B minimum")
cmd.scene("F1","recall")
python end
bg_color white
set ray_shadows, 0
