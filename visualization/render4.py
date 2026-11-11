"""PyMOL: render the A->B path, top view (down the 3-fold axis) and side view.
Chains A/B/C coloured; S1 (orange) and S3 (green) drawn as C->A centroid dashes, as in the
earlier tnf_S1_S3_endpoints.pml (ECD numbering: S1 54-56+124-126, S3 C 75-77+97-99 -> A 64-66+114-116)."""
import os, sys
from pymol import cmd
L=os.environ.get('MOVIE_DIR','.')   # folder with v3_top.pdb, v3_path.dcd; frames go to MOVIE_DIR/png4
XTAL=os.environ.get('XTAL_DIR','examples/tnf/structures')   # 1TNF/7JRA crystal PDBs
cmd.load(f"{L}/v3_top.pdb","tnf"); cmd.load_traj(f"{L}/v3_path.dcd","tnf",state=1)
n=cmd.count_states("tnf"); print("states",n)
cmd.remove("hydro"); cmd.hide("everything"); cmd.show("cartoon","tnf")
cmd.set("cartoon_transparency",0.0); cmd.set("ray_opaque_background",1); cmd.bg_color("white")
cmd.color("slate","tnf and chain A"); cmd.color("grey70","tnf and chain B"); cmd.color("salmon","tnf and chain C")
cmd.load(f"{XTAL}/1TNF_crystal_symmetric.pdb","g1TNF"); cmd.load(f"{XTAL}/7JRA_crystal_asymmetric.pdb","g7JRA"); cmd.remove("g7JRA and not polymer")
for g in ("g1TNF","g7JRA"):
    cmd.super(f"{g} and chain B and name CA","tnf and chain B and name CA",mobile_state=1,target_state=1)
    cmd.hide("everything",g); cmd.show("cartoon",g)
cmd.color("grey40","g1TNF"); cmd.color("gold","g7JRA")
cmd.set("cartoon_transparency",0.6,"g1TNF"); cmd.set("cartoon_transparency",0.45,"g7JRA"); cmd.set("ray_transparency_oblique",1)
cmd.create("loopC","tnf and chain C and resi 102-109",0,0); cmd.hide("everything","loopC"); cmd.show("cartoon","loopC"); cmd.color("red","loopC"); cmd.set("cartoon_loop_radius",0.9,"loopC")
cmd.hide("cartoon","tnf and chain C and resi 103-108")
cmd.set("cartoon_smooth_loops",1); cmd.set("antialias",int(os.environ.get("AA","2"))); cmd.set("ray_shadows",0); cmd.set("ambient",0.45); cmd.set("specular",0.2)
cmd.set("dash_width",5); cmd.set("dash_gap",0.2); cmd.set("dash_radius",0.25); cmd.set("label_size",0); cmd.set("sphere_scale",0.9)
sels={"S1":("chain C and resi 54-56+124-126 and name CA","chain A and resi 54-56+124-126 and name CA","orange"),
      "S3":("chain C and resi 75-77+97-99 and name CA","chain A and resi 64-66+114-116 and name CA","forest")}
import numpy as np
def draw(state):
    for k,(a,b,col) in sels.items():
        for o in (f"{k}a",f"{k}b",f"d{k}"): cmd.delete(o)
        pa=cmd.get_coords(f"tnf and {a}",state).mean(0); pb=cmd.get_coords(f"tnf and {b}",state).mean(0)
        cmd.pseudoatom(f"{k}a",pos=list(map(float,pa))); cmd.pseudoatom(f"{k}b",pos=list(map(float,pb)))
        cmd.distance(f"d{k}",f"{k}a",f"{k}b"); cmd.hide("labels",f"d{k}"); cmd.color(col,f"d{k}")
        cmd.show("spheres",f"{k}a or {k}b"); cmd.color(col,f"{k}a or {k}b")
cmd.dss("tnf",state=1)
cmd.set("cartoon_highlight_color","-1"); cmd.set("ray_trace_mode",1); cmd.set("ray_trace_color","grey30"); cmd.set("cartoon_discrete_colors",0)
cmd.set_view((1,0,0, 0,1,0, 0,0,1, 0,0,-200, 0,0,0, 150,250,-20))
cmd.zoom("tnf",buffer=4,state=1,complete=1); cmd.clip("slab",200)
topview=cmd.get_view()
import math
m=(cmd.get_coords("tnf and chain C and resi 54-56+124-126+75-77+97-99 and name CA",1).mean(0)+cmd.get_coords("tnf and chain A and resi 54-56+124-126+64-66+114-116 and name CA",1).mean(0))/2
phi=math.degrees(math.atan2(m[1],m[0])); SIDE_SIGN=float(os.environ.get("SIDE_SIGN","1"))
cmd.turn("z",SIDE_SIGN*(-90-phi)); cmd.turn("x",-90); cmd.zoom("tnf",buffer=4,state=1,complete=1); cmd.clip("slab",200); sideview=cmd.get_view()
W,H=900,820
NMAX=int(os.environ.get("NMAX",n))
for s in range(int(os.environ.get("START","1")),NMAX+1,int(os.environ.get("STEP","1"))):
    cmd.frame(s); draw(s)
    cmd.delete("loopC"); cmd.create("loopC","tnf and chain C and resi 102-109",s,1); cmd.hide("everything","loopC"); cmd.show("cartoon","loopC"); cmd.color("red","loopC"); cmd.set("cartoon_loop_radius",0.9,"loopC")
    cmd.set_view(topview); cmd.png(f"{L}/png4/top_{s:04d}.png",width=W,height=H,ray=1)
    cmd.set_view(sideview); cmd.png(f"{L}/png4/side_{s:04d}.png",width=W,height=H,ray=1)
    if s%30==0: print("rendered",s,flush=True)
