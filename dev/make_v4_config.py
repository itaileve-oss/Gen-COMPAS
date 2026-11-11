s=open("/workspace/gc/gencompas/tnf.workflow.yaml").read()
hdr="""# v4 FRESH RUN (2026-11-08, user: "start a fresh run from scratch applying all what we have learned,
# default Gen-COMPAS protocol, 4 iterations"). Changes vs tnf.workflow.yaml, all learned 2026-11-07/08:
#  - iteration-0 diffusion epochs 300 (50 under-trained -> mean collapse; expA 300 matched Trp-cage control)
#  - TMD on ALL HEAVY ATOMS (3570), 500 ps, TMDk 783000 = 1e5 x 3570/456 (same per-atom stiffness, NAMD k/N);
#    heavy selection as in their NANMA/Trialanine examples; templates in NAMD_inputs_v4
#  - RiteWeight patched in the repo (ergodic trimming; shipped code -> NaN)
# Unchanged: nAChR noise schedule 10/10/5/3, epochs 1-3 = 20/10/10 (warm start), 20 targets, 2 ns unbiased.
"""
rep=[('root_dir: "/workspace/gc/Iterations"','root_dir: "/workspace/gc/Iterations_v4"'),
 ("iteration_diffusion_epochs: {0: 50, 1: 20, 2: 10, 3: 10}","iteration_diffusion_epochs: {0: 300, 1: 20, 2: 10, 3: 10}   # v4: 0 -> 300"),
 ('selection: "name CA"          # TMD pulls the 456 CA only','selection: "element != H"     # v4: TMD on all 3570 heavy atoms'),
 ('template_path: "/workspace/gc/gencompas/NAMD_inputs"','template_path: "/workspace/gc/gencompas/NAMD_inputs_v4"'),
 ("tmd_force_constant: 100000.0","tmd_force_constant: 783000.0  # v4: 1e5 x 3570/456 for heavy atoms;")]
for a,b in rep:
    assert s.count(a)==1,a; s=s.replace(a,b)
open("/workspace/gc/gencompas/tnf.v4.workflow.yaml","w").write(hdr+s)
print("ok")
