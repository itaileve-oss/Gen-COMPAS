import sys, glob, mdtraj as md
sys.argv=[""]; exec(open("/workspace/gc_control/saltbridge.py").read().split("T=")[0])
AO, BO = set(f(A-B)), set(f(B-A))
print("A-only:", sorted(AO)); print("B-only:", sorted(BO))
for it in ("0th","1st"):
    fs=sorted(glob.glob(f"/workspace/gc/Iterations_v4/{it}.Iteration/targets/*.pdb"))
    print(f"== {it} iteration targets ({len(fs)})  [A-only kept / B-only formed]")
    for p in fs:
        E=sb(md.load(p)); print(p.split("/")[-1], len(E&AO),"/",len(AO), " ", len(E&BO),"/",len(BO))
