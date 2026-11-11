"""Pick diffusion model+noise closest to the Trp-cage control (user target 2026-11-07).
Eligible: A-like>=20%, B-like>=20%, between>=1%, nearest-train median >=0.5 A (not copies).
Rank: lowest CA-CA frac>4.2, then more 'between'. Prints 'DIR NOISE' on last line; exit 1 if none."""
import re, sys
C='/workspace/gc_control'; rows=[]
for d in ['expA_300ep','expB_iter0data']:
    tag=d.split('_')[0]
    for nz in ['10','2']:
        try:
            g=open(f'{C}/{tag}_gate_noise{nz}.txt').read(); c=open(f'{C}/{d}_copy_noise{nz}.txt').read()
        except OSError as e: print('missing',e); continue
        m=re.search(r'between \(0.2<lam<0.8\): (\d+)\s+A-like \(<0.2\): (\d+)\s+B-like \(>=0.8\): (\d+)',g)
        f=re.search(r'gen: CA-CA ([\d.]+)\+/-([\d.]+) frac>4.2 ([\d.]+) \| pairwise CA RMSD median ([\d.]+)',g)
        cm=float(re.search(r'COPY_MEDIAN ([\d.]+)',c).group(1))
        bt,a,b=map(int,m.groups()); n=bt+a+b
        r=dict(dir=d,nz=nz,between=bt/n,A=a/n,B=b/n,bad=float(f.group(3)),pair=float(f.group(4)),copy=cm)
        r['ok']=r['A']>=.2 and r['B']>=.2 and r['between']>=.01 and cm>=0.5
        rows.append(r); print(r)
ok=sorted([r for r in rows if r['ok']],key=lambda r:(r['bad'],-r['between']))
if not ok: print('NO MODEL PASSED'); sys.exit(1)
print('CHOSEN', ok[0]); print(ok[0]['dir'], ok[0]['nz'])
