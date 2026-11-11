"""MFEP on iteration-3 FES: A min -> B min -> 7JRA crystal CVs (24.16, 14.01). 40 nodes.
For each node a candidate pool of real frames (<=0.25 A in S1/S3, else 25 nearest)."""
import json, heapq, numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
kT=0.0019872*310
import sys
d=pd.read_csv(sys.argv[1] if len(sys.argv)>1 else 'frame_weights.csv')   # RiteWeight frame_weights.csv of the iteration to use; S1,S3,w=d.S1.values,d.S3.values,d.fel_weight.values
e0=np.linspace(17,26,201); e1=np.linspace(7,17,201); c0=.5*(e0[1:]+e0[:-1]); c1=.5*(e1[1:]+e1[:-1])
H,_,_=np.histogram2d(S1,S3,bins=[e0,e1],weights=w); F=-kT*np.log(np.maximum(gaussian_filter(H,(0.15/0.045,0.15/0.05)),1e-14)); F-=F.min()
X,Y=np.meshgrid(c0,c1,indexing='ij'); R=np.exp(np.minimum(F,12)/kT)
def ij(s1,s3): return (int(np.argmin(abs(c0-s1))),int(np.argmin(abs(c1-s3))))
def dijk(a,b):
    dist=np.full(F.shape,np.inf); prev={}; dist[a]=0; pq=[(0,a)]
    while pq:
        dc,u=heapq.heappop(pq)
        if u==b: break
        if dc>dist[u]: continue
        for di in (-1,0,1):
            for dj in (-1,0,1):
                v=(u[0]+di,u[1]+dj)
                if (di or dj) and 0<=v[0]<200 and 0<=v[1]<200:
                    nd=dc+np.hypot(di*.045,dj*.05)*0.5*(R[u]+R[v])
                    if nd<dist[v]: dist[v]=nd; prev[v]=u; heapq.heappush(pq,(nd,v))
    p=[b]
    while p[-1]!=a: p.append(prev[p[-1]])
    return p[::-1]
A=np.unravel_index(np.argmin(np.where(X<21,F,np.inf)),F.shape); B=np.unravel_index(np.argmin(np.where(X>21.8,F,np.inf)),F.shape); C=ij(24.16,14.01)
p=dijk(A,B)+dijk(B,C)[1:]; P=np.array([[c0[i],c1[j]] for i,j in p]); FP=np.array([F[i,j] for i,j in p])
s=np.r_[0,np.cumsum(np.hypot(*np.diff(P,axis=0).T))]; N=40; ss=np.linspace(0,s[-1],N)
nodes=np.c_[np.interp(ss,s,P[:,0]),np.interp(ss,s,P[:,1])]; Fn=np.interp(ss,s,FP)
out=[]
for k,(x,y) in enumerate(nodes):
    r=np.hypot(S1-x,S3-y); pool=np.where(r<=0.25)[0]
    if len(pool)<25: pool=np.argsort(r)[:25]
    if len(pool)>60: pool=pool[np.argsort(-w[pool])][:60]
    out.append(dict(node=k,S1=float(x),S3=float(y),F=float(Fn[k]),pool=[int(i) for i in pool]))
print('nodes',N,'path F max %.2f, end F %.2f; B min (%.2f, %.2f)'%(FP.max(),FP[-1],c0[B[0]],c1[B[1]]))
json.dump(dict(nodes=out,path=P.tolist(),bmin=[float(c0[B[0]]),float(c1[B[1]])]),open('nodes_v3.json','w'))
