"""One row of Table 1: the distribution of p_min over a decade of traces T.

Usage:  python3 decade_scan.py LO HI
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys, json
from scan_parallel import analyse
lo,hi=int(sys.argv[1]),int(sys.argv[2])
dist={}; maxn=0; nmax_at=None; recs={}
for T in range(lo,hi+1):
    p,n=analyse(T)
    if n<2: continue
    k=str(p) if p else 'gt23'
    dist[k]=dist.get(k,0)+1
    if n>maxn: maxn, nmax_at = n, T
    if p is None or p>=11: recs.setdefault(str(p),[]).append(T)
json.dump({'lo':lo,'hi':hi,'dist':dist,'maxnGL':maxn,'maxnGL_at':nmax_at,
           'recs':{k:v[:20] for k,v in recs.items()}},
          open(os.path.join(_R,'data',f'dec_{lo}_{hi}.json'),'w'))
print(lo,hi,dist,'maxnGL',maxn,'at',nmax_at)
