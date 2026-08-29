"""Maximal class number over a range of traces, for the last column of Table 1.

Usage:  python3 maxh.py LO HI
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys, json
from scan_parallel import class_reps, gl_reps
lo,hi=int(sys.argv[1]),int(sys.argv[2])
mx=0; at=None
for T in range(lo,hi+1):
    D=T*T-4; n=len(gl_reps(T,D,class_reps(D)))
    if n>mx: mx,at=n,T
json.dump({'lo':lo,'hi':hi,'max':mx,'at':at},open(os.path.join(_R,'data',f'mh_{lo}_{hi}.json'),'w'))
print(lo,hi,'max nGL',mx,'at',at)
