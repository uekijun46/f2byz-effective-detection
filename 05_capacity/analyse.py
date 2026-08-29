"""Analysis of the cycle-type samples produced by vocab.g, plus the exhaustive
computation at q=3.  Reproduces the figures of Remark 7.1 and Section 5.9."""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys, numpy as np
from collections import Counter
from fractions import Fraction

def read(path):
    """Read a GAP output file, plain or gzipped."""
    if path.endswith('.gz'):
        import gzip
        fh = gzip.open(path, 'rt')
    else:
        fh = open(path)
    meta={}; types=[]
    for line in fh:
        line=line.strip()
        if not line: continue
        if line.startswith('#'):
            parts=line[1:].split(None,1)
            if len(parts)==2: meta[parts[0]]=parts[1]
        else: types.append(line)
    return meta, types

def stats(types):
    cnt=Counter(types); n=len(types); D=len(cnt)
    f1=sum(1 for v in cnt.values() if v==1)
    coll=sum(c*(c-1) for c in cnt.values())
    P=coll/(n*(n-1))
    return dict(n=n, distinct=D, f1=f1, missing_mass=f1/n, P=P,
                invP=(1/P if P>0 else float('inf')))

def exhaust_G3():
    q=3; pts=[(x,y,z) for x in range(q) for y in range(q) for z in range(q)]
    idx={p:i for i,p in enumerate(pts)}
    gens=[np.array([idx[(x,z,(x*z-y)%q)] for (x,y,z) in pts],dtype=np.int32),
          np.array([idx[(z,y,(y*z-x)%q)] for (x,y,z) in pts],dtype=np.int32),
          np.array([idx[(x,y,(x*y-z)%q)] for (x,y,z) in pts],dtype=np.int32)]
    ident=np.arange(27,dtype=np.int32); seen={ident.tobytes()}; fr=[ident]
    while fr:
        nx=[]
        for P in fr:
            for g in gens:
                Q=g[P]; b=Q.tobytes()
                if b not in seen: seen.add(b); nx.append(Q)
        fr=nx
    def ct(P):
        s=bytearray(27); out=[]
        for i in range(27):
            if s[i]: continue
            c=0;j=i
            while not s[j]: s[j]=1;j=int(P[j]);c+=1
            out.append(c)
        return ",".join(map(str,sorted(out)))
    cnt=Counter(ct(np.frombuffer(b,dtype=np.int32)) for b in seen)
    NG=len(seen)
    P3=sum(Fraction(c,NG)**2 for c in cnt.values())
    return NG, cnt, P3

if __name__=="__main__":
    args=[a for a in sys.argv[1:] if a!="--json"]
    jsonpath=None
    if "--json" in sys.argv:
        i=sys.argv.index("--json"); jsonpath=sys.argv[i+1]; args=[a for a in args if a!=jsonpath]
    summary={}
    NG,exact,P3=exhaust_G3()
    print("q=3, exhaustive:  |G_3| = %d,  N(3) = %d,  P_3 = %s = %.7f"
          %(NG,len(exact),P3,float(P3)))
    summary["P3_exact"]=str(P3)
    summary["P3_exact_float"]=float(P3)
    summary["N3_exact"]=len(exact)
    summary["order_G3"]=NG
    for path in args:
        meta,types=read(path)
        s=stats(types)
        summary["q%s"%meta['q']]={"n":s['n'],"distinct":s['distinct'],
                                  "f1":s['f1'],"missing_mass":s['missing_mass'],
                                  "P":s['P'],"invP":s['invP'],"seed":meta.get('seed')}
        print("q=%s  n=%d  distinct=%d  f1=%d  missing mass=%.3f  Phat=%.3g  1/Phat=%.3g"
              %(meta['q'],s['n'],s['distinct'],s['f1'],s['missing_mass'],s['P'],s['invP']))
        if meta['q']=='3':
            obs=Counter(types); n=s['n']
            chi=sum((obs.get(t,0)-n*c/NG)**2/(n*c/NG) for t,c in exact.items())
            df=len(exact)-1
            print("      uniformity self-test against the exhaustive distribution:"
                  " chi2 = %.1f on %d degrees of freedom"%(chi,df))
            print("      P_3 exact = %.7f, sampled = %.7f"%(float(P3),s['P']))
            summary["q3_chi2"]=round(chi,1); summary["q3_df"]=df
    if jsonpath:
        import json
        json.dump(summary,open(jsonpath,"w"),indent=2)
        print("written:",jsonpath)
