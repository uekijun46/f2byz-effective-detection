# Uniform sampling of cycle types in G_q = <sigma_R, sigma_L, sigma_eps> <= Sym(F_q^3).
# Usage:  gap -q -b -c 'q:=5;; nsamples:=20000;; seed:=1;;' < vocab.g > out_q5.txt
if not IsBound(q) then q := 5; fi;
if not IsBound(nsamples) then nsamples := 20000; fi;
if not IsBound(seed) then seed := 1; fi;

SetPrintFormattingStatus("*stdout*", false);;
pts := [];;
idx := NewDictionary([0,0,0], true);;
for x in [0..q-1] do for y in [0..q-1] do for z in [0..q-1] do
  Add(pts,[x,y,z]); AddDictionary(idx,[x,y,z],Length(pts));
od; od; od;
n := Length(pts);;
mk := function(f) local Lst,p; Lst:=[];
  for p in pts do Add(Lst, LookupDictionary(idx, f(p))); od; return PermList(Lst); end;;
R  := mk(p -> [p[1], p[3], (p[1]*p[3]-p[2]) mod q]);;
Lg := mk(p -> [p[3], p[2], (p[2]*p[3]-p[1]) mod q]);;
Ep := mk(p -> [p[1], p[2], (p[1]*p[2]-p[3]) mod q]);;
G  := Group(R,Lg,Ep);;

# verified stabilizer chain: random:=1000 is GAP's deterministic/verified setting
SetStabChainMutable(G, StabChain(G, rec(random := 1000)));;
Print("#q ", q, "\n");
Print("#degree ", n, "\n");
Print("#order ", Size(G), "\n");
Print("#baselength ", Length(BaseOfGroup(G)), "\n");
Print("#nsamples ", nsamples, "\n");
Print("#seed ", seed, "\n");

Reset(GlobalMersenneTwister, seed);;
ctype := function(g) local c;
  c := SortedList(List(Cycles(g,[1..n]), Length));
  return JoinStringsWithSeparator(List(c,String), ",");
end;;
for i in [1..nsamples] do
  Print(ctype(Random(G)), "\n");
od;
QUIT;
