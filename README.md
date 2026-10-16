# Effective detection of $F_2$-by-$\mathbb{Z}$ groups by their finite quotients

Code and data for the paper of that title.

Jun Ueki, Ochanomizu University.

Every numerical statement of the paper appears in the table below, together
with the script that produces it. Each statement is named by its label, and the
number in the current version follows in brackets. Labels do not change when
the paper is renumbered, so the table stays correct.

## Layout

```
src/              shared modules: class enumeration, Vieta action, scans
01_records/       the scans and the records of Section 5
02_period_law/    the collision statistics
03_congruence/    the congruence levels and the small quotients
04_twoblock/      the blindness law and its refinements
05_capacity/      the vocabulary N(q)
06_fixed_locus/   the fixed loci and the counting formulas
data/             recorded outputs
```

## Conventions

The classes of discriminant $D = T^2 - 4$ are the primitive ones, that is,
those whose multiplicator ring is $\mathbb{Z}[X]/(X^2 - TX + 1)$. A class is
recorded by the exponent vector of its positive word $L^{a_1} R^{a_2} \cdots$,
taken up to rotation and reversal. Cycle types are computed on all of
$X(\mathbb{F}_q) = \mathbb{F}_q^3$, unless a level refinement is named.

The enumeration of classes and the Vieta action are implemented twice, and the
two implementations are checked against each other.

* `src/markoff_rigidity.py`, `src/fastscan.py`, `src/scan_parallel.py` work
  through reduced forms, Gaussian composition, and automorph matrices.
* `src/classes.py`, `src/cyc.py`, `src/prim.py` work through positive words in
  $L$ and $R$, with a primitivity filter.
* `06_fixed_locus/gl2.py` treats an arbitrary element of
  $\mathrm{GL}_2\mathbb{Z}$, by factoring it into the Nielsen generators.

## Claim to script

| Statement | Script | Output |
|---|---|---|
| Table 1, $p_{\min}$ by decade of $T$ | `01_records/decade_scan.py`, `01_records/bigscan.py` | `data/win1.json`, `data/win2.json` |
| Table 1, maximal class number per decade | `01_records/maxh.py` | `data/win1.json` |
| Table 1, first row, second implementation | `src/classes.py`, `src/cyc.py` | — |
| `rem:ramdetect` [5.4], ramified primes detect | `01_records/defcheck.py` | — |
| Table 2, least trace blind at $\{3,5,7\}$ and at $\{3,5,7,11\}$ | `01_records/tier1_scan.py`, `01_records/reproduce_window.py` | `data/win2.json` |
| Section 5.5, failures at inert and at ramified primes, $T = 102$ and $258$ | `01_records/defcheck.py` | — |
| `prop:below24` [5.10], groups of order below 24 | `03_congruence/below24_check.py`, `03_congruence/below24_check_alt.py`, `03_congruence/below24_crosscheck.py` | — |
| `ex:fake` [4.12] and `ex:small` [5.11], homomorphism counts | `03_congruence/hom_check.py`, `03_congruence/s4.py`, `03_congruence/levels.py`, `03_congruence/xg.py` | — |
| `prop:order54` [5.15], the groups $(72,41)$ and $(48,33)$ | `03_congruence/order72_check.py` | — |
| `prop:cayley` [5.16], fixed points on the Cayley level | `06_fixed_locus/gl2.py` | — |
| Table 3, collision rates by period | `02_period_law/collision_law.py`, `02_period_law/blocks.py` | — |
| Section 5.8, order of the difference class | `02_period_law/collision_law.py` | — |
| `thm:family` [6.4], blindness law | `04_twoblock/verify_blindness_law.py` | — |
| `cor:Mp` [6.5] and `cor:order` [6.6], the table of $T_p$ | `04_twoblock/tp_table.py` | — |
| Table 4, stratified vocabulary of the two-block stratum | `04_twoblock/twoblock.py` | — |
| `thm:familyE` [6.17], sign-twisted blindness | `04_twoblock/verify_psl_theorem.py` | — |
| Section 6.4, sharpness of $E(q) \mid M$ | `04_twoblock/psl_threshold_scan.py` | — |
| `prop:twist` [6.16], twisted coprime powers | `04_twoblock/cohom_proof.py` | — |
| `rem:vocabsize` [7.1], the vocabulary $N(q)$ and the missing mass | `05_capacity/vocab.g`, `05_capacity/analyse.py`, see `05_capacity/README.md` | `data/summary.json`, `data/out_q*.txt.gz` |
| Section 7.3, fixed loci for $D = 40$ and $D = 60$ | `06_fixed_locus/fixloc.py` | — |
| Section 7.3, counting formulas and the curve $E$ | `06_fixed_locus/gl2.py` | — |

## Requirements

Python 3.11 with NumPy and SymPy, GAP 4.12, and Singular 4.3. The Python
packages are installed by

```
pip install -r requirements.txt
```

A script may be run from any directory. Each one resolves its imports relative
to its own location.

## Reproducing

```
python3 01_records/decade_scan.py 3 1000     # first row of Table 1
python3 src/classes.py 1400                  # independent class enumeration
python3 04_twoblock/tp_table.py              # the table of T_p
python3 04_twoblock/verify_blindness_law.py  # the blindness law
python3 06_fixed_locus/fixloc.py             # the fixed loci, through Singular
python3 05_capacity/analyse.py data/out_q3.txt.gz data/out_q5.txt.gz data/out_q7.txt.gz
python3 03_congruence/below24_crosscheck.py   # prop:below24, both implementations
python3 03_congruence/order72_check.py       # prop:order54, the groups (72,41) and (48,33)
```

The scans behind Table 1 and Table 2 take several hours in total. They write
JSON checkpoints into `data/`, and a scan that stops may be resumed from the
last checkpoint.

## The check for `prop:below24`

The proposition says that every group of order less than 24 is congruence with
$N = \exp G$. `03_congruence/groups_below_24.py` builds the 25 nonabelian
groups of that range and confirms that they have the right orders, that none
is abelian, and that no two of the same order are isomorphic, so the list
agrees with the classification. `below24_check.py` and `below24_check_alt.py`
then decide congruence, by two different criteria, and
`below24_crosscheck.py` compares them.

Both runs add $S_4$ and $\mathrm{SL}_2(\mathbb{F}_3)$ as controls. The
proposition says that these two are not congruence, so a check that reported
them as congruence would be worthless, whatever it said about the other 25.

What the two runs decide is that $\Gamma(\exp G)$ acts trivially on
$\mathcal N_G$. Definition 3.1 of the paper asks for more, namely that the
action of $\mathrm{GL}_2\mathbb{Z}$ extend to an action of
$\mathrm{GL}_2(\mathbb{Z}/N)$; for the 25 groups, which are metabelian, the
extension is given by the theorem of Chen and Deligne, and the runs are a
direct check of the weaker statement.

## The check for `prop:order54`

The proposition says that $(72,41) = (\mathbb{Z}/3)^2 \rtimes Q_8$ is
congruence of level 12 and that it is the smallest two-generated
non-metabelian congruence group. `03_congruence/order72_check.py` builds the
group, confirms that the image of $\mathrm{SL}_2\mathbb{Z}$ in
$\mathrm{SL}_2(\mathbb{Z}/12) \times \mathrm{Sym}(X_G)$ has order 1152, and
then finds, on each orbit of the group generated by the determinant $\pm1$
subgroup and $\mathrm{Aut}(G)$, the involution that extends the action to
$\mathrm{GL}_2(\mathbb{Z}/12)$. For the group $(48,33) = C_4 \circ
\mathrm{SL}_2(\mathbb{F}_3)$ it exhibits a matrix in $\Gamma(24)$ that acts
nontrivially on $X_G$, so that group is not congruence; the construction is
identified with $(48,33)$ through the orbit data of Table 4 of Chen.

## License

The code is under the MIT license. The data in `data/` are under CC BY 4.0.
