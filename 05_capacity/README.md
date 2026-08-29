# The vocabulary $N(q)$: sampling and Good–Turing

These two scripts produce the figures of `rem:vocabsize` [7.1] and the
collision probabilities of Section 5.8.

## Files

| file | role |
|---|---|
| `vocab.g` | GAP. Builds $G_q = \langle \sigma_R, \sigma_L, \sigma_\varepsilon \rangle \le \mathrm{Sym}(\mathbb{F}_q^3)$, computes a verified stabilizer chain, and prints the cycle type of `nsamples` uniform random elements, one per line. |
| `analyse.py` | Python. Exhausts $G_3$ independently of GAP, and from the GAP output computes the number of distinct types, the count $f_1$ of types seen once, the Good–Turing missing mass $f_1/n$, and the collision probability $P_q$. |
| `../data/out_q3.txt.gz`, `../data/out_q5.txt.gz`, `../data/out_q7.txt.gz` | the recorded output of the three runs. |
| `../data/summary.json` | the figures, as computed from those three files. |

## Running

The figures of the paper are those of the three runs recorded in `../data/`.
They are recomputed from the recorded output by

```
python3 analyse.py ../data/out_q3.txt.gz ../data/out_q5.txt.gz ../data/out_q7.txt.gz
```

To repeat the runs themselves,

```
gap -q -b -c 'q:=3;; nsamples:=100000;; seed:=1;;' < vocab.g > out_q3.txt
gap -q -b -c 'q:=5;; nsamples:=600000;; seed:=1;;' < vocab.g > out_q5.txt
gap -q -b -c 'q:=7;; nsamples:=600000;; seed:=1;;' < vocab.g > out_q7.txt
python3 analyse.py out_q3.txt out_q5.txt out_q7.txt --json ../data/summary.json
```

GAP 4.12 and Python 3.11 with NumPy. Building the chain at $q = 7$ takes about
seven seconds. Drawing $6 \times 10^5$ samples takes longer, and the two large
runs are best left to finish on their own.

`analyse.py` runs the $q = 3$ control with no argument, and needs no GAP for
that part.

## Why GAP

SymPy's `schreier_sims` handles $q = 3$ and $q = 5$, but at $q = 7$ it does not
return within 285 seconds. GAP builds the chain in 6.5 seconds. The option
`rec(random := 1000)` in `StabChain` is the verified setting, not the Monte
Carlo one, so `Random(G)` is exactly uniform on $G_q$.

Two samplers were tried first and both are wrong. Random words in the
generators fail because $\sigma_R$, $\sigma_L$, and $\sigma_\varepsilon$ are
all odd permutations, so a word of fixed length lies in the even part; at
$q = 3$ such a sampler reaches only 23 of the 41 cycle types. Product
replacement fails as well, since at $q = 7$ it does not mix, and it
overestimates $P_7$ by a factor of about 500.

## The control at $q = 3$

Two quantities are known exactly at $q = 3$, and `analyse.py` recomputes them
by exhausting the 36864 elements of $G_3$, with an implementation independent
of GAP. If either figure comes out differently, then the sampler is at fault.

| quantity | exact value |
|---|---|
| $\lvert G_3 \rvert$ | 36864 |
| $N(3)$ | 41 |
| $P_3$ | 32014139/679477248 = 0.0471158 |

Against $10^5$ draws, the sampled distribution gave $N(3) = 41$ and
$P_3 = 0.0470747$, with chi-squared 28.2 on 40 degrees of freedom.

## The recorded figures

The three runs in `../data/` were made with GAP 4.12.1 and seed 1, and they
give the figures quoted in `rem:vocabsize` [7.1] and in Section 5.8.

| quantity | $q = 5$ | $q = 7$ |
|---|---|---|
| samples | $6 \times 10^5$ | $6 \times 10^5$ |
| distinct cycle types | 30953 | 575613 |
| types seen once | 5391 | 556290 |
| missing mass | 0.009 | 0.927 |
| $\widehat{P_q}$ | $3.98 \times 10^{-4}$ | $1.81 \times 10^{-7}$ |
| $1/\widehat{P_q}$ | $2.51 \times 10^3$ | $5.53 \times 10^6$ |

A sampled figure depends on the seed and on the version of GAP, since the
stabilizer chain that supplies the uniform distribution is built by the
version in use. The count of distinct types varies by a fraction of a percent
between runs. The exact quantities at $q = 3$ do not vary, and neither do the
conclusions of the paper, which use the counts only as lower bounds for
$N(q)$.
