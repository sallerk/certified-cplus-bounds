# Four fifths of the `C_+(A)` gap was slack in the upper bound, and CMS's own dual closes it

Every earlier round of this project pushed the **lower** bound on `C_+(A)`. This
round attacks the **upper** bound, which nobody had touched since 2017. At
`A = 28` the interval that had to contain `C_+(28)` was `[1.09405, 1.20995]`,
about 10.6% wide. It is now `[1.09405, 1.1079]`, about 1.9% wide. **88.9% of
that gap was looseness in the upper bound, not room in the lower bound.**

All new upper bounds are certified in ARB ball arithmetic, the same way the
lower bounds were.

## 1. What `C_+(A)` is

Carneiro–Milinovich–Soundararajan (CMS), *Fourier optimization and prime gaps*,
arXiv:1708.04122, Extremal Problem 2, equation (1.3). Definitions first.

* `F` is an even, continuous, real function on the line with `int |F| < infinity`.
* `Fhat(t) = int e^{-2 pi i x t} F(x) dx` is its Fourier transform.
* `(y)_+` means `max(y, 0)`.
* `||F||_1` means `int |F(x)| dx`.

```
C_+(A) := sup over F   [ F(0) - A int_{|t|>1} (Fhat(t))_+ dt ] / ||F||_1
```

`A` is a fixed number bigger than 1. A larger `C_+(A)` gives a better
conditional bound on gaps between primes (CMS Corollary 4) and between primes
represented by a quadratic form (Chirre–Quesada-Herrera Corollary 5). Since
`C_+` appears in the denominator of those bounds, a smaller **upper** bound on
`C_+` does not improve them — it tells you how much improvement is still
available. That is what this round measures.

## 2. Search first: the dual is CMS's, not new

**Read from the paper before computing.** CMS prove their upper bound
`C_+(A) <= min{1.2/(1 - 0.222/(A-1)), 2}` (their Theorem 2(c.2), eq. 1.8) in
their **Section 4.5**, and they prove it
**by exactly the duality used here**. They take a function `Phi` in
`L^infinity` whose distributional Fourier transform equals `1` on `(-1,1)` and
satisfies `-A <= Phihat(t) - 1 <= 0` everywhere, and conclude
`C_+(A) <= ||Phi||_infinity`. That is the same object as the `m` below, with
`m = Phihat` and `mcheck = Phi`.

So **the dual is not new**. What is new is only that it is *optimised
numerically* instead of being fed one hand-built example.

Nobody appears to have done that. Chirre–Pereira–de Laat (*Primes in arithmetic
progressions and semidefinite programming*, arXiv:2005.02393), Chirre–
Quesada-Herrera (arXiv:2012.07781) and Quesada-Herrera's 2024 survey
(arXiv:2411.05095) all construct test functions `F`, which can only push
`C_+(A)` up from below; their semidefinite programming produces certified
*lower* bounds. No improvement of the upper bound appears anywhere I could find.
See Caveats for how hard that negative was checked.

## 3. The dual, re-derived line by line

Fourier inversion gives `F(0) = int Fhat(t) dt`. So the thing being maximised is

```
J(F) = int Fhat(t) dt - A int_{|t|>1} (Fhat(t))_+ dt.
```

Let `m` be any even real function with

```
m(t) = 1              for |t| <= 1
1 - A <= m(t) <= 1    for |t| >  1.
```

**Claim.** `J(F) <= int Fhat(t) m(t) dt` for every admissible `F`.

*Proof.* The claim is the same as
`int Fhat (1-m) dt <= A int_{|t|>1}(Fhat)_+ dt`. On `|t| <= 1` we have `m = 1`,
so `Fhat(1-m) = 0` there and that stretch contributes nothing to either side.
On `|t| > 1` check the sign of `Fhat(t)`:

* `Fhat(t) > 0`: the right side is `A Fhat(t)`, so we need `1 - m <= A`, that is
  `m >= 1 - A`. That is the lower constraint.
* `Fhat(t) < 0`: the right side is `0`, so we need `Fhat (1-m) <= 0`, that is
  `1 - m >= 0`, that is `m <= 1`. That is the upper constraint.
* `Fhat(t) = 0`: both sides are `0`.

Integrating the pointwise inequality gives the claim. ∎

Next, the multiplication formula. Writing
`mcheck(x) = int e^{2 pi i x t} m(t) dt` (for even real `m` this is the same as
the forward transform),

```
int Fhat(t) m(t) dt = int F(x) mcheck(x) dx <= ||F||_1 * sup_x |mcheck(x)|.
```

Dividing by `||F||_1` and taking the supremum over `F`:

```
C_+(A) <= inf over admissible m of  sup_x |mcheck(x)|.
```

**The mandatory sanity check.** Take `m = 1` on `[-1,1]` and `0` outside. This
is admissible whenever `A >= 1`, because then `1 - A <= 0`. Its transform is
`mcheck(x) = sin(2 pi x)/(pi x)`, whose largest absolute value is `2`, at
`x = 0`. That is exactly CMS's trivial bound `C_+(A) <= 2`. The code reproduces
`2.000000` from this `m`, and the ARB certificate reproduces it as a *proof*.

## 4. Why `1.2/(1 - 0.222/(A-1))` is loose

CMS's `Phi` is written out in their Section 4.5:

```
Psi(x) = sin(2 pi x)/(pi x)
       + 2 sin(a pi x)/(pi x) cos(3 pi x)     a = 0.018
       + 2 sin(b pi x)/(pi x) cos(4 pi x)     b = 0.027
       + 2 sin(c pi x)/(pi x) cos(10 pi x)    c = 0.002
       - 0.888 cos(2 pi x) - 0.01 cos(6 pi x)
```

Transformed, their multiplier is

```
m = 1 on [-1,1]
  + boxes of height 1 and widths 0.018, 0.027, 0.002 at t = +-3/2, +-2, +-5
  - Dirac deltas of mass 0.444 at t = +-1
  - Dirac deltas of mass 0.005 at t = +-3.
```

Evaluating `Psi` on a grid of 40 million points gives `sup |Psi| = 1.196530`,
attained at `x = 0.1987`, with `Psi(0) = 1.196000`. So their `1.2` is their
example's own sup-norm, essentially tight — the looseness is not in the
rounding. It is in two other places.

1. **The deltas are illegal at finite `A`.** A delta is the `A = infinity` limit
   of a spike of depth `1-A`. To make it legal CMS smear each delta into a box,
   which needs width about `0.444/(A-1)`, and they pay the multiplicative factor
   `1/(1 - 0.222/(A-1))` for doing so. That factor is entirely an artefact of
   repairing an `A = infinity` object. Solving the problem *at the actual `A`*
   produces the right spike width for free and pays nothing.
2. **`1.196530` is itself far from optimal.** It is the sup-norm of one
   hand-built function with six terms and three tuned decimals. A linear program
   over 150 free cells certifies `1.1079` at `A = 28`, 7.90% lower.

## 5. The minimax, as a linear program

Take `m` piecewise constant: `m = 1` on `[0,1]`, `m = m_j` on the `j`-th cell of
a partition `1 = tau_0 < tau_1 < ... < tau_N = T`, and `m = 0` beyond `T`
(legal, since `1 - A <= 0 <= 1`). Then, with `c_j = m_j - m_{j+1}`, `m_0 = 1`,
`m_{N+1} = 0`,

```
mcheck(x) = (1/(pi x)) sum_{j=0}^{N} c_j sin(2 pi x tau_j)      for x != 0
mcheck(0) = 2 + 2 sum_j m_j (tau_j - tau_{j-1}).
```

This is **linear in the `m_j`**, so minimising `max_x |mcheck(x)|` is a linear
program: minimise `s` subject to `-s <= mcheck(x_i) <= s` on a grid of `x`, and
`1 - A <= m_j <= 1`. Solved with HiGHS through `scipy.optimize.linprog`.

Two facts from the telescoped form are used throughout:

* `|mcheck(x)| <= V/(pi |x|)`, where `V = sum_j |c_j|` is the total variation of
  `m`. So the search over `x` can stop at `X0 = V/(pi s)`. This is a rigorous
  tail bound, not an observation about decay.
* Because of that, the `x`-grid is grown by **cutting planes**: solve, find every
  local maximum of `|mcheck|` on `[0, X0]` that exceeds `s`, add those points,
  re-solve. A fixed grid fine enough to be safe would be far too large.

One trap, found the hard way: if the initial `x`-grid has fewer points than `m`
has cells, the discrete maximum is meaningless and the LP returns an `m` whose
true sup-norm is nowhere near `s` — observed once as `s = 1.1383` against a true
`1.4508`. The grid is now sized at roughly ten sample points per cell and at
least twelve per oscillation period `1/T`. **The number reported is always the
honestly recomputed sup, never the LP objective.**

## 6. Results

**Bandwidth.** `T` is where the multiplier is cut off. The bound falls as `T`
grows and is still falling at the largest `T` tried, so these numbers are not
the best the method can give. At `A = 28`, uniform cells, `N = 30 T`:

| T | 3 | 4 | 5 |
|---|---|---|---|
| cells `N` | 90 | 120 | 150 |
| `sup \|mcheck\|` | 1.13816270 | 1.11949589 | 1.11426206 |

(`T = 2` with `N = 80` gives `1.16470097`.) The successive drops are `-0.0187`
and `-0.0052`, a ratio of `0.28`. Continuing that ratio geometrically would put
the `T -> infinity` value near `1.112`. That is an extrapolation, not a
computation.

**The shape of the optimal multiplier.** It is `1` on `[-1,1]`, plunges in the
single cell just past `t = 1`, and afterwards wanders between about `-1` and
`+1` out to `T`, touching `+1` exactly on sixteen of them at `A = 28`. At
`T = 5`, `N = 150`,
so cells of width `0.02667`:

| A | 2 | 5 | 15 | 28 | 34.5 |
|---|---|---|---|---|---|
| clamp `1-A` | -1 | -4 | -14 | -27 | -33.5 |
| deepest cell | -1.0000 | -4.0000 | -14.0000 | -23.7400 | -23.7399 |
| clamp active? | yes | yes | yes | no | no |
| certified bound | 1.3102 | 1.1852 | 1.1250 | 1.1079 | 1.1079 |

`A = 28` and `A = 34.5` give the **same** multiplier and the same bound. That
looks like `A`-independence but is a resolution effect: the deep cell has mass
`23.74 x 0.02667 = 0.633`, so the clamp would force width at least
`0.633/(A-1) = 0.0234` at `A = 28`, which is just *below* the cell width
`0.02667`. The constraint therefore stops biting at exactly this resolution. A
finer partition near `t = 1` would make it bite again and would separate the two
rows. See Caveats.

**The brackets.** Every entry in the two bound columns is an ARB certificate,
not a float. The lower bounds are from `EXTREMAL2.md`; the upper bounds are from
this round. `CMS U(A)` is `min{1.2/(1 - 0.222/(A-1)), 2}`.

| A | certified lower | certified upper | CMS U(A) | old bracket | new bracket | slack removed |
|---|---|---|---|---|---|---|
| 2 | 1.2982302204 | 1.3102 | 1.542416 | 18.81% | 0.92% | 95.1% |
| 5 | 1.1543329392 | 1.1852 | 1.270513 | 10.06% | 2.67% | 73.4% |
| 15 | 1.1047549197 | 1.1250 | 1.219335 | 10.37% | 1.83% | 82.3% |
| 28 | 1.0940451159 | 1.1079 | 1.209948 | 10.59% | 1.16% | 88.9% |
| 34.5 | 1.0916893276 | 1.1079 | 1.208005 | 10.65% | 2.08% | 80.5% |

`old bracket` and `new bracket` are the width of the interval as a percentage of
the lower bound. `slack removed` is `(CMS - new upper)/(CMS - lower)`: the
fraction of the old interval this round cut away.

**The number the task asked for.** At `A = 28` the gap was `10.59%` and is now
`1.16%`, so **at least 88.9% of it was slack in the upper bound**. "At least",
because the upper bound is still falling in `T`.

## 7. Certification in ARB

The lower-bound certificates in `EXTREMAL2.md` round every quantity so as to
*understate* the answer. An upper bound needs every rounding to go the other
way. `cplus_dual_cert.py` does that. It reads a multiplier file whose cell edges
and cell values are **exact rationals** and proves a statement of the form
`sup_x |mcheck(x)| <= s`.

**Admissibility is checked exactly.** `1 - A <= m_j <= 1` is verified in
rational arithmetic, so there is no question of a floating-point value sitting a
hair outside the allowed range. Cell values are snapped to multiples of `2^-20`
before the file is written, precisely so this check can be exact. Any `m` obeying
the constraints is admissible, so snapping cannot invalidate the bound; it only
changes the value, which is then re-measured.

**A maximum over a grid is not a maximum.** The whole line is covered in three
pieces, none of which is a grid search.

1. `[0, delta]`. From `cos(u) = 1 - 2 sin^2(u/2)` and `|sin(pi x t)| <= pi x t`,
   ```
   |mcheck(x) - mcheck(0)| <= 4 pi^2 x^2 M2,   M2 = int_0^T |m(t)| t^2 dt,
   ```
   and `mcheck(0) = 2 int_0^T m` and `M2` are both exact rationals. `delta` is
   chosen so this is already below the target.
2. `[delta, X0]`. Adaptive bisection with a **second-order** test. On a cell of
   centre `c` and half-width `h`,
   ```
   |mcheck(x)| <= |mcheck(c)| + |mcheck'(c)| h + (1/2) Mb2 h^2,
   Mb2 = 8 pi^2 M2 >= sup |mcheck''|,
   ```
   with `mcheck(c)` and `mcheck'(c)` computed as ARB balls. A cell whose test
   fails is bisected. Nothing is assumed about where the peaks are. The
   second-order term is what makes this affordable: a plain ball enclosure would
   need cells about `10^4` times smaller.
3. `[X0, infinity)`. `|mcheck(x)| <= V/(pi x)` from the telescoped form, with `V`
   the exact total variation, so `X0 = V/(pi s)` ends the search. This is a
   proof, not an appeal to observed decay.

Two evaluation forms are used, each safe where it is used:

```
mcheck(x)  = 2 sum_j m_j d_j cos(pi x (tau_j + tau_{j-1})) sinc(x d_j)   (safe at x = 0)
mcheck'(x) = -(4 pi / w^2) sum_j m_j [ g(w tau_j) - g(w tau_{j-1}) ],
             w = 2 pi x,  g(z) = sin z - z cos z                         (used for x >= delta)
```

**Cost.** `A = 28`, 150 cells, target `1393/1250 = 1.1079`:
`delta = 1.19e-3`, `X0 = 49.89`, 4022 certified cells, 8043 ARB evaluations,
finest half-width `1.9e-4`, **6.0 seconds**.

**Checks built into every run.**

* `mcheck(0)` is computed twice — once from the cell form, once from the
  telescoped form — and the two exact rationals must be equal.
* The ARB value of `mcheck(0)` must overlap the exact rational.
* `mcheck'(0.7)` from the closed form must agree with a central difference. It
  agrees to fifteen digits.

**Positive control.** Running the certificate on `m = 1_{[-1,1]}`
(`multiplier_trivial.txt`) proves `C_+(A) <= 2.000001`, reproducing CMS's
trivial bound as an actual proof.

**Negative control.** Asked to certify `C_+(28) <= 1.1140`, which is *below* the
multiplier's true sup-norm `1.1142830`, the certificate refuses: it reports an
undecided cell at `x ~ 0.0515` after bisecting to width `9e-11`. A certificate
that could not fail would be worthless.

That negative control also exposed something worth stating: the float sup
reported by the LP driver for the `A = 28` multiplier is `1.1142816`, while a
6-million-point grid finds `1.1142830`. The float scan under-resolves by about
`1.4e-6`. This does not affect any certified number — the target is set about
`2e-4` above the float scan, and the ARB proof would simply fail if the target
were too low — but it is why the float sup is never the quoted number.

## Why

Three things were true at once, and together they explain a 10.6% gap.

* The upper bound had **one data point**. CMS wrote down a single multiplier by
  hand and reported its sup-norm. Nobody optimised over multipliers afterwards.
* The lower bound had **three independent efforts** (CMS, Chirre–Pereira–de Laat
  by semidefinite programming, and this project), which between them moved it by
  well under 1%. A quantity three groups cannot move is probably close to
  something.
* The two facts point the same way. When one side of a bracket has been attacked
  three times and the other side never, the slack is on the side nobody
  attacked. That was the bet, and it paid.

There is also a concrete mechanism. CMS's multiplier lives at `A = infinity`: it
uses Dirac deltas, which no finite `A` allows, and the factor
`1/(1 - 0.222/(A-1))` is the price of repairing an object built for the wrong
problem. The linear program never pays that price, because it picks the legal
spike depth and width itself.

## Novelty

**Not new: the dual.** It is CMS's own Section 4.5. The derivation in Section 3
was done independently and then checked against the paper — a successful
search-first outcome, not a discovery.

**Not new: the certification technique.** ARB ball arithmetic for bounds of this
kind is what Chirre–Quesada-Herrera did and what earlier rounds of this project
did. Adaptive bisection with a second-order enclosure is standard.

**New here: the numbers, and the direction.** As far as I could find, nobody has
optimised the dual multiplier numerically, and no upper bound on `C_+(A)` better
than CMS's has been published. Five certified upper bounds, each between 6.71% and
15.06% below the published one, and each with its multiplier saved as exact
rationals so it can be rechecked on its own.

**Not large in the usual sense.** These bounds do not improve any prime-gap
constant — `C_+` sits in a denominator, so only lower bounds do that. What they
buy is knowing how much is left: at `A = 28`, at most 1.9%.

## Caveats

* **The negative literature claim is a search result, not a proof of absence.** I
  read CMS's Section 4.5 and checked Chirre–Pereira–de Laat, Chirre–
  Quesada-Herrera and the 2024 Quesada-Herrera survey, plus targeted searches. All
  the follow-up work I found produces lower bounds. I did not read every paper
  citing CMS.
* **The dual value need not equal `C_+(A)`.** Everything here is a valid upper
  bound whatever the duality gap. But whether `inf_m sup_x |mcheck|` equals
  `C_+(A)` — that is, whether pushing this computation harder would converge to
  the true constant — is not proved here. If there is a gap, the remaining
  bracket cannot be attributed to either side.
* **These are not the best this method can do.** Three restrictions are still
  active. The bandwidth `T = 5` (the bound was still falling in `T`); the uniform
  cell width (which is why `A = 28` and `A = 34.5` collapse to the same number);
  and the cutting-plane loop, which was capped at 20 iterations. At `A = 5` the
  loop had clearly not converged — its LP objective was `1.1719` while the
  reported honest sup was `1.1851`, which is why that row removes only 73% of the
  gap rather than 89%. An earlier run at `A = 28` with a different grid setting
  reached an uncertified `1.1134975`, below the certified `1.1079` reported here.
* **ARB is trusted.** As in every earlier certificate here, the proof is only as
  good as the library's ball arithmetic for `sin`, `cos` and `sinc_pi`. The
  15-digit agreement between the closed-form derivative and a central difference
  is a check on my formulas, not on ARB.
* **Reading the problem is a human step.** That the ratio bounded is `C_+(A)`
  with CMS's normalisation `e^{-2 pi i x t}` and class `A_+` was checked against
  the paper by eye.

## Reproduce

```
python cplus_dual.py 28 5 150                              # one LP, prints the float sup
python cplus_dual_scan.py T                                # bandwidth scan at A = 28
python cplus_dual_dump.py 28 5 150 multiplier_A28.txt      # exact-rational multiplier
python cplus_dual_cert.py multiplier_A28.txt 1.1079        # the ARB proof (~7 s)
python cplus_dual_cert.py multiplier_trivial.txt 2.000001  # positive control: reproduces 2
python cplus_dual_cert.py multiplier_A28.txt 1.1140        # negative control: must FAIL
python cplus_dual_table.py 5 150 2,5,15,28,34.5            # solve + dump every row
python cplus_dual_certall.py 2,5,15,28,34.5                # certify every row, print the table
```

Files: `cplus_dual.py` (dual, LP, cutting planes), `cplus_dual_dump.py`
(exact-rational multiplier files), `cplus_dual_cert.py` (the ARB certificate),
`cplus_dual_scan.py`, `cplus_dual_table.py`, `cplus_dual_certall.py`, and the
multipliers `multiplier_A{2,5,15,28,34.5}.txt` and `multiplier_trivial.txt`.


## Update: a tighter multiplier at A = 28

The bound at `A = 28` was improved from `1.1144` to **`C_+(28) <= 696/625 = 1.1079`**,
certified by `certificates/cplus_dual_cert.py` on `functions/dual/multiplier_A28_g300.txt`.

The gain came from **300 cells on a geometric grid** (denser near `t = 1`, exponent 2.2)
in place of 150 uniform cells, with the `x`-grid refined in step.

Two things that did **not** work, recorded because both looked promising:

* **Enlarging `T` from 5 to 9 with a coarse `x`-grid.** The linear program minimises
  `max |mcheck(x_i)|` over a finite sample, which is a *relaxation* of the true `sup`.
  Extra freedom then buys a smaller sampled max and a larger true sup. The reported
  figure was `1.1033`; the true sup was `1.1682`, worse than the value it replaced.
  Since `m` is supported on `[0,T]`, `mcheck` has feature scale `~1/(2T)`, and the sup
  is attained far out (`x = 17.76` for the original multiplier) - so the `x`-grid must
  be refined **and extended** whenever `T` grows.
* **A bang-bang multiplier**, as used by Carneiro-Milinovich-Quesada-Herrera-Ramos
  (arXiv:2404.08380, eq. 5.1) for their own extremal problem. Projecting the certified
  multiplier onto its bounds makes it **663% worse** (`1.1143 -> 8.5067`), and 134 of
  150 cells are strictly interior. Their box carries an `e^{pi t}` weight; ours is the
  flat, asymmetric `[1-A, 1] = [-27, 1]`, where saturating produces violent oscillation
  that the sup-norm punishes. Seeding the switch points from the primal's sign changes
  also failed (`1.6613`).

## Update: 300-cell geometric grids at every A

All five certified upper bounds were re-solved with **300 cells on a geometric grid**
(denser near `t = 1`) and **`T = 7`** in place of 150 uniform cells on `[1,5]`, with the
`x`-grid refined in step. `T = 7` won at every `A`, so the earlier fixed `T = 5` was
leaving value on the table throughout.

| A | certified lower | CMS `U(A)` | old UB | **new UB** | exact | bracket old -> new | CMS slack removed |
|---|---|---|---|---|---|---|---|
| 2 | 1.2987479784 | 1.542416 | 1.3102 | **1.3101** | `13101/10000` | 0.882% -> **0.874%** | 95.3% |
| 5 | 1.1545963838 | 1.270513 | 1.1852 | **1.1685** | `2337/2000` | 2.651% -> **1.204%** | 88.0% |
| 15 | 1.1057244615 | 1.219335 | 1.1250 | **1.1189** | `11189/10000` | 1.743% -> **1.192%** | 88.4% |
| 28 | 1.0951787825 | 1.209948 | 1.1144 | **1.1079** | `11079/10000` | 1.755% -> **1.162%** | 88.9% |
| 34.5 | 1.0929515526 | 1.208005 | 1.1136 | **1.1050** | `221/200` | 1.889% -> **1.102%** | 89.5% |

Each is certified by `certificates/cplus_dual_cert.py` on the corresponding
`functions/dual/multiplier_A*_g300.txt`. At `A = 28` the negative control refuses
`1.1077` and names `x = 9.5918`, which is where the supremum is attained.

The improvement is largest where the old bound was loosest: `A = 5` and `A = 34.5`
both had their brackets cut by about half. `A = 2` barely moved, which is consistent
with it already being the tightest row.
