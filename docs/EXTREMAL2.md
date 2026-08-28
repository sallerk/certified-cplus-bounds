# A candidate improvement to the Chirre–Quesada-Herrera constant C_+(28)

Second run of the `EXTREMAL.md` pipeline. Same extremal problem, different parameter.
`EXTREMAL.md` attacked `C_+(36/11)`, which controls gaps between ordinary primes.
This one attacks `C_+(28)`, which controls gaps between primes represented by a
quadratic form. The extremal problem is the same one, so the linear program
carries over; only the parameter `A` and the choice of basis change.

## The problem (verified against the primary source)

Carneiro–Milinovich–Soundararajan, *Fourier optimization and prime gaps*, arXiv:1708.04122,
**Extremal Problem 2**, equation (1.3). Read from the PDF:

```
Fhat(t) = int e^{-2 pi i x t} F(x) dx

C_+(A) := sup_{F in A_+, F != 0}   (1/||F||_1) ( F(0) - A int_{[-1,1]^c} (Fhat(t))_+ dt )
```

`A_+` is the class of even, continuous, real-valued `F`. `(y)_+` means `max(y, 0)`.
`||F||_1` means `int |F(x)| dx`. There is no sign condition and no support condition.
Every function used below is in `L^1` and has `Fhat` in `L^1`, so it is admissible
under either reading of the membership condition.

**Why A = 28.** Chirre–Quesada-Herrera, *Fourier optimization and quadratic forms*,
arXiv:2012.07781 (Q. J. Math. 73 (2022) 539–577). Their Corollary 3 is a
Brun–Titchmarsh-type bound for a positive definite binary quadratic form `f` of
discriminant `-D`, and 28 is the constant in it. Their Section 7 says plainly:
"to optimize the value of the constant in Theorem 4, we must find `C_+(28)`".
Their **Corollary 5** then reads, on GRH for Hecke L-functions,

```
limsup ( p_{n+1,f} - p_{n,f} ) / ( sqrt(p_{n,f}) log p_{n,f} )  <  2 h(-D) / C_+(28)
```

where `p_{n,f}` is the n-th prime represented by `f` and `h(-D)` is the class number.
They report the constant as `1.837 h(-D)`, which is `2 x 0.91833` rounded up.
So a larger `C_+(28)` directly shrinks that constant.

**Published lower bounds at A = 28** (CQH Table 1, whose three columns are three
different ways of building a test function):

| | value | source |
|---|---|---|
| SDP, `P(x) exp(-x^2)`, deg P <= 82 | 1.0818 | CQH Table 1, computed by D. de Laat |
| SDP, deg P <= 122 | 1.0865 | CQH Table 1, same |
| bandlimited, 3 terms | **1.0889** | CQH eq. (6.10)–(6.11), rigorous in ARB |
| upper bound | 1.2/(1 − 0.222/(A−1)) = 1.20995 | CMS Thm 2(c.2), eq. (1.8) |

The record is the third row. Their function is
`F(x) = H(x/0.98644)` with
`H(x) = cos(2 pi x) [ 68/(1 − 16x^2) + 5/(9 − 16x^2) + 1/(25 − 16x^2) ]`,
found by a greedy search over three integers and one dilation. Their rigorous
statement is `||F||_1 / (numerator) < 0.91833`, i.e. `C_+(28) >= 1.0889332`.

## Method

Put `Fhat` on a symmetric interval `[-B, B]` with `B >= 1`, zero outside, and expand
it in cosine harmonics of half-integer order:

```
Fhat(t) = sum_{m=0}^{K} a_m cos(pi m t / (2B))      on [-B, B],  0 outside
F(x)    = B sum_m a_m [ sinc(m/2 + 2Bx) + sinc(m/2 - 2Bx) ],   sinc(u) = sin(pi u)/(pi u)
F(0)    = 2 B sum_m a_m sinc(m/2)
Fhat(B) = sum_m a_m cos(pi m / 2)    (imposed = 0; otherwise F is not in L^1)
```

Two sub-families sit inside this one basis, and each is by itself complete for even
functions on `[0, B]`:

* **even m** gives `cos(pi k t / B)`, the basis used in `EXTREMAL.md`;
* **odd m** gives `cos(pi (2j-1) t / (2B))`, which is *exactly* the CQH family. Their
  `H` is entire of exponential type `2 pi` because `cos(2 pi c) = 0` at `c = (2j-1)/4`,
  and using `p.v. int e^{-2 pi i x t}/(x-a) dx = -i pi sgn(t) e^{-2 pi i a t}` gives
  `Hhat(t) = (pi/4) sum_j A_j (-1)^{j-1}/(2j-1) cos(pi (2j-1) t/2)` on `[-1,1]`, zero
  outside. Dilating by `lambda` puts it in our basis with `B = 1/lambda` and
  `a_{2j-1} = lambda (pi/4) A_j (-1)^{j-1}/(2j-1)`.

So the published test function is a 3-coefficient member of a family we solve over
with up to 71 coefficients. Using even and odd harmonics together is 2-fold
redundant, the Gram matrix is ill-conditioned, and HiGHS reports numerical failure
at most `B`; so the two families are scanned separately and both are reported.

The objective is homogeneous of degree 0, so fix the numerator and minimise the
denominator:

```
minimise ||F||_1   s.t.   F(0) - 2A int_1^B (Fhat)_+ dt >= 1
```

which is a linear program (`u_i >= |F(x_i)|`, `v_j >= max(Fhat(t_j), 0) >= 0`, plus
two variables `w_1 >= |Fhat'(B)|` and `w_2 >= |Fhat''(B)|` that price the tail of
`||F||_1` past the end of the grid, since `F ~ Fhat'(B) cos(2 pi B x)/(2 pi^2 x^2)`
there). Then `C_+(A) >= 1/||F||_1`. The LP value itself is never quoted: the bound
is always recomputed from the returned coefficients with proper quadrature. That
recomputation deliberately rounds against itself — the penalty is taken on a grid
fine enough to over-estimate rather than under-estimate it, which costs about
`1e-8` in the reported value.

## Result

Scanning `B` at `K = 140` (`cplus_refine.py`). Both bases, twelve bandwidths
(nine of them shown):

| B | 1.10 | 1.14 | 1.18 | 1.22 | **1.26** | 1.30 | 1.40 | 1.50 | 1.75 |
|---|---|---|---|---|---|---|---|---|---|
| even m | 1.0939728 | 1.0940372 | 1.0940396 | 1.0940425 | **1.0940433** | 1.0940408 | 1.0940412 | 1.0940426 | 1.0940335 |
| odd m | 1.0939731 | 1.0940375 | 1.0940388 | 1.0940426 | 1.0940432 | 1.0940408 | 1.0940411 | 1.0940427 | 1.0940335 |

The two bases are built from different functions and agree to six decimals at every
bandwidth. The answer is flat in `B` from 1.14 to 1.75, so the band limit is no
longer the binding restriction.

Convergence in the degree `K`, at `B = 1.1`:

| K | 28 | 40 | 56 | 72 | 100 | 140 | 180 |
|---|---|---|---|---|---|---|---|
| even m | 1.0928862 | 1.0934220 | 1.0939534 | 1.0939666 | 1.0939705 | 1.0939731 | 1.0939732 |
| odd m | — | — | — | 1.0939606 | 1.0939716 | 1.0939734 | 1.0939756 |

Re-evaluated from the coefficients **entirely in mpmath at 50 digits**
(`cplus_mp.py`), with `||F||_1` and the penalty computed by locating every sign
change and integrating each smooth piece separately:

```
B        = 1.26,  K = 140,  71 even harmonics

F(0)     = 1.01053143328256986
penalty  = 0.0103944811444205879
||F||_1  = 0.914164255533725168     (of which 6.9e-9 is the tail past x = 120)

C_+(28) >= 1.0940451303844
   1/value = 0.914039075927925      (against CQH's 0.91833)
```

float64 gives 1.0940451304 for the same coefficients, so the two agree to eleven
digits. The margin over the published record is `5.1e-3`.

**This number has since been certified in interval arithmetic.** The certified
value is `C_+(28) >= 1.09404511`, which is `1.4e-8` below the line above. See
the Certification section. It has since been *beaten*, by a different family of
test functions: the current certified value for this row is
`C_+(28) >= 1.0951787825`. See *The full table, certified*. Everything from here to that section is the
uncertified pipeline that found the function; the certificate is what makes the
claim stand.

| | value | ratio |
|---|---|---|
| CQH SDP deg 122 | 1.0865 | **1.00695** |
| CQH bandlimited, as they state it | 1.0889 | **1.00473** |
| CQH bandlimited, rigorous form 1/0.91833 | 1.0889332 | **1.00469** |
| CQH's own function, re-evaluated here | 1.0889972 | **1.00463** |
| CMS upper bound | 1.20995 | 0.904 (consistent) |

**Consequence.** In CQH Corollary 5 the constant is `2 h(-D)/C_+(28)`. The
band-limited function of this section certifies `1.8280782 h(-D)`, improving CQH's
`1.837 h(-D)`. The **current** best, from the consolidated table, is
`2/C_+(28) <= 1.8261858538`, so the certified constant is `1.8261859 h(-D)`. For
`f(u,v) = u^2 + 27 v^2`, where `D = 108` and `h(-D) = 3`, CQH's `5.511` becomes
`5.4786` (this section's function alone gives `5.4843`).

The test function is in `testfunction_A28.txt` (71 cosine coefficients and `B`).
Anyone can check the claim from that file alone — that is the whole content of a
lower bound here.

**The same LP against the rest of CQH's Table 1.** For each `A` the published
record is the best of their three columns. `record` is that best; `LP` is the
honest recomputed value of the best function this solver found (`cplus_table.py`).

| A | F82 | F122 | PW | record | LP | ratio | B* |
|---|---|---|---|---|---|---|---|
| 2.0 | 1.2900 | **1.2933** | 1.2417 | 1.2933 | 1.297783 | 1.00347 | 2.35 |
| 4.0 | 1.1653 | **1.1673** | 1.1439 | 1.1673 | 1.174643 | 1.00629 | 2.35 |
| 6.0 | 1.1320 | **1.1339** | 1.1198 | 1.1339 | 1.141443 | 1.00665 | 1.70 |
| 8.0 | 1.1159 | **1.1192** | 1.1091 | 1.1192 | 1.125667 | 1.00578 | 1.57 |
| 14.0 | 1.0959 | **1.1000** | 1.0966 | 1.1000 | 1.106385 | 1.00580 | 1.70 |
| 22.0 | 1.0852 | 1.0905 | **1.0909** | 1.0909 | 1.097404 | 1.00596 | 1.70 |
| 28.0 | 1.0818 | 1.0865 | **1.0889** | 1.0889 | 1.094043 | 1.00472 | 1.70 |
| 34.5 | 1.0796 | 1.0841 | **1.0875** | 1.0875 | 1.091689 | 1.00385 | 1.70 |

The bold entry says which published method held the record. In CQH's own table the
two methods cross at `A = 21` (there F122 and PW are both 1.0914; below that SDP
wins, above it the bandlimited family does), which is the behaviour they predict in
their Section 7 and state as their Conjecture 11.

The `A = 28` row here is a free cross-check on the headline number: this scan chose
`B = 1.70` and `K = 180`, quite different settings from the `B = 1.26`, `K = 140`
chosen by the careful scan, and returned 1.094043 against 1.094045.

Three caveats on this table. The bandwidth grid for it is coarse (five values,
then two refinements), and at `A = 2` and `A = 4` the winner sat at the top of
that grid, so those two rows are probably not the best this method can do. At
small `A` the LP's own grid for `||F||_1` is only about five points per
oscillation of `F`, so the LP is optimising a slightly wrong objective there — the
reported value is still a valid lower bound, because it is recomputed properly,
just not necessarily the best one available. And only the `A = 28` row was scanned
carefully and re-evaluated in mpmath; the rest are float64 only.

## Checks performed

- **Problem statement read from the CMS and CQH PDFs**, not from a summary. Fourier
  normalisation `e^{-2 pi i x t}` matches. The class `A_+` carries no extra conditions.
- **Fejér calibration.** `F(x) = (sin(pi x)/(pi x))^2` has `Fhat(t) = (1 - |t|)_+`,
  is non-negative, and must give exactly 1. Its cosine series has coefficients
  decaying like `1/k^2`; truncated at K = 400 it gives 0.9998277, and the `1.7e-4`
  deficit is the truncation of that series, not the method.
- **CMS `c_0` calibration.** CMS remark (1.7): `H(x) = cos(2 pi x)/(1 - 16 x^2)` has
  `H(0)/||H||_1 = c_0 = 1.07995...`. That function is `a_1 = pi/4`, `B = 1` in this
  basis. Got **1.0799503136**, and an independent mpmath integration of the closed
  form `cos(2 pi x)/(1 - 16 x^2)` agrees to 12 digits.
- **The published function reproduced.** Feeding CQH's `{68, 5, 1}` and
  `lambda = 0.98644` through this code gives **1.0889972**, i.e.
  `||F||_1/(numerator) = 0.9182760`, consistent with and slightly sharper than their
  rigorous `< 0.91833`. Pointwise, this basis and their closed form agree to `1e-14`.
- **Their whole search reproduced, not just their answer.** Restricting to three odd
  harmonics and scanning `B` reproduces their optimum at their dilation:

  | B | 1.005 | 1.010 | **1.0137** | 1.020 | 1.030 | 1.050 |
  |---|---|---|---|---|---|---|
  | 3-term value | 1.0860894 | 1.0884529 | **1.0889502** | 1.0875280 | 1.0801667 | 1.0524971 |

  Their `lambda = 0.98644` is `B = 1.01375`. This also explains where the gain comes
  from: a 3-coefficient family collapses if the band is widened, so their search
  could not leave `B ≈ 1.014`. With 71 coefficients the optimum moves to `B ≈ 1.26`.
- **Gorbachev's constant recovered, unplanned.** At `B = 1` the penalty vanishes and
  the problem becomes CMS's `C(infinity)`, for which CMS (1.5) quote Gorbachev's
  `1.08185 <= C(infinity) <= 1.09769`. The LP at `B = 1` returns **1.0818063**.
- **Regression against `EXTREMAL.md`.** The even basis at `A = 36/11`, `B = 1.7`,
  `K = 56` returns 1.1986220 against the 1.1986243 recorded there.
- **Consistency with the proven upper bound** 1.20995. The value sits below it.
- **Monotonicity, for free.** CMS prove `C_+(A)` is non-increasing in `A`. The eight
  table rows were each solved independently and come out strictly decreasing:
  1.297783, 1.174643, 1.141443, 1.125667, 1.106385, 1.097404, 1.094043, 1.091689.
- **The pipeline is reproducible.** Re-running `cplus_final.py` returns coefficients
  identical to the saved ones (`max |a_new - a_old| = 0`), so HiGHS is deterministic
  here and the quoted function can be rebuilt from scratch.
- **Admissibility.** `Fhat(B) = 1.8e-15` and `Fhat'(B) = -3.0e-16`, so `F` decays like
  `1/x^3`; measured `x^3 |F(x)|` is about `1e-4` at `x = 400` and at `x = 1600`.
- **The tail is computed, not assumed.** Past `x = 120`, `int |F|` is done two ways:
  by replacing `|sin|` with its mean `2/pi` against the exact envelope
  `R = sqrt(S^2 + C^2)`, and by exact piecewise integration on `[120, 240]` plus
  averaging only beyond that. The two differ by `6e-4` of the tail, and the tail is
  `6.9e-9` of `||F||_1`, so the choice moves the answer by about `4e-12`.
- **Discretisation cannot be flattering the answer.** `||F||_1` recomputed on grids
  finer and wider than the LP used is 0.914164255533 at every setting tried
  (`X` from 60 to 1600, Gauss–Legendre degree 24 to 60, scan density 8 to 40 per
  panel). Refining the penalty grid moves it *down* by `8e-8`, i.e. the quoted
  penalty is an over-estimate, which makes the bound conservative.
- **Three independent formulas for `F`, agreeing to 1e-51.** The raw sinc sum; a
  "two-trig" form that uses `sin(pi(k ± 2Bx)) = ±(-1)^k sin(2 pi B x)` and
  `sin(pi(m/2 ± 2Bx)) = (-1)^{(m-1)/2} cos(2 pi B x)` to replace 142 trig calls by
  two; and a fully collapsed form `F = sin(2 pi B x) S(x) + cos(2 pi B x) C(x)` with
  `S`, `C` sums of rational terms. The third has removable poles at `x = k/(2B)`
  (last one at `K/(4B) = 27.8`) where `S` reaches `2e30` while `F` stays at `1e-9`,
  so it is used only for the tail envelope; the second is the workhorse.
- **Two quadrature rules agree.** `int_0^120 |F|` is `0.45708212461821037704` under
  mpmath's adaptive tanh–sinh and, to all twenty printed digits, the same under
  Gauss–Legendre. It is also unchanged when the sign-change scan is refined
  (1440 / 2880 / 5760 points, 296 roots each time) and when the bisection tolerance
  is tightened.
- **A trap that was caught.** The first pointwise decay check probed
  `x = 25, 50, 100, 200, 400`. With `B = 1.26 = 63/50` every one of those is an exact
  multiple of `1/(2B)`, hence an exact zero of `F`, and the check appeared to show
  `F` collapsing to `1e-48`. Re-probing off the lattice showed ordinary `1/x^3` decay.
  The float64 check had the same blind spot and reported `7e-17` (noise) at all five.

## Certification

Run `python cplus_cert.py testfunction_A28.txt`. It takes about 3.5 seconds.

**The certified result.**

```
C_+(28) >= 1.0940451158772697          (rigorous)
record  =  1.0889331721712239 = 1/0.91833   (CQH, also rigorous)
margin  =  +0.4694%
1/C_+(28) <= 0.9140390880              (CQH: 0.91833)
2/C_+(28) <= 1.8280781761              (CQH: 1.83666, printed as 1.837 h(-D))
```

**What "certified" means.** Every number that feeds the answer is an *interval*,
not a single value. An interval is a pair `[lo, hi]` that is guaranteed to
contain the true value. Adding, multiplying and taking sines of intervals gives
intervals that still contain the true answer, with the rounding error absorbed
into the width. So the final `[lo, hi]` is a proof, not an estimate, and we
report `lo`. The tool is **ARB**, through the `python-flint` package, at 60
working digits. ARB is the same library CQH used to certify the bound this one
displaces.

**Why a lower bound is all that is needed.** `C_+(A)` is a supremum over test
functions. Any admissible `F` gives a valid lower bound. So we do not have to
prove `F` is optimal. We only have to prove two things about this one `F`: that
the top of the fraction is at least some number, and that the bottom is at most
some number.

**Step 1 - exact coefficients.** The 71 numbers in `testfunction_A28.txt` are
`float64` values, so each one is exactly a fraction with a power of two on the
bottom. They are read as exact fractions. One linear condition must hold
*exactly*, or `F` is not in `L^1` and the whole problem does not apply:

```
Fhat(B) = sum_k (-1)^k a_k = 0
```

In the file it is `1.795e-15`, i.e. zero only to rounding. It is forced to zero
by the projection `a_k -> a_k - (-1)^k r/71`, where `r` is that residual. This
is the closest point of the constraint plane, and it moves each coefficient by
at most `2.53e-17`. The second condition, `Fhat'(B) = 0`, needs nothing: in this
basis `Fhat'(B) = -sum_k a_k (pi k/B) sin(pi k)`, and every `sin(pi k)` is zero
because `k` is a whole number. So it is exactly zero for free.

*How much the projection moves the answer:* re-running the whole certification
with the projection switched off (`NOPROJ=1`, a diagnostic only) changes
`int_0^{x_c}|F|` from `0.45708212050602078862` to `0.45708212050602135856`, a
shift of `5.7e-19`, and leaves the penalty unchanged in all 15 printed digits.
So the projection costs nothing. It is not a cosmetic step, though: without it
`Fhat` jumps at `t = B`, `F` decays like `1/x`, and `||F||_1` is infinite.

**Step 2 - the top of the fraction.** `F(0) = 2 B a_0` exactly, because in the
sinc sum every term except `k = 0` vanishes at `x = 0`. With `B = 63/50`,
`F(0) = 1.0105314332825698294`. The subtracted penalty is
`2A int_1^B (Fhat)_+ dt`, where `(y)_+` means `max(y,0)`. It is subtracted, so
it needs an *upper* bound.

**Step 3 - the penalty, from above.** `[1, B]` is swept adaptively. Each
sub-interval gets an interval enclosure of `Fhat`. If the enclosure is entirely
below zero, that piece contributes nothing. If it is entirely above zero, the
piece is integrated in closed form, because `Fhat` is a cosine sum and its
antiderivative is a sine sum. If the enclosure straddles zero, the piece
contributes `width x (upper end of |Fhat|)`, which is an over-estimate and
therefore safe. Only 5 sub-intervals ever straddle. Result:

```
int_1^B (Fhat)_+ dt <= 0.000185615734721797
2A x that           <= 0.0103944811444206
numerator           >= 1.0001369521381492125
```

Two details made this work. First, near a point where `Fhat` almost touches zero
from below, a plain interval enclosure never resolves the sign, because its
width grows like the width of the sub-interval while the function value shrinks
like the *square* of it. The fix is a second-order enclosure,
`Fhat(c) + Fhat'(c)(t-c) + (1/2)Fhat''(interval)(t-c)^2`, which does resolve it.
Second, at `t = B` itself, no enclosure of `Fhat` can ever work, because `Fhat`
and `Fhat'` are both exactly zero there. Instead, Taylor's theorem with those
two exact zeros gives `Fhat(t) = (1/2) Fhat''(xi) (t-B)^2` for some `xi` between
`t` and `B`. So it is enough to certify `Fhat'' < 0` on `[1.246, B]`, which a
sweep does in 437 evaluations, and then `Fhat <= 0` there and that stretch
contributes nothing.

**Step 4 - `||F||_1`, from above.** This is the part that needed the most care.
Write `u = 2Bx`. Two exact facts do the work.

*Fact one: the integral of `F` has a closed form.* `F` is a sum of sinc
functions, and the antiderivative of `sinc` is the sine integral `Si`, which ARB
computes as a certified interval. So

```
int F dx = (1/(2 pi)) sum_k a_k [ Si(pi(u-k)) + Si(pi(u+k)) ]
```

There is therefore **no quadrature error anywhere**. The only thing left to
decide is the sign of `F`, since `int |F| = |int F|` on any stretch where the
sign does not change.

*Fact two: `F` can be written without cancellation.* Because `sum_k (-1)^k a_k`
is exactly zero,

```
F = (2B/(pi u^3)) sin(pi u) W(u),   W(u) = T2 + sum_k (-1)^k a_k k^4/(u^2-k^2)
T2 = sum_k (-1)^k k^2 a_k = 2.8205371847e-03
```

The naive sinc sum loses six or seven digits to cancellation at large `u`, which
made the sweep hopeless. `W` does not, because the exact zero has already been
divided out. Both forms are evaluated at four off-lattice points and agree.

The sweep runs over `u` in `[0, 400]`, split first at half-integers so that no
interval ever contains two of the poles of `W`. Each piece is bisected until its
interval enclosure of `F` has a definite sign, or until its width falls below
`2^-41`. Result: 32678 pieces, 393 sign changes, 797 undecided pieces of total
length `3.6e-10`. Runs of pieces with the same sign are merged and integrated in
closed form. Undecided pieces contribute `width x (upper end of |F|)`.

Note that this needs **no assumption about where the roots are, and no claim
that they were all found**. A missed sign change would show up as a piece whose
enclosure never settles, and would be charged at the over-estimate rate.

*The tail past `u = 400`.* From the `W` form and `|sin| <= 1`,

```
|F(x)| <= (2B/(pi u^3)) |W(u)|,   |W(u)| <= |T2| + sum_k |a_k| k^4/(u^2-k^2)
```

and the right-hand side falls as `u` grows, so for `u >= u_c > 70` it can be
frozen at `u_c`. Integrating `1/u^3`:

```
int_{x_c}^inf |F| dx <= Wb/(2 pi u_c^2),   Wb = 0.0078962809 at u_c = 400
                     <= 7.85458e-9
```

This is the promised `|F| <= K/x^3` with an explicit `K`, derived from the
coefficients and from the exact constraint, not assumed from the observed decay.
Adding everything:

```
int_0^{x_c} |F| dx <= 0.45708212597323277468   (x_c = 158.7302)
tail               <= 7.85458e-9
||F||_1 = 2 x (both) <= 0.91416426765561723464
```

**Cross-checks.**

- The penalty comes out as `0.0103944811444206`. The independent 50-digit mpmath
  pipeline (`cplus_mp.py`, different library, different algorithm) gave
  `0.0103944811444205879`. Sixteen digits.
- Set the cutoff to `u_c = 302`, so `x_c = 119.8413`. The certificate gives
  `int_0^{x_c}|F| = 0.45708212461095745`. mpmath at `x = 120` gave
  `0.45708212461821038`. The gap is `7.3e-12`, and `int_{119.84}^{120}|F|`
  estimated from the measured `1/x^3` decay is about `6e-12`. So the two agree
  to eleven digits and the residue is accounted for.
- The number of sign changes tracks the structure exactly. Above `u = 71`, `F`
  vanishes at every whole number `u`, because every `sinc(u +- k)` in the sum is
  then a sinc of a nonzero whole number. Below it, `F` never vanishes at a whole
  number: `F(k/(2B)) = B a_k`, and no `a_k` is zero. The sweep finds 64 sign
  changes below `u = 71` at every cutoff tried, and then exactly one per whole
  number above it: 295 at `u_c = 302`, 393 at `u_c = 400`. Both equal
  `64 + (u_c - 71)`.
- Independent structural check on that count: clearing denominators,
  `F = (2Bu/pi) sin(pi u) P(u)/D(u)` with `P` a polynomial. Its top coefficient
  is `sum_k (-1)^k a_k = 0` and the next is `T2 != 0`, so `deg P = 138`, `P` is
  even, and `P` can have at most 69 roots with `u > 0`. Observed: 64. Consistent,
  with five roots off the positive real axis. **The certificate does not use
  this**; it is stated because it is a free check that the sweep is not missing
  a large family of sign changes.
- Both directions come out right: the certified `||F||_1` upper bound
  `0.914164268` sits above the uncertified mpmath value `0.914164256`, and the
  certified numerator lower bound sits below the uncertified one.

**Where rigour could still be lost.** Four places, listed worst first.

1. **ARB itself.** The certificate is only as good as the library's ball
   arithmetic, in particular `Si`, `sin_pi`, `cos_pi` and `sinc_pi`. This is the
   same assumption CQH make. It is not independently checked here, beyond the
   16-digit agreement with mpmath noted above.
2. **The reading of the problem.** The certificate proves a statement about a
   ratio. That this ratio is `C_+(28)` as CMS define it, with their Fourier
   normalisation and their class `A_+`, was checked against the PDFs but is a
   human reading, not a machine one.
3. **The projection is part of the claim.** The certified function is *not* the
   raw file; it is the file's 71 float64 values shifted by `-(-1)^k r/71`.
   Anyone reproducing must apply the same shift. The rule is exact and stated
   above, and the shift is `2.5e-17`, but it is a step someone could skip.
4. **The remaining `1.4e-8`.** The certified value is below the uncertified
   `1.0940451303844` by `1.45e-8`, which is the tail bound at `u_c = 400`.
   Raising `u_c` shrinks it. It is a real loss, not an error, and it is 350000
   times smaller than the margin over the record.

None of these is a heuristic step inside the computation. There is no step in
the chain that is "computed numerically and assumed accurate".

## The full table, certified

The section above certifies one number, `C_+(28)`.  This section certifies
the whole published table, twice.  Every row is solved by both pipelines,
both answers are certified in ARB, and the larger of the two is kept.

**What the table is.** Chirre-Quesada-Herrera, arXiv:2012.07781, Table 1,
lists lower bounds for `C_+(A)` at 68 values of `A`.  The values run from
`A = 1` to `A = 34.5` in steps of `0.5`.  For each row this project builds
two test functions of its own, saves the coefficients of each to a file,
and runs an ARB certificate on each file.  So every number in the two
`certified` columns is a rigorous enclosure, not an estimate.

**Novelty, stated first.** The extremal problem is Carneiro-Milinovich-
Soundararajan's.  The constants and the application are Chirre-Quesada-
Herrera's.  The semidefinite program is Chirre-Pereira Junior-de Laat's.
Using ARB ball arithmetic to certify a bound of this kind is what both
groups did.  **New here is only the pair of test functions on each row and
the certificate for each**: 136 coefficient files and 136 rigorous
enclosures.

**The two pipelines, in plain words.**

* **LP** - the band-limited linear program.  `Fhat` is set to zero outside
  a band `[-B, B]` and expanded inside it in cosine harmonics.  Choosing
  the harmonic coefficients is then a linear program, because `||F||_1` and
  the penalty are both maxima of linear functions.  The test function is a
  sum of `sinc`s.  Files: `testfunction_A<A>_r28.txt`.
* **SDP** - the semidefinite program.  `F(x) = P(x) e^{-pi x^2}` with `P` a
  polynomial, and the sign conditions are imposed exactly by a
  sum-of-squares certificate rather than on a grid.  There is no band.  A
  second free parameter `c` dilates the function.  Files:
  `sdpfin_A<A>_d<d>_c<c>.txt`.

The two families are disjoint: nothing band-limited is a Gaussian times a
polynomial.  So the two columns are independent attempts, and taking the
larger one per row is legitimate - any admissible `F` gives a valid lower
bound on a supremum, so the best of two valid bounds is a valid bound.

**Which published number counts as the record.** Their table has three
columns, and all three are lower bounds on the same quantity; they differ
only in which family of test functions was searched.  `F82` and `F122` are
semidefinite programming over `P(x) exp(-x^2)` with `deg P <= 82` and
`<= 122`, computed for them by David de Laat.  `PW` is their band-limited
family, eq. (7.1)-(7.2).  **The record for a row is the largest of the
three, and which column that is changes as `A` grows**: `PW` at `A = 1`,
`F122` for `1.5 <= A <= 20.5`, a dead tie at `A = 21`, and `PW` again for
`21.5 <= A <= 34.5`.  Comparing against a fixed column would invent an
improvement that is not there - at `A = 8.5` `F122` wins (1.1166 against
1.1073) and at `A = 25` `PW` wins (1.0898 against 1.0883).  The `record`
column below is `max(F82, F122, PW)` on every row and the `from` column
says which one it came from.

**The table.** `LP certified` and `SDP certified` are the two ARB lower
bounds.  `best` is the larger of them, `won` says which pipeline produced
it, and `margin` is `best/record - 1` as a percentage.  Bold marks the
winning published column.

| A | F82 | F122 | PW | record | from | LP certified | SDP certified | best | margin | won |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 1.9016 | 1.9307 | **1.9602** | 1.9602 | PW | 1.9797820564 | 1.9988807616 | **1.9988807616** | +1.9733% | SDP |
| 1.5 | 1.4070 | **1.4089** | 1.3430 | 1.4089 | F122 | 1.4127418167 | 1.4124486135 | **1.4127418167** | +0.2727% | LP |
| 2.0 | 1.2900 | **1.2933** | 1.2417 | 1.2933 | F122 | 1.2987479784 | 1.2980200113 | **1.2987479784** | +0.4212% | LP |
| 2.5 | 1.2346 | **1.2378** | 1.1972 | 1.2378 | F122 | 1.2441416449 | 1.2438715799 | **1.2441416449** | +0.5123% | LP |
| 3.0 | 1.2025 | **1.2049** | 1.1719 | 1.2049 | F122 | 1.2117900400 | 1.2112987671 | **1.2117900400** | +0.5718% | LP |
| 3.5 | 1.1807 | **1.1830** | 1.1555 | 1.1830 | F122 | 1.1903272644 | 1.1900058235 | **1.1903272644** | +0.6194% | LP |
| 4.0 | 1.1653 | **1.1673** | 1.1439 | 1.1673 | F122 | 1.1750163119 | 1.1746461035 | **1.1750163119** | +0.6610% | LP |
| 4.5 | 1.1538 | **1.1555** | 1.1355 | 1.1555 | F122 | 1.1635357948 | 1.1630859606 | **1.1635357948** | +0.6954% | LP |
| 5.0 | 1.1448 | **1.1467** | 1.1290 | 1.1467 | F122 | 1.1545963838 | 1.1542845477 | **1.1545963838** | +0.6886% | LP |
| 5.5 | 1.1378 | **1.1396** | 1.1239 | 1.1396 | F122 | 1.1474479961 | 1.1471450899 | **1.1474479961** | +0.6887% | LP |
| 6.0 | 1.1320 | **1.1339** | 1.1198 | 1.1339 | F122 | 1.1415894489 | 1.1412930987 | **1.1415894489** | +0.6781% | LP |
| 6.5 | 1.1271 | **1.1294** | 1.1164 | 1.1294 | F122 | 1.1366961911 | 1.1362111397 | **1.1366961911** | +0.6460% | LP |
| 7.0 | 1.1228 | **1.1255** | 1.1079 | 1.1255 | F122 | 1.1325541012 | 1.1321691877 | **1.1325541012** | +0.6268% | LP |
| 7.5 | 1.1191 | **1.1222** | 1.1112 | 1.1222 | F122 | 1.1290017018 | 1.1285293311 | **1.1290017018** | +0.6061% | LP |
| 8.0 | 1.1159 | **1.1192** | 1.1091 | 1.1192 | F122 | 1.1259220825 | 1.1255066697 | **1.1259220825** | +0.6006% | LP |
| 8.5 | 1.1131 | **1.1166** | 1.1073 | 1.1166 | F122 | 1.1232177322 | 1.1230465945 | **1.1232177322** | +0.5927% | LP |
| 9.0 | 1.1107 | **1.1142** | 1.1058 | 1.1142 | F122 | 1.1208381955 | 1.1206251233 | **1.1208381955** | +0.5958% | LP |
| 9.5 | 1.1086 | **1.1121** | 1.1044 | 1.1121 | F122 | 1.1187127980 | 1.1188316689 | **1.1188316689** | +0.6053% | SDP |
| 10.0 | 1.1067 | **1.1101** | 1.1031 | 1.1101 | F122 | 1.1168110779 | 1.1170610631 | **1.1170610631** | +0.6271% | SDP |
| 10.5 | 1.1049 | **1.1084** | 1.1020 | 1.1084 | F122 | 1.1151023320 | 1.1153827494 | **1.1153827494** | +0.6300% | SDP |
| 11.0 | 1.1033 | **1.1068** | 1.1010 | 1.1068 | F122 | 1.1135467425 | 1.1139221361 | **1.1139221361** | +0.6435% | SDP |
| 11.5 | 1.1019 | **1.1054** | 1.1001 | 1.1054 | F122 | 1.1121356078 | 1.1126215844 | **1.1126215844** | +0.6533% | SDP |
| 12.0 | 1.1005 | **1.1041** | 1.0993 | 1.1041 | F122 | 1.1108483025 | 1.1114469381 | **1.1114469381** | +0.6654% | SDP |
| 12.5 | 1.0992 | **1.1030** | 1.0985 | 1.1030 | F122 | 1.1096665167 | 1.1103482841 | **1.1103482841** | +0.6662% | SDP |
| 13.0 | 1.0980 | **1.1019** | 1.0978 | 1.1019 | F122 | 1.1085760058 | 1.1090110943 | **1.1090110943** | +0.6453% | SDP |
| 13.5 | 1.0969 | **1.1009** | 1.0972 | 1.1009 | F122 | 1.1075637468 | 1.1081577614 | **1.1081577614** | +0.6593% | SDP |
| 14.0 | 1.0959 | **1.1000** | 1.0966 | 1.1000 | F122 | 1.1066355370 | 1.1073745030 | **1.1073745030** | +0.6704% | SDP |
| 14.5 | 1.0949 | **1.0992** | 1.0960 | 1.0992 | F122 | 1.1057706814 | 1.1062637965 | **1.1062637965** | +0.6426% | SDP |
| 15.0 | 1.0940 | **1.0984** | 1.0955 | 1.0984 | F122 | 1.1049629769 | 1.1057244615 | **1.1057244615** | +0.6668% | SDP |
| 15.5 | 1.0931 | **1.0976** | 1.0951 | 1.0976 | F122 | 1.1042083383 | 1.1049996707 | **1.1049996707** | +0.6742% | SDP |
| 16.0 | 1.0922 | **1.0969** | 1.0946 | 1.0969 | F122 | 1.1035017530 | 1.1039744866 | **1.1039744866** | +0.6450% | SDP |
| 16.5 | 1.0915 | **1.0962** | 1.0942 | 1.0962 | F122 | 1.1028371068 | 1.1037076108 | **1.1037076108** | +0.6849% | SDP |
| 17.0 | 1.0907 | **1.0956** | 1.0938 | 1.0956 | F122 | 1.1022128276 | 1.1030960216 | **1.1030960216** | +0.6842% | SDP |
| 17.5 | 1.0900 | **1.0950** | 1.0935 | 1.0950 | F122 | 1.1016274568 | 1.1025323513 | **1.1025323513** | +0.6879% | SDP |
| 18.0 | 1.0893 | **1.0944** | 1.0931 | 1.0944 | F122 | 1.1010754832 | 1.1019911074 | **1.1019911074** | +0.6936% | SDP |
| 18.5 | 1.0887 | **1.0938** | 1.0928 | 1.0938 | F122 | 1.1005552019 | 1.1014718646 | **1.1014718646** | +0.7014% | SDP |
| 19.0 | 1.0881 | **1.0933** | 1.0925 | 1.0933 | F122 | 1.1000614471 | 1.1010085878 | **1.1010085878** | +0.7051% | SDP |
| 19.5 | 1.0875 | **1.0928** | 1.0922 | 1.0928 | F122 | 1.0995912820 | 1.1005804499 | **1.1005804499** | +0.7120% | SDP |
| 20.0 | 1.0870 | **1.0923** | 1.0919 | 1.0923 | F122 | 1.0991451766 | 1.1000835957 | **1.1000835957** | +0.7126% | SDP |
| 20.5 | 1.0865 | **1.0918** | 1.0917 | 1.0918 | F122 | 1.0987203726 | 1.0995818248 | **1.0995818248** | +0.7128% | SDP |
| 21.0 | 1.0860 | **1.0914** | **1.0914** | 1.0914 | F122+PW | 1.0983166448 | 1.0993840223 | **1.0993840223** | +0.7315% | SDP |
| 21.5 | 1.0856 | 1.0909 | **1.0912** | 1.0912 | PW | 1.0979323423 | 1.0990236828 | **1.0990236828** | +0.7170% | SDP |
| 22.0 | 1.0852 | 1.0905 | **1.0909** | 1.0909 | PW | 1.0975661350 | 1.0986651015 | **1.0986651015** | +0.7118% | SDP |
| 22.5 | 1.0848 | 1.0901 | **1.0907** | 1.0907 | PW | 1.0972146966 | 1.0982427870 | **1.0982427870** | +0.6916% | SDP |
| 23.0 | 1.0845 | 1.0897 | **1.0905** | 1.0905 | PW | 1.0968809003 | 1.0979658741 | **1.0979658741** | +0.6846% | SDP |
| 23.5 | 1.0841 | 1.0893 | **1.0903** | 1.0903 | PW | 1.0965619548 | 1.0976871718 | **1.0976871718** | +0.6775% | SDP |
| 24.0 | 1.0838 | 1.0890 | **1.0901** | 1.0901 | PW | 1.0962555580 | 1.0973558414 | **1.0973558414** | +0.6656% | SDP |
| 24.5 | 1.0835 | 1.0886 | **1.0900** | 1.0900 | PW | 1.0959631055 | 1.0972829464 | **1.0972829464** | +0.6689% | SDP |
| 25.0 | 1.0832 | 1.0883 | **1.0898** | 1.0898 | PW | 1.0956811707 | 1.0970508460 | **1.0970508460** | +0.6653% | SDP |
| 25.5 | 1.0830 | 1.0880 | **1.0896** | 1.0896 | PW | 1.0954114575 | 1.0967874328 | **1.0967874328** | +0.6596% | SDP |
| 26.0 | 1.0827 | 1.0876 | **1.0895** | 1.0895 | PW | 1.0951512934 | 1.0963380206 | **1.0963380206** | +0.6276% | SDP |
| 26.5 | 1.0825 | 1.0873 | **1.0893** | 1.0893 | PW | 1.0948998987 | 1.0958948980 | **1.0958948980** | +0.6054% | SDP |
| 27.0 | 1.0823 | 1.0871 | **1.0892** | 1.0892 | PW | 1.0946585538 | 1.0956153198 | **1.0956153198** | +0.5890% | SDP |
| 27.5 | 1.0820 | 1.0868 | **1.0890** | 1.0890 | PW | 1.0944270345 | 1.0953656944 | **1.0953656944** | +0.5845% | SDP |
| 28.0 | 1.0818 | 1.0865 | **1.0889** | 1.0889 | PW | 1.0942039168 | 1.0951787825 | **1.0951787825** | +0.5766% | SDP |
| 28.5 | 1.0816 | 1.0863 | **1.0888** | 1.0888 | PW | 1.0939880327 | 1.0949477092 | **1.0949477092** | +0.5646% | SDP |
| 29.0 | 1.0814 | 1.0860 | **1.0886** | 1.0886 | PW | 1.0937793808 | 1.0947004759 | **1.0947004759** | +0.5604% | SDP |
| 29.5 | 1.0812 | 1.0858 | **1.0885** | 1.0885 | PW | 1.0935770275 | 1.0944995281 | **1.0944995281** | +0.5512% | SDP |
| 30.0 | 1.0810 | 1.0856 | **1.0884** | 1.0884 | PW | 1.0933813304 | 1.0943802573 | **1.0943802573** | +0.5495% | SDP |
| 30.5 | 1.0809 | 1.0854 | **1.0883** | 1.0883 | PW | 1.0931914580 | 1.0941923325 | **1.0941923325** | +0.5414% | SDP |
| 31.0 | 1.0807 | 1.0852 | **1.0882** | 1.0882 | PW | 1.0930081821 | 1.0941720124 | **1.0941720124** | +0.5488% | SDP |
| 31.5 | 1.0805 | 1.0850 | **1.0881** | 1.0881 | PW | 1.0928307671 | 1.0940430535 | **1.0940430535** | +0.5462% | SDP |
| 32.0 | 1.0804 | 1.0848 | **1.0880** | 1.0880 | PW | 1.0926588314 | 1.0938734896 | **1.0938734896** | +0.5398% | SDP |
| 32.5 | 1.0802 | 1.0847 | **1.0879** | 1.0879 | PW | 1.0924934508 | 1.0937109644 | **1.0937109644** | +0.5341% | SDP |
| 33.0 | 1.0800 | 1.0845 | **1.0878** | 1.0878 | PW | 1.0923331194 | 1.0936208702 | **1.0936208702** | +0.5351% | SDP |
| 33.5 | 1.0799 | 1.0844 | **1.0877** | 1.0877 | PW | 1.0921774255 | 1.0934700803 | **1.0934700803** | +0.5305% | SDP |
| 34.0 | 1.0797 | 1.0842 | **1.0876** | 1.0876 | PW | 1.0920273639 | 1.0931156543 | **1.0931156543** | +0.5071% | SDP |
| 34.5 | 1.0796 | 1.0841 | **1.0875** | 1.0875 | PW | 1.0918813665 | 1.0929515526 | **1.0929515526** | +0.5013% | SDP |

**The score.**

```
rows in the table            : 68
rows certified by BOTH       : 68
rows beating the record      : 68
rows NOT beating it          : 0
margin, smallest             : +0.2726%   (A = 1.5)
margin, largest              : +1.9733%   (A = 1.0)
   largest excluding A = 1   : +0.7315%   (A = 21.0)
margin, mean over all rows   : +0.6427%
LP won the row               : 16 rows
SDP won the row              : 52 rows
SDP minus LP, over the rows  : -7.28e-04 to +1.91e-02
```

**The margins survive the rounding.**  CQH print four decimals, so the
record on a row could really be as much as `5e-5` above the printed value.
Measuring every margin against `printed + 5e-5` instead — the worst case
the rounding allows — still leaves **all 68 rows ahead**, with margins from
`+0.2691%` to `+1.9707%` and a mean of `+0.6368%`.  So no row's win is an
artefact of reading a rounded table.  The one row where the unrounded
record is known independently is `A = 28`, where CQH's rigorous
`1/0.91833 = 1.0889331721712239` gives `+0.5736%` against the `+0.5766%`
computed from the printed `1.0889`.

**The two methods cross once, at `A = 9.5`.**  Below that the band-limited
linear program wins every row; from there up the semidefinite program
wins every row.  The single exception is `A = 1` itself, where the SDP
wins by a wide margin, for the reason given below.  A crossing of this
kind is what Chirre-Quesada-Herrera describe in their Section 7 and state
as their Conjecture 11.  In their own table it sits at `A = 21`; here it
sits at `A = 9.5`.  Both of this project's pipelines are stronger than the
published version of the same method, and they are not stronger by the
same amount, so the crossing moves.

**The headline row, `A = 28`.**

```
LP  certified   1.0942039167907789
SDP certified   1.0951787824934853      <- the row winner
CQH's rigorous record 1/0.91833 = 1.0889331721712239
certified upper bound (UPPER.md) 1.1079
margin over the record          +0.5766%
```

**The row that used to fail, `A = 1`.**  The previous round recorded
`A = 1` as the one row of 68 that lost to the published record, by 3.08%.
Both pipelines now clear it.  The reason the linear program used to lose
is that at `A = 1` the penalty is cheap, so the best function wants a very
wide band, and a wide band spread over a capped number of harmonics
cannot resolve it.  Widening the band and paying for the grid fixes it:

```
B = 6,  K = 500   certified 1.9284575941
B = 8,  K = 560   certified 1.9572273260
B = 12, K = 600   certified 1.9797820564      <- used
```

That is still climbing with `B`, so the `A = 1` LP number is not
converged; it is simply the best one paid for.  The Gaussian family has no
band at all, so the obstruction does not exist there, and the SDP
certifies `1.9988807616`.  CMS Theorem 2(b) proves `C_+(1) = 2`
exactly, so this row was never a standing record; the certified value is
99.944% of the exact answer.

### Sanity checks

**1. The three regressions.**  Three saved files must certify to the same
digits as before.  Re-run in fresh processes during this consolidation,
they do:

```
cplus_certA.py    testfunction_A28.txt                 -> 1.0940451158772697
cplus_sdp_cert.py sdpfun_A28.0_d42_c1.2.txt   28 1.2   -> 1.0951787824934853
cplus_sdp_cert.py sdpfun_A1.0_d50_c0.45.txt   1  0.45  -> 1.9903625363954200
```

All three match bit for bit.  Note that the first is the older, smaller
`A = 28` test function of the section above, not a table row; the table's
`A = 28` LP entry is a different and better function.

**2. The tripwire.**  `UPPER.md` proves certified *upper* bounds on
`C_+(A)` at `A = 2, 5, 15, 28, 34.5`.  `C_+` never increases with `A`,
because the penalty term grows with `A`, so an upper bound proved at some
`A0` is also an upper bound at every `A >= A0`.  Every certified lower
bound in the table was checked against the bound at the nearest such `A0`
below it.  **No row violates it.**  A value above that line would have
been a bug, not a record.

**3. Certified below the ordinary float64 computation.**  A rigorous
lower bound must come out below a plain float64 evaluation of the same
function.  It does on every LP row.  On a handful of SDP rows the
certificate sits a few parts in `1e12` *above* the float64 number; that
is the float64 number being wrong, not the certificate.  At `A = 1` the
gap is larger, `6.9e-5`, and it was checked against an independent
40-digit mpmath evaluation of the same file
(`cplus_sdp_mp.py sdpfin_A1.0_d50_c0.15.txt 1 0.15`), which returns
`1.9988807616146563797` against the certificate's `1.9988807616146564`.
The float64 recomputation uses a fixed scan grid and is known to be
unreliable at small `A` with a small dilation; the certificate is not.

**4. Nine rows re-certified from scratch.**  The JSONL files are not
trusted on their own.  Nine rows spanning the whole range were re-run in
**fresh processes**, reading only the saved coefficient file.  Every one
reproduces the recorded value to the last printed digit:

| row | file | recorded | re-certified | match |
|---|---|---|---|---|
| A = 1.0 SDP | `sdpfin_A1.0_d50_c0.15.txt` | 1.9988807616146564 | 1.9988807616146564 | yes |
| A = 1.5 LP | `testfunction_A1.5_r28.txt` | 1.4127410698406371 | 1.4127410698406371 | yes |
| A = 4.0 LP | `testfunction_A4.0_r28.txt` | 1.1750163119424673 | 1.1750163119424673 | yes |
| A = 9.0 LP | `testfunction_A9.0_r28.txt` | 1.1208381955374593 | 1.1208381955374593 | yes |
| A = 12.0 SDP | `sdpfin_A12.0_d42_c1.2.txt` | 1.1114469381465950 | 1.1114469381465950 | yes |
| A = 21.0 SDP | `sdpfin_A21.0_d42_c1.1.txt` | 1.0993840222943752 | 1.0993840222943752 | yes |
| A = 28.0 LP | `testfunction_A28.0_r27raw.txt` | 1.0942039167907789 | 1.0942039167907789 | yes |
| A = 28.0 SDP | `sdpfun_A28.0_d42_c1.2.txt` | 1.0951787824934853 | 1.0951787824934853 | yes |
| A = 34.5 SDP | `sdpfin_A34.5_d42_c1.15.txt` | 1.0929515526238738 | 1.0929515526238738 | yes |

**Nothing failed to reproduce.**  An earlier draft of this section claimed
the `A = 28` LP value `1.0942039167907789` did not reproduce.  That was a
file mix-up, not a discrepancy: there are two `A = 28` LP files at the same
`B = 1.7, K = 190`, and they certify two different numbers, each
deterministically.  `testfunction_A28.0_r28.txt` gives
`1.0942039156939909` and `testfunction_A28.0_r27raw.txt` gives
`1.0942039167907789`.  The table uses the larger, and it is the one
re-certified above.  The two differ by `1.1e-9`, which is `1.7e-7` of the
margin over the record, so nothing depends on the choice.  Both are
rigorous lower bounds, so either is safe to quote.

**5. Every row is independently checkable.**  The consolidation script
drops any certified value whose coefficient file is no longer on disk, so
no number in the table rests on a file that cannot be re-certified.  Eleven
LP values were dropped by that rule (the `testfunction_A<A>w.txt` and
`w2.txt` variants of an earlier round).  **None of them would have won its
row**: each is below the value actually kept, so the rule removed nothing
from the table.  In particular the `A = 1` LP value `1.8998323` is not in
the table because a later, wider-band run certifies `1.9797820564`, which
is larger — not because its file is missing.  A copy of that value with a
surviving coefficient file (`testfunction_A1.0.txt`) is present and is
superseded on value alone.

**6. No table entry rests on an overwritten file.**  Taking the maximum
over every JSONL ever written is only sound if a filename still holds the
coefficients it held when its number was produced.  Some do not: six
filenames carry two different certified values in two different JSONLs,
because the file was rewritten between rounds —

```
testfunction_A28.0_r28.txt   1.0942038986  and  1.0942039157
testfunction_A1.0.txt        1.7653367680  and  1.8998323076
testfunction_A1.5.txt        1.4082760586  and  1.4116512919
testfunction_A2.0.txt        1.2972676321  and  1.2982302204
testfunction_A3.0.txt        1.2109673955  and  1.2116201425
testfunction_A4.5.txt        1.1632592060  and  1.1633719867
```

Every one of those is a generic, round-keyed name.  **None of them backs
any entry in the table** — not a winning column and not a losing one.
Every entry above comes from a name that encodes what makes it unique
(`_r28`, `_r27raw`, `_r28B12`, `sdpfin_A<A>_d<d>_c<c>`,
`sdpfun_A<A>_d<d>_c<c>`), and no such name carries more than one value.
So the ambiguity is real but it is confined to superseded rows.

### Caveats

* **The margins are small.**  Except at `A = 1`, every improvement over
  the published record is between about `+0.27%` and `+0.75%`.
* **These are lower bounds only.**  Nothing here says what `C_+(A)` is.
  The certified brackets in `UPPER.md` are still 1.8%-2.7% wide.
* **Neither pipeline was pushed to its limit.**  The LP is noise-limited,
  not degree-limited: its coefficients plateau near `1.5e-8` past about
  the 48th harmonic, so raising `K` far beyond `190` buys nothing.  The
  SDP is solver-limited: CLARABEL reports `optimal_inaccurate` above
  block degree `d ~ 45`, so the runs here stop at `d = 42` (`d = 50` at
  `A <= 2`).  Both ceilings are tooling, not mathematics.
* **The search settings were interpolated, not scanned, on most rows.**
  The LP bandwidth `B` and grid density, and the SDP dilation `c`, were
  measured at a few `A` and interpolated in between.  Rows away from those
  anchors are probably not the best either method can do.  Measured cost
  of this where it was checked: a few parts in `1e6` for the LP.
* **The SDP candidate ranking is float64.**  Six candidates per row are
  ranked by a float64 recomputation and only the top two are certified.
  As sanity check 3 shows, that ranking is not always right.  A bad
  ranking can only cost a better answer; it cannot produce a wrong one,
  because the number reported is always the certificate.
* **`A = 1` is not a record.**  CMS prove `C_+(1) = 2`.
* **What "certified" rests on.**  ARB ball arithmetic - `Si`, `sin_pi`,
  `cos_pi`, `sinc_pi` for the LP; `erf`, `exp`, `sqrt` and the incomplete
  gamma for the SDP - and a human reading of the problem statement out of
  the CMS and CQH PDFs.  Both assumptions are the ones the published work
  makes too.

### Reproduce

```
# the two searches, both resumable, both writing one JSONL row per A
python -u cplus_final_lp.py      # -> final_lp*.jsonl, testfunction_A<A>_r28.txt
python -u cplus_final_sdp.py     # -> final_sdp*.jsonl, sdpfin_A<A>_d<d>_c<c>.txt
python -u cplus_final_sdp2.py    # walks the SDP dilation c off the grid edge
python -u lp_A1.py               # the A = 1 LP row, which needs a wide band

# consolidate and write this section
python cplus_final_table.py          # -> final_table.json, the table + the score
python cplus_final_report.py --write # splices the section into EXTREMAL2.md

# independent re-merge of every certified value on disk, from the shards up
python consolidate.py                # -> consolidated.jsonl, the table + the score

# re-certify any single row from its file alone
python cplus_certA.py testfunction_A28.0_r28.txt
python cplus_sdp_cert.py sdpfin_A28.0_d42_c1.2.txt 28 1.2
```

Set `SHARD=i/n` to split either search across `n` processes.  Each shard
needs its own output file (`LPOUT`, `SDPOUT`, `SDP2OUT`).

## What this is NOT

- **Not verified to be the current record.** That CQH's 1.0889 still stands comes
  from searches plus a November 2024 expository article (Quesada-Herrera,
  arXiv:2411.05095) that still quotes it, not from a systematic check of everything
  published since.
- **Not the sharp constant.** This is a lower bound. The certified bracket for
  `C_+(28)` is now `[1.0951787825, 1.1079]`, the upper end from `UPPER.md`. The
  flatness in `B` suggests the bandlimited family is close to exhausted, but that
  is an observation, not an argument.
- **Not a new method.** Bandlimited test functions for this problem are CQH's own
  idea, and they say in Section 7 that they expect bandlimited functions to beat
  SDP-with-Gaussians as `A` grows. This is that expectation carried further with a
  bigger solver, not a new approach.
- **Not large.** A 0.58% improvement to one constant in a conditional bound on gaps
  between primes represented by a quadratic form.

## Reproduce

```
python cplus_calib.py                       # Fejer, CMS c_0, and the CQH function
python cplus_scan.py 28                     # coarse scan over B and parity
python cplus_refine.py 28 28                # fine scan over B at K = 140
python cplus_final.py                       # the one LP quoted; writes cplus_final_28.npz
python cplus_verify.py cplus_final_28.npz   # admissibility and grid-independence
python cplus_mp.py cplus_final_28.npz 120   # definitive 50-digit evaluation
python cplus_dump.py cplus_final_28.npz testfunction_A28.txt
python cplus_table.py                       # the same LP across the published table
python cplus_cert.py testfunction_A28.txt   # the ARB certificate (3.5 s)
python cplus_certtable.py                   # the same, for all 68 rows (~40 min)
python cplus_certreport.py                  # the certified table and its checks
```

`cplus_cert.py` takes two environment variables: `UEXACT` (the cutoff `u_c`,
default 400) and `NOPROJ=1` (skip the projection, diagnostic only).

## Novelty

**New here:** the numbers. `C_+(28) >= 1.0951787825`, certified in ARB, against
the published certified 1.0889332, and — see *The full table, certified* — **two**
certified lower bounds for every one of the 68 rows of CQH's Table 1, one from a
bandlimited linear program and one from a semidefinite program, beating the best
of their three columns on all 68 rows by 0.27% to 1.97%, with the coefficient file
for each of the 136 test functions saved so it can be rechecked on its own. Also
the observation of *why* the 3-parameter family stopped where it did — its value
collapses when the band widens, so the published search was locked to
`B ≈ 1.014`.

**Not new:** the extremal problem (CMS 2017), the application to quadratic forms and
the constant 28 (CQH 2022), the idea of using bandlimited test functions here and the
prediction that they beat Gaussian-times-polynomial SDP at large `A` (CQH Section 7,
including their Conjecture 11), and the linear-programming discretisation itself,
which is standard for Beurling–Selberg problems and was already used in `EXTREMAL.md`.
The only thing done differently is solving over 71 coefficients and the bandwidth at
once instead of over three integers and a dilation.

Also **not new: the certification method.** Ball arithmetic in ARB for exactly
this kind of extremal-problem bound is what CQH did. Nothing in the Certification
section is a new technique; the closed-form `Si` antiderivative, the
cancellation-free `W` form and the `Fhat''<0` argument at `t=B` are
conveniences specific to this basis, not new mathematics.
