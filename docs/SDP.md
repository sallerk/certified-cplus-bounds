# The published SDP column was scored with a loose relaxation; scored exactly, the same function class beats the linear program on 6 of 9 rows — and fixes the one row this project was losing

> **Note on the numbers in this file.** This documents the semidefinite-programming
> round on its own. Several values were improved in the final consolidation, where
> both solvers were run on every row and the better kept. In particular `A = 1` reads
> `1.99036` here but the package value is `1.9988807616` (99.944% of the exact
> `C_+(1) = 2`, CMS Theorem 2(b)). **The authoritative values are `final_table.json`
> and the table in `EXTREMAL2.md`**; this file is kept for the method and its
> calibration, not for its headline figures.


**Novelty, stated first.** The semidefinite program is **not new and not mine**. It
is Chirre–Pereira Júnior–de Laat, *Primes in arithmetic progressions and
semidefinite programming*, Math. Comp. **90** (2021) 2235–2246
(arXiv:2005.02393), Lemma 4 and Section 4. I searched first, found their
formulation in print, and used it rather than inventing one. The `F82` and
`F122` columns of Chirre–Quesada-Herrera's Table 1 were computed with exactly
this machinery by David de Laat. Certifying such a bound in ARB ball arithmetic
is also theirs.

What is new here is three things:

1. **Ten certified numbers.** A rigorous lower bound on `C_+(A)` at ten values of
   `A`, each above the best published value, and six of them above this
   project's own certified linear-programming value.
2. **The row this project loses is now won.** At `A = 1` the linear program
   certifies `1.89983` against a published `1.96020` — a 3.08% loss recorded in
   `EXTREMAL2.md`. The SDP certifies **`1.99036`**, which is 99.52% of the exact
   answer `C_+(1) = 2` and beats the published record by 1.54%.
3. **A diagnosis of why the published SDP column looks weak.** The number the
   SDP maximises is not the value of the function it returns. It is the value of
   a *relaxation* of that value. At `A = 28` and degree 170 the relaxation reads
   `1.08810` while the same function is honestly worth `1.09518` — a gap of
   `7.1e-3`, or 0.65%. So the published comparison "bandlimited beats
   Gaussian-times-polynomial at large `A`" is a comparison between a relaxed
   score and an exact one.

None of this is large. It is a fraction of a percent on a lower bound for a
constant inside a conditional bound on prime gaps.

---

## 1. What the problem is

Let `F` be an even, continuous, real function on the line whose absolute value
has a finite integral. Its **Fourier transform** is

```
Fhat(t) = int F(x) e^{-2 pi i x t} dx
```

Write `(y)_+` for `max(y, 0)`, and `||F||_1` for `int |F(x)| dx`. For a number
`A > 1`,

```
C_+(A) = sup over all such F of  [ F(0) - A int_{|t|>1} (Fhat(t))_+ dt ] / ||F||_1
```

`sup` means "the largest value reachable". Because it is a largest value, any
single `F` you can write down gives a **lower bound**. That is all this file
produces.

The source is Carneiro–Milinovich–Soundararajan, arXiv:1708.04122, Extremal
Problem 2, equation (1.3). `A = 28` is the value Chirre–Quesada-Herrera
(arXiv:2012.07781) need for gaps between primes represented by a binary
quadratic form. Their Table 1 lists lower bounds at 68 values of `A`.

## 2. What the literature actually does (searched first)

Two papers matter, and both were read from the PDF, not from a summary.

**Chirre–Pereira Júnior–de Laat, Math. Comp. 90 (2021) 2235–2246.** Their
Lemma 4 rewrites the problem so that a computer can attack it. It replaces the
single function `F` by four **non-negative** functions:

```
C_+(A) = sup { f1(0) - f2(0) - A int_{|t|>1} f3(t) dt }
```

over tuples `(f1, f2, f3, f4)` of even, non-negative, continuous, integrable
functions with

```
int f1 + int f2 = 1        (a normalisation)
f1 - f2 = f3hat - f4hat    (an identity between functions)
```

The idea: `f1` is the positive part of `F`, `f2` the negative part, `f3` the
positive part of `Fhat`, `f4` the negative part. Then `||F||_1 = int f1 + int f2`
and `(Fhat)_+ = f3`, so the two problems agree.

Why bother? Because "non-negative" is a condition a semidefinite program can
enforce *exactly*, whereas "`|F|`" is not linear. Their Section 4 then restricts
every `f_i` to

```
f_i(x) = p_i(x^2) e^{-pi x^2},   p_i a polynomial, p_i >= 0 on [0, infinity)
```

and writes the non-negativity as a **sum-of-squares certificate**:

```
p(u) = v(u)^T Q v(u) + u w(u)^T R w(u),     Q, R positive semidefinite
```

This is the Markov–Łukasiewicz theorem: a polynomial is non-negative on
`[0, infinity)` **if and only if** it has this form. It is an exact finite
description of an infinite family of inequalities. That is the whole reason to
use an SDP rather than a linear program: a linear program can only impose
`p(u_j) >= 0` at sample points `u_j`, which is a relaxation, and this is not.

They report `C_+(36/11) >= 1.1965` with matrices of size 90, and they certify
their numbers in ARB.

**Chirre–Quesada-Herrera, arXiv:2012.07781, Section 7.** Their Table 1 has three
columns. `F82` and `F122` are the CPdL construction with `deg P <= 82` and
`deg P <= 122` where `F(x) = P(x) e^{-pi x^2}`; both were computed for them by
de Laat. `PW` is their own hand-built bandlimited function. Their Proposition 10
and Conjecture 11 predict the polynomial-times-Gaussian family gets worse as `A`
grows.

**So the formulation used below is theirs, line for line.** Section 3 records the
one place I had to work something out for myself, and Section 6 records the two
places where I extended it.

## 3. The two bases, and why the choice matters

Everything lives in the space of "even polynomial times Gaussian" functions. The
Gaussian is `e^{-pi x^2}`, which is the one that is its own Fourier transform
under the convention above. That is forced: the family has to be closed under
the Fourier transform for the identity `f1 - f2 = f3hat - f4hat` to make sense
inside it.

Two different bases are needed, and using the wrong one wrecks the arithmetic.

**Basis for the functions: the even Hermite functions.**

```
h_k(x) = eta_k L_k^{-1/2}(2 pi x^2) e^{-pi x^2},    eta_k = sqrt( sqrt(2 pi) k! / Gamma(k+1/2) )
```

`L_k^{-1/2}` is the generalised Laguerre polynomial with parameter `-1/2`. These
satisfy two things at once:

```
int h_k h_m dx = 1 if k = m, else 0        (orthonormal)
h_k-hat = (-1)^k h_k                       (Fourier eigenfunctions)
```

The second line is the key. **In this basis the Fourier transform is a sign
flip.** Nothing is computed; the transform costs nothing and introduces no
error. The constraint `f1 - f2 = f3hat - f4hat` becomes, on coefficient vectors,

```
c1 - c2 = ((-1)^k) * (c3 - c4)     componentwise
```

**Basis for the Gram matrices: orthonormal for the Gaussian weight.**

```
v_a(x) = L_a^{-1/2}(pi x^2) / sqrt( Gamma(a+1/2) / (a! sqrt(pi)) )
```

with `int v_a v_b e^{-pi x^2} dx = 1 if a = b, else 0`. The Laguerre *argument*
is `pi x^2` here and `2 pi x^2` above; that is not a typo, it is the whole
point. If the same argument is used for both, the matrix of moments
`int v_a v_b e^{-pi x^2} dx` has entries of size about `16^d`, and at `d = 30`
that is `1e36` — the SDP is unsolvable in double precision. With the `pi x^2`
argument that matrix is the identity.

De Laat says only "for the numerical conditioning a good choice for the basis is
Laguerre". This section is what that sentence means in practice.

**The map from Gram matrices to Hermite coefficients** is then

```
c[k] = <Q, G_k> + <R, H_k>,
G_k[a,b] = int v_a v_b e^{-pi x^2} h_k dx,   H_k[a,b] = int (pi x^2) v_a v_b e^{-pi x^2} h_k dx
```

computed once by Gauss–Hermite quadrature, which is exact here because every
integrand is a polynomial times `e^{-2 pi x^2}`. Every entry is of size about 1.

**Sizes.** With `Q` and `R` both `(d+1) x (d+1)`, the polynomial `p` has degree
`2d+1` in `x^2`, so `F(x) = p(x^2) e^{-pi x^2}` has degree `4d + 2` in `x`.
`d = 20` gives degree 82 and `d = 30` gives degree 122 — exactly CQH's two
columns.

## 4. The tooling

Already installed, nothing added: **cvxpy 1.9.2**, **CLARABEL 0.11.1**,
**SCS 3.2.11**, **python-flint 0.9.0** (which is ARB), numpy, scipy, mpmath.
Not installed and not needed: mosek, cvxopt, picos, sdpa-gmp.

**CLARABEL is the solver.** SCS was tried and is not competitive: at `A = 28`,
`d = 20`, asked for `eps = 1e-11` it ran to its iteration cap and returned
`1.0204` where CLARABEL returns `1.0813` in 0.8 seconds. All numbers below are
CLARABEL. Above about `d = 45` CLARABEL starts reporting `optimal_inaccurate`,
and at `d = 60` on one row it raised a `SolverError`. That is the ceiling of
this setup, and it is a double-precision ceiling: CPdL used `sdpa-gmp`, which is
multiprecision.

The problem sizes are small: eight positive-semidefinite blocks of side `d+1`,
plus `2d + 3` linear equations. At `d = 42` (degree 170) one solve takes about
12 seconds.

## 5. Calibration — five checks, done before any claim

**5.1 The basis identities.** `sdp_check.py` verifies five things numerically at
`d = 8, 20, 30`: that the `h_k` are orthonormal; that `h_k-hat = (-1)^k h_k`, by
integrating `2 int_0^inf h_k(x) cos(2 pi x t) dx` and comparing; that the
Gram-matrix-to-coefficient map reproduces `p(x) e^{-pi x^2}` pointwise; that
`F(0)` and `int F` come out right; and that the "tail" vector
`int_{|t|>c} h_k dt` is right. Worst error over all of them, at `d = 30`:
`2.9e-13` relative.

**5.2 The Gaussian, in closed form.** `F(x) = e^{-pi x^2}` has `F(0) = 1`,
`||F||_1 = 1`, `Fhat = F`, so its ratio at `A = 28` is exactly
`1 - 28 erfc(sqrt(pi))`. Pushed through the **certificate**:

```
certified   0.6587112988255192
closed form 0.658711298825519167   (mpmath, 30 digits)
```

Sixteen digits, and the certificate correctly reports zero sign changes and
`||F||_1 = 1`. This tests the whole chain — evaluation, sign sweep, closed-form
integration, tail bound — against a number that is known exactly.

*(The Fejér calibration used in `EXTREMAL.md` cannot be run here: the Fejér
kernel has a compactly supported transform and is not in this family. The
Gaussian plays the same role.)*

**5.3 Against the published SDP numbers — close, but NOT a reproduction.** At
`A = 28`, with the SDP's own objective (which is what the papers report):

| what | degree | this code | published | difference |
|---|---|---|---|---|
| `F82` | 82 (`d = 20`) | 1.0812851 | 1.0818 | `-5.1e-4` |
| `F122` | 122 (`d = 30`) | 1.0860884 | 1.0865 | `-4.1e-4` |

and at `A = 36/11`, where CPdL report `1.1965` with matrices of size 90:

| `d` | 20 | 30 | 45 | 90 |
|---|---|---|---|---|
| this code | 1.1894934 | 1.1920176 | 1.1950901 | (CLARABEL fails) |
| CPdL | | | | 1.1965 |

**I did not reproduce a published SDP number exactly.** I land about `4e-4`
low at matched degree, consistently, in the same direction, and the trend in
`d` is consistent with their value. Two candidate explanations, neither
settled: their solver is multiprecision and mine is not (Section 4 shows the
double-precision ceiling is real); or their nominal degree corresponds to a
slightly different split of the sum-of-squares certificate. The second was
tested and does not fit: shifting the degree by one step gives `d = 21`
`-> 1.0816435`, still *below* their 1.0818, while `d = 31 -> 1.0866847` is
*above* their 1.0865. No single shift matches both. The gap is unexplained and
is stated as such.

**5.4 The value is monotone in the degree.** The SDP objective must not decrease
as the matrices grow, since the feasible set only grows. `A = 28`, `c = 1`:

| `d` | 6 | 10 | 20 | 30 | 40 | 48 | 54 | 60 |
|---|---|---|---|---|---|---|---|---|
| degree | 26 | 42 | 82 | 122 | 162 | 194 | 218 | 242 |
| SDP objective | 1.0565356 | 1.0668467 | 1.0812851 | 1.0860884 | 1.0887839 | 1.0900300 | 1.0917387 | 1.0944168 |

Strictly increasing across the whole range. Nothing enforces that.

**5.5 Against the proven upper bound.** `UPPER.md` certifies
`C_+(28) <= 1.1136`. Every number below is under it. The headline `1.0951788`
sits `4.8%` of the way across the remaining bracket `[1.0942036, 1.1136]`.

## 6. Two changes to the published setup

**6.1 The scale of the Gaussian is a free parameter, and nobody used it.**
`P(x) e^{-pi x^2}` fixes a width. But the extremal problem is not scale
invariant — the interval `[-1, 1]` is fixed — so `F(x/c)` is a genuinely
different test function for each `c > 0`. Writing it out:

```
ratio for F(x/c)  =  [ F(0) - A int_{|s|>c} (Fhat(s))_+ ds ] / ( c ||F||_1 )
```

So the SDP is solved with the penalty taken over `|t| > c` instead of `|t| > 1`,
and the answer divided by `c`. Everything else — the sign conditions, the
Fourier eigenbasis, the certificate — is untouched; only one moment vector
changes. `c` is then scanned.

It matters, and how much depends on `A`. Degree held fixed within each column,
so the only thing that changes is `c`:

| A | 1.0 | 4.0 | 28.0 |
|---|---|---|---|
| degree | 162 (`d = 40`) | 122 (`d = 30`) | 170 (`d = 42`) |
| value at `c = 1` | 1.9463466 | 1.1740413 | 1.0948646 |
| best `c` found | 0.45 | 0.60 | 1.20 |
| value at best `c` | **1.9880708** | **1.1746286** | **1.0951788** |

At `A = 28` it is worth `+3.1e-4`. At `A = 4` it is worth `+5.9e-4`. At `A = 1`
it is worth `+0.0417`, which is most of the margin over the published record.
(The best row in Section 8 at `A = 1` is degree 202, not 162; at `c = 0.45` that
gives 1.9903685.)

**6.2 The SDP objective is not the value of the function it returns.** This is
the finding of the round, and it is a property of CPdL's Lemma 4, not a bug.

The Lemma is an equality *over all* non-negative `f_i`. Restricted to polynomials
of bounded degree, the two inequalities inside its proof stop being equalities:

```
||F||_1  <=  int f1 + int f2        (equal only if f1 = F^+ and f2 = F^- exactly)
(Fhat)_+ <=  f3                     (equal only if f3 = (Fhat)_+ exactly)
```

`F^+` — the positive part of a polynomial times a Gaussian — is **not** a
polynomial times a Gaussian. It has corners where `F` crosses zero. So no finite
degree can represent it, and the two inequalities are strict. The SDP therefore
under-reports the function it hands back.

Measured. For each row, the SDP objective against the **certified** value of the
very same `F = f1 - f2`:

| A | degree | SDP objective | same `F`, certified | gap |
|---|---|---|---|---|
| 1.0 | 202 | 1.9903628 | 1.9903625 | `-2.4e-07` |
| 1.5 | 146 | 1.4109739 | 1.4124486 | `+1.47e-03` |
| 36/11 | 82 | 1.1894934 | 1.1988170 | `+9.32e-03` |
| 4.0 | 122 | 1.1705740 | 1.1746286 | `+4.05e-03` |
| 6.0 | 146 | 1.1376199 | 1.1412931 | `+3.67e-03` |
| 10.0 | 170 | 1.1123587 | 1.1170611 | `+4.70e-03` |
| 14.0 | 170 | 1.1020093 | 1.1073745 | `+5.37e-03` |
| 21.0 | 170 | 1.0931291 | 1.0993840 | `+6.25e-03` |
| 28.0 | 170 | 1.0880995 | 1.0951788 | `+7.08e-03` |
| 34.5 | 170 | 1.0869257 | 1.0928295 | `+5.90e-03` |

The gap is essentially zero at `A = 1` and grows to `7e-3` at `A = 28`. That is
the shape you would expect: when `A` is large the optimal `F` oscillates more,
so `F^+` has more corners, so the polynomial approximation to it is worse.

**The consequence for the published comparison.** At `A = 28`, degree 82, the
SDP objective here is `1.0812851` — but the same function, certified, is
**`1.0904668`**. The bandlimited record at `A = 28` is `1.0889332`. So a
degree-82 Gaussian-times-polynomial function already beats the bandlimited
record at `A = 28`, and the published table says the opposite only because the
two columns are scored by different functionals.

*This is stated about **my** degree-82 function, not de Laat's — I do not have
his coefficients.* His was optimised for the relaxed objective, as mine was, so
the same argument should apply to it; but that is an inference, flagged as one.

## 7. The certificate

A lower bound is only worth what its proof is worth. `cplus_sdp_cert.py` turns
a coefficient file into a rigorous statement using ARB ball arithmetic, at 300
bits. A **ball** is a pair `[centre +/- radius]` guaranteed to contain the true
value; arithmetic on balls returns balls that still contain the true answer.

Three things make this class *easier* to certify than the bandlimited one.

**No admissibility condition.** A polynomial times a Gaussian is even, real,
continuous and integrable no matter what the coefficients are. The bandlimited
certificate `cplus_certA.py` has to force `Fhat(B) = 0` exactly (otherwise `F`
decays like `1/x` and `||F||_1` is infinite) and then project the coefficients
onto that plane. Here there is nothing to enforce and nothing to project. **The
file is the function.**

**The Fourier transform is exact.** Flip the sign of the odd-index coefficients.

**The integral has a closed form.** Writing `psi_n` for the Hermite function
normalised so `int psi_n^2 = 1`, and using
`psi_n = sqrt((n-1)/n) psi_{n-2} - sqrt(2/n) psi_{n-1}'`,

```
J_n(a,b) := int_a^b psi_n = sqrt((n-1)/n) J_{n-2}(a,b) - sqrt(2/n) [psi_{n-1}(b) - psi_{n-1}(a)]
J_0(a,b) = pi^{-1/4} sqrt(pi/2) [ erf(b/sqrt2) - erf(a/sqrt2) ]
```

So **there is no quadrature error anywhere**, exactly as the `Si` antiderivative
does for the bandlimited certificate. All that is left is deciding the *sign* of
the integrand, because `int |F| = |int F|` wherever the sign does not change.

**The one hard part: enclosing `F` on an interval.** Feeding an interval into
the Hermite three-term recurrence does not work. Ball arithmetic does not track
that the same `y` appears in every step, so the radius compounds. Measured at
`d = 20`, `x = 3`:

| input radius | 1e-3 | 1e-5 | 1e-8 | 1e-11 |
|---|---|---|---|---|
| output radius | 1.79e13 | 1.76e11 | 1.76e8 | 1.76e5 |

A constant amplification of `1.8e16`, against a true derivative of order 1. No
amount of precision fixes it and no width floor fixes it: it is an amplification
of the *input width*.

Worse, the two available a-priori bounds on `|psi_n|` are also useless in the
middle range. Both

```
|psi_n(y)| <= sqrt(2) (n + 1/2)^{1/4}                (global; proved below)
|psi_n(y)| <= r_n(y) e^{-y^2/2}, r_n the plus-sign recurrence   (decaying)
```

are far above `|F|` for `2.5 < x < 6`, because `F`'s decay there **is** a
cancellation between Hermite functions that are individually of size 1.

The fix is a **Taylor model**. `F` and its first `J-1` derivatives are evaluated
at the exact midpoint of each interval, where there is no amplification at all
(radius about `1e-90`). Only the `J`-th derivative is bounded crudely, and it is
multiplied by `(w/2)^J / J!`, which destroys the crudeness. Measured at `J = 6`
on the headline function:

| width | 2e-3 | 5e-4 |
|---|---|---|
| total enclosure radius (worst of `x = 0.5, 2.5, 4.5`) | 2.47e-3 | 6.18e-4 |
| of which the `J`-th remainder | 9.6e-14 | 2.3e-17 |

So the radius is just the **true** first-order variation `\|F'(c)\| w/2` — it
halves when the width halves — and the crude part contributes at the `1e-11`
level or below. That is the difference between a bisection that terminates and
one that does not.

The derivatives are exact: in the `psi` basis, differentiation is the
tridiagonal map `(Da)_m = a_{m+1} sqrt((m+1)/2) - a_{m-1} sqrt(m/2)`.

**The global bound is proved, not looked up.** `psi_n(y)^2 = 2 int_{-inf}^y psi psi'`,
so by Cauchy–Schwarz `psi_n(y)^2 <= 2 ||psi_n||_2 ||psi_n'||_2`. Here
`||psi_n||_2 = 1`, and
`||psi_n'||_2^2 = -int psi psi'' = (2n+1) - int y^2 psi^2 = (2n+1) - (n+1/2) = n+1/2`.
Hence `|psi_n| <= sqrt(2)(n+1/2)^{1/4}`. It is loose (4.26 against a true 0.8 at
`n = 82`) but it is elementary, self-contained and only ever multiplied by
something tiny.

**The tail.** For `x >= X`, `|F| <= R(y) e^{-y^2/2}` with `R` a polynomial with
non-negative coefficients, so `R(y) <= R(Y)(y/Y)^N` for `y >= Y`, and
`int_Y^inf y^N e^{-y^2/2} dy = 2^{(N-1)/2} Gamma((N+1)/2, Y^2/2)` is an
incomplete gamma function, which ARB gives as a ball. `X` is chosen so this is
below `1e-20`.

**The sweep.** `[0, X]` and `[c, X]` are swept adaptively. Runs of intervals with
the same certified sign are merged and integrated in closed form. Intervals
whose sign never settles are charged `width x (bound on |F| there)`, which
over-estimates and is therefore safe. **Nothing anywhere assumes that every sign
change was found**: a missed one shows up as an interval that never settles and
is charged at the over-estimate rate. In practice the undecided length is about
`1e-11` and the charge about `1e-25`.

**Cross-checks on the certificate.**

* The Gaussian, against a closed form: 16 digits (Section 5.2).
* An independent 40-digit mpmath re-evaluation of the headline function
  (`cplus_sdp_mp.py`: different library, different root finder, different
  quadrature) gives `1.0951788132` where the certificate gives `1.0951787825`.
  The certificate is below, as it must be, by `3.1e-8`.
* On every row the certified value comes out below the float64 recomputation, as
  it must. Certification costs between `3e-8` and `5e-8`.
* **Independent of the seed partition.** The headline row certified with 4000
  seed intervals and with 2500 gives `1.0951787824934853` both times, digit for
  digit. The adaptive part of the sweep is doing the work, not the starting grid.
* Cost: 1 to 7 seconds per row.

## 8. The result

`record` is the largest of CQH Table 1's three columns, and `from` says which
column that was; at `A = 28` their own rigorous form `1/0.91833` is used
instead of the rounded `1.0889`. `LP certified` is this project's best certified
bandlimited value (`FREEDOMS.md` for eight rows, `EXTREMAL2.md` for `A = 1`).
`deg` is the degree of the polynomial `P` in `F(x) = P(x) e^{-pi x^2}`, and `c`
is the dilation of Section 6.1.

| A | CQH record | from | LP certified | SDP certified | deg | c | SDP vs LP | SDP vs record |
|---|---|---|---|---|---|---|---|---|
| 1.0 | 1.960200 | PW | 1.8998323076 | **1.9903625364** | 202 | 0.45 | **+4.7652%** | **+1.5387%** |
| 1.5 | 1.408900 | F122 | 1.4127406445 | 1.4124486135 | 146 | 0.70 | −0.0207% | +0.2519% |
| 36/11 | — | — | (1.1986243, uncertified) | 1.1988170354 | 82 | 1.00 | — | +0.19% vs CPdL 1.1965 |
| 4.0 | 1.167300 | F122 | 1.1750094563 | 1.1746286326 | 122 | 0.60 | −0.0324% | +0.6278% |
| 6.0 | 1.133900 | F122 | 1.1415736559 | 1.1412930987 | 146 | 0.60 | −0.0246% | +0.6520% |
| 10.0 | 1.110100 | F122 | 1.1168052478 | **1.1170610631** | 170 | 1.00 | +0.0229% | +0.6271% |
| 14.0 | 1.100000 | F122 | 1.1066166708 | **1.1073745030** | 170 | 1.10 | +0.0685% | +0.6704% |
| 21.0 | 1.091400 | F122+PW | 1.0982453728 | **1.0993840223** | 170 | 1.10 | +0.1037% | +0.7315% |
| 28.0 | 1.088933 | PW | 1.0942036229 | **1.0951787825** | 170 | 1.20 | +0.0891% | +0.5736% |
| 34.5 | 1.087500 | PW | 1.0918503685 | **1.0928294757** | 170 | 1.10 | +0.0897% | +0.4901% |

**All nine rows that appear in CQH Table 1 beat its best column**, and the tenth
(`A = 36/11`) beats CPdL's own published `1.1965` by 0.19%. That `A = 36/11` row
was not searched for anything published since CPdL 2021. Against this
project's own linear program it is **6 wins and 3 losses**: the SDP wins at
`A = 1` and for `A >= 10`, and loses by two or three hundredths of a percent at
`A = 1.5, 4, 6`.

**The headline row.**

```
C_+(28)  >=  1.0951787824934853        (certified, ARB)
previous certified (LP, FREEDOMS.md)    1.0942036228559691
CQH's rigorous record 1/0.91833       = 1.0889331721712239
proven upper bound (UPPER.md)           1.1136

2/C_+(28) <= 1.8261858538   (was 1.8278133596; CQH print 1.837 h(-D))
```

In CQH's Corollary 5 the constant `2 h(-D)/C_+(28)` improves from
`1.8278134 h(-D)` to `1.8261859 h(-D)`. For `f(u,v) = u^2 + 27 v^2`, where
`h(-108) = 3`, the constant `5.4834` becomes `5.4786`.

**The row that was losing.** `EXTREMAL2.md` records `A = 1` as the one row of 68
where the bandlimited linear program loses to the published record, by 3.08%,
and diagnoses it: at `A = 1` the best function wants a very wide band, and a
wide band spread over a capped number of harmonics cannot resolve it. The
Gaussian family has no band at all, so that obstruction does not exist. The SDP
certifies

```
C_+(1)  >=  1.9903625363954200      = 99.518% of the exact value
exact (CMS Theorem 2(b))      C_+(1) = 2
published record                     1.9602   (+1.5387%)
LP certified                         1.8998323 (+4.7652%)
```

This row was never a standing record — CMS prove `C_+(1) = 2` — but it is the
row this project was recorded as losing, and it is no longer lost.

## Why

Two sentences.

**Why the SDP wins at large `A`:** it does not, at the objective it optimises —
its objective is still `1.088` at `A = 28` where the linear program certifies
`1.094`. It wins because the *function* it returns is worth `1.095`, and the
`7e-3` difference between those two numbers is the slack in Lemma 4's two
inequalities, which nobody was recovering because the published pipeline reports
the relaxed number. Recovering it costs one certificate run per row.

**Why the SDP wins at `A = 1`:** the bandlimited family has a bandwidth `B`, and
at small `A` the optimum wants `B` large, which spreads a fixed number of
harmonics too thin. The Gaussian family has no bandwidth. The dilation parameter
of Section 6.1 replaces it and is much cheaper to push: `c = 0.45` costs the
same as `c = 1`.

**Why the linear program still wins at `A = 1.5, 4, 6`:** the losses there are
`2e-4` to `4e-4`, about the size of the scatter between neighbouring degrees, so
this is a statement about how hard each side was pushed, not a structural one.
Those rows were solved at degree 122–146 because the solver degrades above that;
the winning rows are at degree 170–202.

## Novelty

**Not new.** The extremal problem (CMS 2017). The application and the constant 28
(CQH 2022). The four-function convex reformulation, the Gaussian-times-polynomial
restriction, the sum-of-squares certificate, the Laguerre basis, and the use of
ARB to certify the result — all Chirre–Pereira Júnior–de Laat 2021, Section 4.
Using this machinery on CQH's Table 1 is what de Laat already did for them.

**New.** The ten certified numbers in Section 8, and the coefficient file behind
each. The measurement of Section 6.2 — that the SDP's own objective sits `5e-3`
to `7e-3` below the certified value of the function it returns, and that this,
not the function class, is why the published `F` columns look weak at large `A`.
The dilation parameter of Section 6.1. And the certificate itself: the closed-form
Hermite antiderivative and the Taylor-model sign sweep are, as far as I can find,
not written down for this problem — but they are conveniences for this basis, not
new mathematics, and the *idea* of certifying in ARB is CPdL's.

## Caveats

* **I did not reproduce a published SDP value exactly.** At matched degree I am
  `4e-4` to `5e-4` low, consistently (Section 5.3). The cause is unexplained.
  This weakens the claim that my formulation is identical to theirs; what it
  does not weaken is the certified numbers, which stand on the certificate and
  not on the formulation being anyone's in particular.
* **Ten rows, not 68.** CQH's Table 1 has 68 rows. Nine of them are done here
  (plus `A = 36/11`). The other 58 are untouched, and `EXTREMAL2.md`'s certified
  values there still stand as the project's best.
* **Small.** The wins over the linear program are between `+0.02%` and `+0.10%`
  except at `A = 1`. The wins over the published record are between `+0.25%` and
  `+0.73%`.
* **The solver is the ceiling, not the method.** CLARABEL reports
  `optimal_inaccurate` above about `d = 45` and crashed once at `d = 60`. The
  SDP objective is still climbing there (Section 5.4: `1.0944` at `d = 60`
  against `1.0881` at `d = 42`), and the honest value of the returned function
  scatters by about `3e-4` between neighbouring degrees, which is the same size
  as the three losses. A multiprecision solver (`sdpa-gmp`, as CPdL used) would
  settle whether those three rows are really losses.
* **The float64 ranking is not always reliable.** The scan ranks candidates by a
  float64 recomputation on a grid. At a few small-`A` settings the coarse and
  fine grids disagreed by up to `3e-2` (for example `A = 1`, `d = 30`,
  `c = 1.4`: 1.86616 against 1.89594). Candidates are therefore ranked by the
  smaller of the two, and the top two per row are certified and the better
  certified one kept. So a bad ranking can only cost a better answer; it cannot
  produce a wrong one.
* **The `c` grid is coarse.** Five to nine values per row, and on several rows
  the winner sat at or next to the edge of the grid. Those rows are probably not
  the best this method can do.
* **What "certified" rests on.** ARB's ball arithmetic, in particular `erf`,
  `exp`, `sqrt` and the incomplete gamma function; and the human reading of the
  problem statement out of the CMS and CQH PDFs. Both assumptions are the same
  ones `EXTREMAL2.md` and CPdL make.
* **The SDP objective alone is not rigorous.** CLARABEL satisfies the equality
  constraints only to about `1e-7`, so its reported objective can sit a hair
  above the true relaxed value (visible at `A = 1`, where the gap in Section 6.2
  is `-2.4e-7`). CPdL handle this with a Gershgorin argument that repairs the
  identity exactly. That step is **not** implemented here, and it is not needed:
  nothing in Section 8 uses the SDP objective. Every reported number is the
  certificate's evaluation of the returned `F` itself, and that is valid for any
  coefficients whatever.

## Reproduce

All in `weil-form/`. Write long runs to a file and read the file; piping a
`python -u` run into `tail` loses the output.

```
python sdp_check.py 30                          # the five basis identities
python cplus_sdp.py 20 28                       # one solve: d=20, A=28, c=1
python cplus_sdp.py 42 28 1.2                   # the headline solve

# scans (append to a .jsonl and write one coefficient file per candidate)
python -u cplus_sdp_scan.py 28 20,30,40,42 1.0        > sdp_A28.log 2>&1
python -u cplus_sdp_scan.py 28 42 0.9,1.0,1.1,1.2,1.35 sdp_mid.jsonl
python -u cplus_sdp_scan.py 1.0 40,50 0.45,0.55,0.7   sdp_smallA.jsonl

# certify one file, and the calibration
python cplus_sdp_cert.py sdpfun_A28.0_d42_c1.2.txt 28 1.2
python cplus_sdp_cert.py sdpfun_gauss.txt 28 1.0      # -> 0.6587112988255192

# an independent 40-digit evaluation of the same file
python cplus_sdp_mp.py sdpfun_A28.0_d42_c1.2.txt 28 1.2

# pick and certify the best candidate per A, then build the table
python cplus_sdp_best.py cplus_sdp_runs.jsonl,sdp_smallA.jsonl,sdp_mid.jsonl 2
python cplus_sdp_report.py
```

**Files.** `sdp_basis.py` (the two bases and the moment matrices),
`sdp_check.py` (validation), `cplus_sdp.py` (the SDP and the float64
recomputation), `cplus_sdp_scan.py` (degree and dilation scans),
`cplus_sdp_cert.py` (the ARB certificate), `cplus_sdp_mp.py` (independent
mpmath evaluation), `cplus_sdp_best.py` (rank and certify),
`cplus_sdp_report.py` (the tables above). Test functions:
`sdpfun_A<A>_d<d>_c<c>.txt`, one Hermite coefficient per line; each file is a
complete, self-contained claim.

**Environment variables** for the certificate: `PREC` (bits, default 300),
`JORD` (Taylor order, default 6), `NSEED` (seed intervals, default 4000),
`FLOORBITS`, `XCUT`, `TAILTGT`, `QUIET`.
