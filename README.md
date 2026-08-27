# Certified lower bounds for `C_+(A)`

https://github.com/sallerk/certified-cplus-bounds

This package contains everything needed to check, independently, that

```
C_+(28) >= 1.0951787824934853
```

and the corresponding bound at each of the other 67 values of `A` tabulated by
Chirre and Quesada-Herrera. Every bound here is **certified** — computed end to
end in interval (ball) arithmetic, so each is a proved inequality rather than a
numerical estimate.

The consequence for prime gaps in binary quadratic forms:

```
2/C_+(28) <= 1.8261858538 * h(-D)        (previously 1.837)
```

## The problem

Carneiro, Milinovich and Soundararajan (*Fourier optimization and prime gaps*,
arXiv:1708.04122), Extremal Problem 2, eq. (1.3):

```
C_+(A) = sup over even, continuous, real F in L^1(R), F nonzero, of

             F(0) - A * integral over |t|>1 of (Fhat(t))_+ dt
             -------------------------------------------------
                              ||F||_1
```

with `Fhat(t) = integral of exp(-2*pi*i*x*t) F(x) dx`. Any admissible `F`
gives a valid lower bound, so **a lower bound is verifiable from the test
function alone**. That is what this package is: 136 test functions and the
programs that evaluate them rigorously.

## Requirements

```
pip install python-flint numpy
```

`python-flint` provides the ARB ball-arithmetic bindings. Developed against
python-flint 0.9.0 and numpy 2.x. Nothing else is needed — the three
certificate programs import only these and the standard library.

## Verify one row

```
python certificates/cplus_certA.py functions/lp/testfunction_A28.0_r27raw.txt        # LP function
python certificates/cplus_sdp_cert.py functions/sdp/sdpfin_A28.0_d42_c1.2.txt 28 1.2  # SDP function
```

Each prints its intermediate enclosures and then the certified value. A run
takes a few seconds.

## Verify everything

```
python verify_all.py          # all 68 rows
python verify_all.py 28       # just A = 28
```

For each row this re-certifies the bound from its coefficient file and checks
three things: that it matches the recorded value, that it strictly exceeds the
published record, and that it does **not** exceed the certified upper bound.
It reports a problem count and exits non-zero if anything fails.

## What the certificate actually proves

Three places where a certification of this kind usually leaks, and how each is
handled:

- **No quadrature error exists anywhere.** All the integrals are closed-form
  (`integral of sinc = Si/pi` for the band-limited functions; a Hermite
  antiderivative for the Gaussian ones), so there is nothing to estimate. The
  only real work is deciding the sign of a function.
- **Sign decisions are conservative.** A subinterval is declared sign-definite
  only when its enclosure excludes zero. Subintervals that remain undecided are
  charged at `width * max|F|`. So no claim about having located every root is
  ever required — a missed sign change can only weaken the bound, never
  invalidate it.
- **The tail is proved, not observed.** An explicit `|F(x)| <= K/x^3` bound is
  derived from the coefficients, with `K` written down.

The certificates also refuse false claims. `cplus_dual_cert.py` takes a target
and declines if it cannot prove it, reporting the subinterval where it failed:

```
python certificates/cplus_dual_cert.py functions/dual/multiplier_A28.txt 1.1136   # certifies
python certificates/cplus_dual_cert.py functions/dual/multiplier_A28.txt 1.1140   # correctly refuses
```

## Comparing against the published table

`cqh_table1.txt` holds all 68 rows of Chirre–Quesada-Herrera Table 1, read from
their PDF. It has three columns of published lower bounds: `F82` and `F122`
(semidefinite programming, polynomial degree at most 82 and 122) and `PW` (their
band-limited function).

**The record for each row is the maximum across those three columns, and which
column wins changes three times** — `PW` at `A = 1`, `F122` for `A = 1.5` to
`20.5`, an exact tie at `A = 21.0`, and `PW` from `A = 21.5` to `34.5`.
Comparing against a single fixed column would misstate the comparison on 28 rows.

## Files

| path | contents |
|---|---|
| `verify_all.py` | re-certifies and checks every row; start here |
| `final_table.json` | the consolidated results, one record per row |
| `cqh_table1.txt` | the published records being compared against |
| `full_verify.txt` | log of a complete 68-row run |
| `certificates/` | the three ARB certificate programs |
| `functions/lp/` | band-limited test functions, one per row |
| `functions/sdp/` | Gaussian-times-polynomial test functions, one per row |
| `functions/dual/` | dual multipliers, giving the upper bounds |
| `docs/note.html` | the short note; open in a browser |
| `docs/EXTREMAL2.md` | full method and the complete certified table |
| `docs/UPPER.md` | certified upper bounds, and what they say about the remaining gap |
| `docs/RECORD_CHECK.md` | the literature check, and what could not be checked |
| `docs/SDP.md` | the semidefinite formulation and its calibration |

## What is and is not new

**None of the mathematics is new.** The extremal problem is Carneiro,
Milinovich and Soundararajan's. The constants, the table and the application to
primes represented by quadratic forms are Chirre and Quesada-Herrera's. The
semidefinite formulation is Chirre, Pereira Júnior and de Laat's. Certifying in
ARB is what both groups already did.

What is new here is only the 136 test functions and their certificates — found
by running a linear program and a semidefinite program over 100–200 coefficients
each, in place of the hand-tuned three- and four-parameter families that produced
the published records, and keeping whichever of the two did better on each row.

## Known limitations

- The SDP pipeline does **not** exactly reproduce published SDP values at
  matched degree; it lands 4–5e-4 low, and an off-by-one degree convention was
  tested and rejected. The discrepancy is unexplained. It does not affect the
  bounds, which are certified by evaluating the resulting functions directly
  rather than by trusting any solver's own objective. Details in `SDP.md`.
- A solver's reported objective is not the value of the function it returns. At
  `A = 28` the SDP objective reads 1.0881 while the function it produces is
  worth 1.0952.
- The literature check found no improvement since 2022 across three citation
  indexes, author sweeps and eleven downloaded papers — but MathSciNet, zbMATH
  and Google Scholar were unreachable, and the indexes demonstrably undercount.
  See `RECORD_CHECK.md`.
- Correctness rests on ARB, as does the record being displaced.
