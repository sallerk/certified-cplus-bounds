# Verdict: 1.0889 is still the published record for C_+(28), and 1.1965 is still the record for C_+(36/11)

Checked on **2026-08-26**. Nothing published between the Chirre–Quesada-Herrera paper
(2022) and today improves `C_+(A)` at any value of `A`. The lower bounds in
`EXTREMAL.md` and `EXTREMAL2.md` are therefore still improvements over the literature,
so far as three independent citation indexes and a direct sweep of the authors' arXiv
listings can show.

This file replaces the "Not verified to be the current record" caveat in the
"What this is NOT" section of `EXTREMAL2.md`, and the matching caveat in `EXTREMAL.md`,
with a documented check. It does not make the caveat vanish — see
[What I could not check](#what-i-could-not-check).

## 1. What the record is

The quantity is Carneiro–Milinovich–Soundararajan's **Extremal Problem 2**
(arXiv:1708.04122, eq. 1.3):

```
C_+(A) = sup over F in A_+, F != 0, of
         [ F(0) - A * int_{|t|>1} (Fhat(t))_+ dt ] / ||F||_1

where 1 <= A < infinity and A_+ is the class of even, continuous
F : R -> R with F in L^1(R).  The paper writes the domain of the
integral as [-1,1]^c, which is the set |t| > 1.
```

Three values of `A` matter here.

| constant | published record | who, when | what it controls |
|---|---|---|---|
| `C_+(4)` | **1.17233** | Chirre–Pereira–de Laat, *Primes in arithmetic progressions and semidefinite programming*, Math. Comp. **90** (2021) 2235–2246, Theorem 1 | gaps between primes in an arithmetic progression, under GRH; the constant is `phi(q)/C_+(4)`, giving their `0.8531 phi(q)` |
| `C_+(36/11)` | **1.1965**, stated rigorously as `1/0.8358 = 1.1964585` | Chirre–Pereira–de Laat, op. cit., Section 1.3 | gaps between ordinary primes under RH |
| `C_+(28)` | **1.0889**, stated rigorously as `1/0.91833 = 1.0889332` | Chirre–Quesada-Herrera, *Fourier optimization and quadratic forms*, Q. J. Math. **73** (2022) 539–577, eq. (6.10)–(6.11) | gaps between primes represented by a binary quadratic form, under GRH |

All three are lower bounds. Anyone can beat them by exhibiting a better test
function. This package does, at all three: `C_+(4) >= 1.1750163119`,
`C_+(36/11) >= 1.1991189437902487` and `C_+(28) >= 1.0951787825`, each certified
in ARB. The `36/11` row is not part of CQH Table 1 and is kept in
`extra_rows.json`; it lowers the ordinary-prime gap constant from `0.8358` to
`0.8340` under RH.

`C_+(4)` is also the `J(F)` of CMS eq. (6.11) — their `22/25 = 0.88` interval for
ordinary primes comes from `A = 4` — but for ordinary primes that route was
superseded by the `36/11` one above. Of the three, only `A = 4` and `A = 28` are
rows of CQH Table 1; `36/11 = 3.2727...` falls between grid points. **At `A = 4`
the CPdL value exceeds anything in Table 1** (`1.1673`), so the record for that
row is theirs, not CQH's. This is recorded in `external_records.json` and is
applied by `verify_all.py`.

## 2. How I searched

The key structural fact is that **any paper improving `C_+(A)` would have to cite CMS**,
because CMS is where the problem is defined and where its name comes from. So the
citation list of CMS is a superset of the candidates. I pulled that list from three
indexes that build their citation graphs differently, then downloaded and machine-read
every plausible candidate.

| # | Source | Query | Result |
|---|---|---|---|
| 1 | Semantic Scholar API | citations of `arXiv:1708.04122` (CMS) | 43 citing works listed, 2017–2026 |
| 2 | Semantic Scholar API | citations of `arXiv:2012.07781` and of `DOI:10.1093/qmath/haab041` (CQH) | 3 citing works |
| 3 | Semantic Scholar API | citations of `arXiv:2005.02393` and `DOI:10.1090/mcom/3638` (Chirre–Pereira–de Laat) | 3 citing works |
| 4 | OpenAlex API | citations of CMS (`W2749290691`) with `from_publication_date:2025-01-01` | 11 works, all listed below |
| 5 | OpenCitations / Crossref | citations of `10.4171/cmh/467` (CMS journal version) | 30 citing DOIs |
| 6 | OpenCitations / Crossref | citations of `10.1093/qmath/haab041` (CQH journal version) | 5 citing DOIs |
| 7 | OpenCitations / Crossref | citations of `10.1090/mcom/3638` (CPdL journal version) | 2 citing DOIs |
| 8 | arXiv API | `all:"Fourier optimization"`, newest first, 100 results | 11 hits; the newest relevant one is Feb 2025. **No 2026 arXiv paper has "Fourier optimization" in its title or abstract.** |
| 9 | arXiv API | author listings for Chirre, Quesada-Herrera, de Laat, Carneiro, Milinovich, Gonçalves, Ortega-Cerdà, and Antonio Pedro Ramos | full lists to Aug 2026; nothing on `C_+(A)` |
| 10 | arXiv API | `(extremal OR Beurling-Selberg OR semidefinite) AND (prime gaps OR prime numbers OR Riemann hypothesis)`, newest first | 20 entries from 2024 on; none relevant |
| 11 | arXiv API | `cat:math.NT AND "quadratic form" AND "gaps between"` | one hit since 2022 (Żmija, below), not this problem |
| 12 | Web search | direct searches on the strings `C_+(28)`, `1.0889`, `1.837`, `0.8358` | no improvement found; results were noise |
| 13 | Direct read | 11 candidate PDFs downloaded and grepped for `C+(`, `C_+(`, `1.0889`, `0.91833`, `0.8358`, `1.837` | see next section |

## 3. What the candidates actually say

I downloaded each plausible paper and searched its text. "No `C_+`" below means the
string `C+(` does not occur anywhere in the paper — CQH or CMS appear only in the
bibliography.

| Paper | Date | Verdict |
|---|---|---|
| Quesada-Herrera, *Fourier optimization and consequences of GRH*, arXiv:2411.05095 | Nov 2024 | **The key witness.** An expository article by one of the two CQH authors. It states "we have the lower bound `C_+(28) >= 1.0889`" and prints the same three-term test function with `{68, 5, 1}` and `lambda = 0.98644`. For ordinary primes it says Chirre–Pereira–de Laat "replaced 21/25 by 0.8358" and that "this is still the best approach to date". |
| Carneiro–Milinovich–Quesada-Herrera–Ramos, *Fourier optimization, the least quadratic non-residue, and the least prime in an AP*, arXiv:2404.08380, Math. Comp. (2025) | v2 Aug 2025 | No `C_+`. It sets up its own, different extremal problems for a different application. |
| Carneiro–Milinovich–Ramos, *Fourier optimization and Montgomery's pair correlation conjecture*, arXiv:2310.01913, Math. Comp. (2024) | Oct 2023 | No `C_+`. Different extremal problem. |
| Das–Ismoilov–Ramos, *Fourier optimization and pair correlation problems*, arXiv:2502.05106 | Feb 2025 | No `C_+`. Two extremal problems, neither is this one. |
| Chirre–Dimitrov–Quesada-Herrera–Sousa, *An extremal problem and inequalities for entire functions of exponential type*, arXiv:2304.05337, Proc. AMS (2024) | Apr 2023 | No `C_+`. |
| Ortega-Cerdà, *Some Fourier Extremal Problems*, ICM 2026 proceedings, SIAM, DOI 10.1137/25m1804261 | Jul 2026 | **The most recent survey of the field.** It is about point evaluation in Paley–Wiener spaces and the Hörmander–Bernhardsson function. It cites CMS as an application but gives no numerical bound for `C_+(A)`. |
| Gonçalves–Radchenko–Ramos, *The Hörmander–Bernhardsson function in higher dimensions*, arXiv:2608.22198 | Aug 2026 | CMS in the bibliography only. Different extremal problem. |
| Bondarenko–Ortega-Cerdà–Radchenko–Seip, *The Hörmander–Bernhardsson extremal function: a preliminary study*, arXiv:2504.05205 | 2024/2025 | CMS in the bibliography only. Different extremal problem. |
| Chamberland–Straub, *Weakening the Legendre Conjecture*, arXiv:2602.22502, Amer. Math. Monthly (2026) | 2026 | Uses CMS. Quotes the **unimproved** CMS explicit constant `22/25` for "a prime in `[x, x + (22/25) sqrt(x) log x]`". Does not know of anything better. |
| Balanzario, *Clusters and deserts of prime numbers*, arXiv:2411.05932 | Nov 2024 | No `C_+`, no constants of this kind. |
| Johnston–Kerr, *The infinitude of square-free palindromes*, arXiv:2601.07097 | Jan 2026 | CMS in the bibliography only. |
| Żmija, *Large gaps between values of several binary quadratic forms*, arXiv:2509.15365 | Sep 2025 | Gaps between the **values** of quadratic forms, unconditional, Erdős-type. Not primes, not GRH, not this constant. |

The remaining works in the CMS citation list are visibly off-topic by title
(Bernstein–Nikolskii inequalities, sphere packing, cyclotomic polynomial coefficients,
Mills' constant, elliptic curves over Hasse pairs, explicit prime-counting estimates,
Jacobi transforms, sign uncertainty principles, and a general RH survey). None of them
is a Fourier-optimization computation of `C_+(A)`.

**One negative result worth stating positively.** Between them, Chirre,
Quesada-Herrera, Carneiro, Milinovich, de Laat, Gonçalves, and Ramos have published
about thirty papers since 2022. None of them revisits `C_+(A)`. The people who own
this problem have moved on to other extremal problems — pair correlation, the least
quadratic non-residue, de Branges spaces, the Hörmander–Bernhardsson function.

## 4. Semidefinite programming specifically

The task asked whether a new SDP computation of these constants exists. David de Laat
computed the `F82` and `F122` columns of CQH's Table 1. Since then his number-theory
SDP work is Gonçalves–de Laat–Leijenhorst, *Multiplicity of nontrivial zeros of
primitive L-functions via higher-level correlations* (arXiv:2303.01095, Mar 2023) —
a different problem. His other recent papers are about spherical codes, the `D_4`
root system, equiangular lines, and sphere packing. There is no new SDP run at
`C_+(A)` in the literature.

Note also that a bigger SDP would not obviously help here. CQH's own Conjecture 11
says the `P(x) exp(-x^2)` family is asymptotically the *wrong* family as `A` grows;
their Table 1 shows the bandlimited column overtaking the SDP column at `A = 21`.
At `A = 28` the SDP column is already 0.2% behind.

## Verdict

1. **`C_+(28) >= 1.0889331722` (Chirre–Quesada-Herrera, 2022) is still the published
   record.** No improvement found anywhere in 2022–2026. This project's certified
   `C_+(28) >= 1.0940451158772697` beats it by 0.469%, and the claim to have beaten
   the record stands.

2. **`C_+(36/11) >= 1.1964585` (Chirre–Pereira–de Laat, 2021) is still the published
   record.** No improvement found. A November 2024 survey calls it "still the best
   approach to date", and a 2026 *Monthly* article quotes the older, weaker CMS
   constant. This project's `EXTREMAL.md` value 1.1986243 beats it by 0.18%.

3. **`C_+(A)` has not been improved at any other value of `A` either.** The best
   published value for each row of CQH's Table 1 is still the maximum of that row's
   three columns, which is what `weil-form/cqh_table1.txt` now records.

4. **The prime-gap constant `1.837 h(-D)` for binary quadratic forms has not been
   improved** by anyone. It is `2/C_+(28)`, so any improvement would have had to come
   through `C_+(28)`, and none did.

5. **Nothing new was found that is relevant.** The one genuinely new item is
   Ortega-Cerdà's ICM 2026 proceedings survey *Some Fourier Extremal Problems*
   (July 2026), which is more recent than the Quesada-Herrera survey the earlier
   claim rested on. It is about a neighbouring extremal problem and does not update
   any `C_+(A)` number.

## What I could not check

Listed honestly, worst first.

1. **MathSciNet and zbMATH.** MathSciNet is behind a paywall I have no access to.
   zbMATH Open returned HTTP 403 to my fetcher. These are the two indexes a
   professional would check first for a "has this been superseded" question, and I
   checked neither. Google Scholar's "cited by" was also unreachable.

2. **One Russian-language paper I could not read.** A. D. Manov, *On some extremal
   problems for entire functions of exponential type* (Чебышевский сборник /
   Chebyshevskii Sbornik **26** (2025), no. 1, 47–61, DOI
   `10.22405/2226-8383-2025-26-1-47-61`) cites CMS. I could not retrieve its abstract
   or full text from mathnet.ru or by search. Manov's other work is on positive
   definite functions and Turán-type problems, so it is very unlikely to contain a
   `C_+(A)` computation — but that is an inference from an author's track record, not
   a reading of the paper. **I cannot rule it out.**

3. **One low-profile journal article.** D. Lattanzi, *Statistical Distributions of
   Prime Number Gaps*, J. Adv. Math. Comput. Sci. **39** (2024),
   DOI `10.9734/jamcs/2024/v39i11861`, cites CMS. Not retrieved. The title is about
   empirical statistics, not Fourier optimization.

4. **Citation indexes disagree with each other, and all of them undercount.**
   Semantic Scholar lists only 3 works citing CQH, but I verified by direct reading
   that Chirre–Dimitrov–Quesada-Herrera–Sousa (arXiv:2304.05337) cites CQH and is
   **not** in that list. OpenAlex reports CQH's citation count as 1, which is plainly
   wrong. So no single index here is complete, which is exactly why I used three plus
   direct author sweeps. The union is much better than any one of them, but I cannot
   claim it is exhaustive.

5. **Grey literature.** PhD theses, conference talks, unpublished SDP runs, and
   papers not yet posted. A. P. Ramos and others in this group could easily have
   sharper numbers in a thesis. Not searched.

6. **Papers that improve `C_+(A)` without citing CMS or CQH.** Structurally almost
   impossible, since the problem has no other name, but not formally excluded by
   anything I did.

7. **Errata.** I read the arXiv v2 of CQH, not the published Q. J. Math. version, and
   I did not check for a published correction to either paper.

## Appendix: CQH Table 1

Extracted from the arXiv PDF with `pdftotext -layout` and saved as
`weil-form/cqh_table1.txt`. All 68 rows, `A = 1.0` to `A = 34.5` in steps of 0.5.

The three columns are three different families of test functions, all giving lower
bounds on the same `C_+(A)`:

* **F82** — semidefinite programming with `F(x) = P(x) exp(-x^2)`, `deg P <= 82`,
  computed for CQH by David de Laat.
* **F122** — the same with `deg P <= 122`.
* **PW** — bandlimited, `F(x) = H(x/lambda)` with
  `H(x) = cos(2 pi x)[a1/(1-16x^2) + a2/(9-16x^2) + a3/(25-16x^2)]`; the parameters
  are CQH's Table 2.

**The record for a row is the largest of the three, and which column that is varies by
row.** `F122 >= F82` always, so the SDP record is always `F122`. `F122` and `PW` cross
at `A = 21.0`, where both equal 1.0914. Below that the SDP column wins (at `A = 8.5`,
1.1166 against 1.1073); above it the bandlimited column wins (at `A = 25.0`, 1.0898
against 1.0883; at `A = 28.0`, 1.0889 against 1.0865). The `A = 1.0` row is not a
record at all — CMS prove `C_+(1) = 2` exactly, and CQH include the row only for
comparison.

## Sources

Primary:

* Carneiro, Milinovich, Soundararajan, *Fourier optimization and prime gaps*,
  arXiv:1708.04122, Comment. Math. Helv. **94** (2019).
  <https://arxiv.org/abs/1708.04122> — DOI `10.4171/cmh/467`
* Chirre, Quesada-Herrera, *Fourier optimization and quadratic forms*,
  arXiv:2012.07781, Q. J. Math. **73** (2022) 539–577.
  <https://arxiv.org/abs/2012.07781> — DOI `10.1093/qmath/haab041`
* Chirre, Pereira Júnior, de Laat, *Primes in arithmetic progressions and semidefinite
  programming*, arXiv:2005.02393, Math. Comp. **90** (2021) 2235–2246.
  <https://arxiv.org/abs/2005.02393> — DOI `10.1090/mcom/3638`

The two surveys that pin the dates:

* Quesada-Herrera, *Fourier optimization and consequences of the generalized Riemann
  hypothesis*, arXiv:2411.05095 (7 Nov 2024). <https://arxiv.org/abs/2411.05095>
  — states `C_+(28) >= 1.0889` and that CPdL's 0.8358 is still the best approach.
* Ortega-Cerdà, *Some Fourier Extremal Problems*, ICM 2026 proceedings, SIAM for the
  IMU (July 2026), DOI `10.1137/25m1804261`.
  <https://epubs.siam.org/doi/10.1137/25M1804261> — no `C_+(A)` numerics.

Candidates read and cleared:

* <https://arxiv.org/abs/2404.08380> — Carneiro, Milinovich, Quesada-Herrera, Ramos
* <https://arxiv.org/abs/2310.01913> — Carneiro, Milinovich, Ramos
* <https://arxiv.org/abs/2502.05106> — Das, Ismoilov, Ramos
* <https://arxiv.org/abs/2304.05337> — Chirre, Dimitrov, Quesada-Herrera, Sousa
* <https://arxiv.org/abs/2608.22198> — Gonçalves, Radchenko, Ramos
* <https://arxiv.org/abs/2504.05205> — Bondarenko, Ortega-Cerdà, Radchenko, Seip
* <https://arxiv.org/abs/2602.22502> — Chamberland, Straub
* <https://arxiv.org/abs/2411.05932> — Balanzario
* <https://arxiv.org/abs/2601.07097> — Johnston, Kerr
* <https://arxiv.org/abs/2509.15365> — Żmija
* <https://arxiv.org/abs/2303.01095> — Gonçalves, de Laat, Leijenhorst (by title only)

Indexes used:

* Semantic Scholar Graph API — <https://api.semanticscholar.org/graph/v1/>
* OpenAlex — <https://api.openalex.org/>
* OpenCitations COCI / Crossref — <https://api.opencitations.net/index/v1/>
* arXiv API — <http://export.arxiv.org/api/query>
