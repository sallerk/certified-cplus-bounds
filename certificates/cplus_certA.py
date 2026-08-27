"""Certified lower bound on C_+(A) for ANY saved test function file.

This is `cplus_cert.py` with the three hard-coded constants (A = 28, B = 63/50,
71 coefficients) read from the file instead, so one script certifies every row
of the table.  The mathematics is identical; see cplus_cert.py's docstring and
the Certification section of EXTREMAL2.md for what each step proves.

Every number that feeds the answer is an ARB ball (python-flint), so the printed
bound is a rigorous enclosure and we report its lower end.

  C_+(A) = sup_F [ F(0) - A int_{|t|>1} (Fhat(t))_+ dt ] / ||F||_1

The file holds even harmonics m = 2k only, so with k = m/2:

  Fhat(t) = sum_k a_k cos(pi k t / B)   on |t| <= B,  0 outside
  F(x)    = B sum_k a_k [ sinc(u-k) + sinc(u+k) ],  u = 2Bx

Exact identities used (all independent of A, B and the degree):

  F(0)     = 2 B a_0
  Fhat(B)  = sum_k (-1)^k a_k        -> forced to exactly 0 by a rational shift
  Fhat'(B) = 0                       -> automatic: a sum of sin(pi k)
  F        = (2 B u / pi) sin(pi u) sum_k (-1)^k a_k/(u^2-k^2)
           = (2 B sin(pi u)/(pi u^3)) [ T2 + sum_k (-1)^k a_k k^4/(u^2-k^2) ]
  int F dx = (1/(2 pi)) sum_k a_k [ Si(pi(u-k)) + Si(pi(u+k)) ]

B is taken as the EXACT decimal in the file (1.26 -> 63/50), not as the float64
nearest it.  That is a definition of the function being certified, not an
approximation: any admissible F gives a valid lower bound, so we are free to
declare the exact one.

Usage:   python cplus_certA.py testfunction_A28.txt
Env:     UEXACT (cutoff u_c, default 400), NOPROJ=1 (diagnostic), DELTA (force
         the width of the endpoint interval [B-DELTA, B]), QUIET=1.
"""
import sys, os, time, json
from fractions import Fraction as Fr
from flint import arb, ctx, fmpq

ctx.dps = 60
PI = arb.pi()
U_EXACT = int(os.environ.get('UEXACT', 400))
NOPROJ = os.environ.get('NOPROJ') == '1'
QUIET = os.environ.get('QUIET') == '1'
MINW = Fr(1, 2**41)
# widths of the endpoint interval [B-delta, B], tried in this order.  0.014 is
# the value hard-coded in cplus_cert.py, so an A = 28 run reproduces it exactly.
DELTAS = [Fr(14, 1000), Fr(30, 1000), Fr(60, 1000), Fr(120, 1000)]
if os.environ.get('DELTA'):
    DELTAS = [Fr(os.environ['DELTA'])]


def say(*a, **k):
    if not QUIET:
        print(*a, **k)


def sweep(lo0, hi0, fun, minw, evcap=None):
    """Adaptive ball sweep.  Leaves (lo, hi, s), s = +1/-1 certified sign,
    s = 0 undecided at width minw.  Returns (leaves, nevals) or (None, nevals)
    if the evaluation budget is blown."""
    stack, out = [(lo0, hi0)], []
    n = 0
    while stack:
        lo, hi = stack.pop()
        v = fun(lo, hi)
        n += 1
        if evcap is not None and n > evcap:
            return None, n
        if v.lower() > 0:
            out.append((lo, hi, 1))
        elif v.upper() < 0:
            out.append((lo, hi, -1))
        elif hi - lo <= minw:
            out.append((lo, hi, 0))
        else:
            mid = (lo + hi) / 2
            stack.append((mid, hi))
            stack.append((lo, mid))
    out.sort()
    return out, n


def Q(f):
    return arb(fmpq(f.numerator, f.denominator))


def fq(f):
    return fmpq(f.numerator, f.denominator)


# ------------------------------------------------------------- 1. coefficients
def load(path):
    """Returns (a dict k->Fraction, B exact Fraction, A exact Fraction)."""
    a, Bf, Af = {}, None, None
    for line in open(path):
        s = line.strip()
        if 'eq. (1.3).' in s and 'A =' in s:
            Af = Fr(s.split('A =')[1].strip())
        if s.startswith('B ='):
            Bf = Fr(s.split('=')[1].strip())        # exact decimal, e.g. 63/50
        if s.startswith('a['):
            m = int(s[s.index('[') + 1:s.index(']')])
            assert m % 2 == 0, 'this certificate handles the even-harmonic basis only'
            a[m // 2] = Fr(float(s.split('=')[1]))  # float64 -> exact dyadic
    return a, Bf, Af


path = sys.argv[1] if len(sys.argv) > 1 else 'testfunction_A28.txt'
a_raw, B, A_PARAM = load(path)
KMAX = max(a_raw)
NK = KMAX + 1
assert set(a_raw) == set(range(NK)), 'coefficient indices must be 0..KMAX'
assert B > 1 and U_EXACT > KMAX

say('== 1. exact coefficients ==')
say('file %s :  A = %s,  B = %s (exact),  %d even harmonics m = 0..%d'
    % (path, A_PARAM, B, NK, 2 * KMAX))
resid = sum((-1)**k * a_raw[k] for k in range(NK))
say('Fhat(B) before projection  = %.3e' % float(resid))
a = dict(a_raw) if NOPROJ else {k: a_raw[k] - Fr((-1)**k) * resid / NK for k in range(NK)}
if NOPROJ:
    say('*** NOPROJ diagnostic run: constraint NOT imposed, F is not in L^1 ***')
else:
    assert sum((-1)**k * a[k] for k in range(NK)) == 0
    say('Fhat(B) after  projection  = 0  (exact in rational arithmetic)')
say("Fhat'(B) = -sum_k a_k (pi k/B) sin(pi k) = 0  (exact, structural)")
say('max |a_new - a_old|        = %.3e' % max(abs(float(a[k] - a_raw[k])) for k in range(NK)))
T2 = sum(Fr((-1)**k * k * k) * a[k] for k in range(NK))
S4 = sum(abs(a[k]) * Fr(k**4) for k in range(NK))
say('T2 = sum (-1)^k k^2 a_k    = %.10e' % float(T2))

aA = [Q(a[k]) for k in range(NK)]
BA = Q(B)
AA = Q(A_PARAM)
T2A = Q(T2)
sg = [arb((-1)**k) for k in range(NK)]
k4 = [arb(k**4) for k in range(NK)]
USE_W = 8


# ------------------------------------------------------------- 2. F on a ball
def Fu(u, m, big):
    s = T2A if big else arb(0)
    for j in range(NK):
        if j == m:
            continue
        t = sg[j] * aA[j] / (u * u - j * j)
        s += t * k4[j] if big else t * (2 * u)
    if big:
        r = 2 * BA * u.sin_pi() / (PI * u**3) * s
        if m >= 1:
            r += 2 * BA * aA[m] * k4[m] / (u**3 * (u + m)) * (u - m).sinc_pi()
        return r
    r = BA * u.sin_pi() / PI * s
    if m == 0:
        return 2 * BA * aA[0] * u.sinc_pi() + r
    if m >= 1:
        return BA * aA[m] * (2 * u / (u + m)) * (u - m).sinc_pi() + r
    return r


def mfor(lo, hi):
    for n in (int(lo), int(lo) + 1):
        if lo <= n <= hi and 0 <= n <= KMAX:
            return n
    return -1


def ball(lo, hi):
    c, r = (lo + hi) / 2, (hi - lo) / 2
    return arb(fmpq(c.numerator, c.denominator), fmpq(r.numerator, r.denominator))


say('\n== 2. spot checks on F ==')
say('F(0) = 2*B*a_0 = %s' % (2 * BA * aA[0]).str(18))
for uu in ('8.31', '20.77', '55.13'):
    b = arb(uu)
    r1, r2 = Fu(b, -1, False), Fu(b, -1, True)
    say('  u=%-8s raw/W overlap %s' % (uu, r1.overlaps(r2)))


# ------------------------------------------------------------- 3. sign sweep
def Psi(u):
    s = arb(0)
    for k in range(NK):
        s += aA[k] * ((PI * (u - k)).si() + (PI * (u + k)).si())
    return s / (2 * PI)


say('\n== 3. adaptive sign sweep on u in [0,%d] ==' % U_EXACT)
t0 = time.time()
leaves, nev = [], 0
for cell in range(2 * U_EXACT):
    lo0, hi0 = Fr(cell, 2), Fr(cell + 1, 2)
    m = mfor(lo0, hi0)
    big = lo0 >= USE_W
    out, n = sweep(lo0, hi0, lambda l, h: Fu(ball(l, h), m, big), MINW)
    nev += n
    leaves.extend(out)
say('leaves %d, F-ball evaluations %d, %.1f s' % (len(leaves), nev, time.time() - t0))
unc = [L for L in leaves if L[2] == 0]
unc_len = float(sum(h - l for l, h, _ in unc))
say('undecided leaves: %d, total length %.3e' % (len(unc), unc_len))
sq = [s for _, _, s in leaves if s]
nsign = sum(1 for i in range(len(sq) - 1) if sq[i] != sq[i + 1])
say('sign changes: %d   (structure predicts 64-below-KMAX-ish + one per integer above)' % nsign)

# ------------------------------------------------------------- 4. int |F|
t0 = time.time()
half = arb(0)
psi_cache = {}


def PsiAt(f):
    if f not in psi_cache:
        psi_cache[f] = Psi(Q(f))
    return psi_cache[f]


i, runs = 0, 0
while i < len(leaves):
    lo, hi, s = leaves[i]
    if s == 0:
        v = Fu(ball(lo, hi), mfor(lo, hi), lo >= USE_W)
        half += Q(hi - lo) * v.abs_upper()
        i += 1
        continue
    j = i
    while j < len(leaves) and leaves[j][2] == s:
        j += 1
    half += (PsiAt(leaves[j - 1][1]) - PsiAt(lo)).abs_upper()
    runs += 1
    i = j
say('\n== 4. |F| integrated in closed form ==')
say('constant-sign runs: %d,  Psi evaluations: %d,  %.1f s'
    % (runs, len(psi_cache), time.time() - t0))
say('int_0^{x_c} |F| dx <= %s     (x_c = %.4f)'
    % (arb(half.upper()).str(20), U_EXACT / (2 * float(B))))

uc = arb(U_EXACT)
Wb = Q(abs(T2)) + sum(Q(abs(a[k]) * Fr(k**4)) / (uc**2 - k * k) for k in range(NK))
tail = Wb / (2 * PI * uc**2)
say('int_{x_c}^inf |F| dx <= %s   (analytic 1/x^3 tail)' % tail.str(6))
L1 = arb((2 * (arb(half.upper()) + tail)).upper())
say('||F||_1 <= %s' % L1.str(20))

# ------------------------------------------------------------- 5. penalty
t0 = time.time()
kb = [arb(k) / BA for k in range(NK)]
c1 = [-aA[k] * PI * k / BA for k in range(NK)]
c2 = [-aA[k] * (PI * k / BA)**2 for k in range(NK)]


def Fhat(t):
    s = arb(0)
    for k in range(NK):
        s += aA[k] * (kb[k] * t).cos_pi()
    return s


def Fhat1(t):
    s = arb(0)
    for k in range(1, NK):
        s += c1[k] * (kb[k] * t).sin_pi()
    return s


def Fhat2(t):
    s = arb(0)
    for k in range(1, NK):
        s += c2[k] * (kb[k] * t).cos_pi()
    return s


def PsiH(t):
    s = aA[0] * t
    for k in range(1, NK):
        s += aA[k] * BA / (PI * k) * (kb[k] * t).sin_pi()
    return s


def FhatT(lo, hi):
    """2nd-order Taylor (mean-value) enclosure of Fhat on [lo,hi]."""
    c, r = (lo + hi) / 2, (hi - lo) / 2
    h = r * r / 4
    return (Fhat(Q(c)) + Fhat1(Q(c)) * arb(fmpq(0), fq(r))
            + Fhat2(ball(lo, hi)) * arb(fq(h), fq(h)))


say('\n== 5. penalty ==')
# Fhat(B) = Fhat'(B) = 0 exactly, so Fhat(t) = (1/2)Fhat''(xi)(t-B)^2 for some
# xi in (t,B).  A certified constant sign of Fhat'' on [T1,B] therefore fixes
# the sign of Fhat there without any sweep of Fhat itself -- which is necessary,
# because at distance d from B, |Fhat| ~ c d^2 while any ball form has radius
# ~d, so subdivision alone never resolves the sign.
T1 = endsign = None
for delta in DELTAS:
    cand = B - delta
    if cand <= (1 + B) / 2:
        continue
    d2, n2 = sweep(cand, B, lambda l, h: Fhat2(ball(l, h)), delta / 2**20, evcap=20000)
    if d2 is None:
        continue
    if all(s == -1 for _, _, s in d2):
        T1, endsign, nd2 = cand, -1, n2
        break
    if all(s == 1 for _, _, s in d2):
        T1, endsign, nd2 = cand, +1, n2
        break

pen = arb(0)
if T1 is None:
    # Fallback.  The Taylor argument above needs Fhat'' to keep one sign next to
    # t = B, and it can fail: if the LP also drives Fhat''(B) to 0 then Fhat
    # vanishes to 4th order there and Fhat'' has its own sign changes nearby.
    # No sign argument is actually required.  The Taylor enclosure FhatT has
    # radius O(r^2) on a leaf of half-width r, so bisecting into t = B costs one
    # extra leaf per level, and the leaf that touches B is charged
    # width x sup|Fhat|, which shrinks like r^3.  At the width floor
    # (B-1)/2^40 that charge is below 1e-30.  So sweep the whole of [1, B].
    T1, endsign = B, 0
    say("Fhat'' has no certified constant sign next to t = B "
        "(Fhat vanishes to higher than 2nd order); sweeping all of [1, B] instead")
else:
    say("Fhat'' has constant sign %+d on [%s, B]  (%d evaluations); "
        "so Fhat %s 0 there by Taylor" % (endsign, T1, nd2, '<=' if endsign < 0 else '>='))
    if endsign > 0:                   # Fhat >= 0 on [T1,B]: integrate it exactly
        pen += PsiH(BA) - PsiH(Q(T1))

hl, nh = sweep(Fr(1), T1, FhatT, (T1 - 1) / 2**40)
i, nmix, nrun, pen_unc = 0, 0, 0, arb(0)
while i < len(hl):
    lo, hi, s = hl[i]
    if s == -1:
        i += 1
    elif s == 1:
        j = i
        while j < len(hl) and hl[j][2] == 1:
            j += 1
        pen += PsiH(Q(hl[j - 1][1])) - PsiH(Q(lo))
        nrun += 1
        i = j
    else:
        c = Q(hi - lo) * FhatT(lo, hi).abs_upper()
        pen += c
        pen_unc += c
        nmix += 1
        i += 1
penUB = arb(pen.upper())
say('leaves %d, Fhat evaluations %d, positive runs %d, undecided leaves %d, %.1f s'
    % (len(hl), nh, nrun, nmix, time.time() - t0))
say('   of the penalty, %.3e is charged to undecided leaves' % float(pen_unc.upper()))
say('int_1^B (Fhat)_+ dt <= %s' % penUB.str(15))

# ------------------------------------------------------------- 6. result
num = arb((2 * BA * aA[0] - 2 * AA * penUB).lower())
say('\n== 6. certified value ==')
say('F(0)              = %s' % (2 * BA * aA[0]).str(20))
say('2A*penalty       <= %s' % (2 * AA * penUB).str(15))
say('numerator        >= %s' % num.str(20))
say('||F||_1          <= %s' % L1.str(20))
C = arb((num / L1).lower())
cval = float(C.lower())
say('\nC_+(%s) >= %.16f   (certified)' % (A_PARAM, cval))
print('CERTJSON ' + json.dumps({
    'file': path, 'A': float(A_PARAM), 'B': str(B), 'nk': NK,
    'certified': repr(cval),
    'F0': repr(float((2 * BA * aA[0]).lower())),
    'penalty_ub': repr(float(penUB.upper())),
    'L1_ub': repr(float(L1.upper())),
    'u_exact': U_EXACT, 'T1': str(T1), 'endsign': endsign,
    'sign_changes': nsign, 'undecided_F': len(unc), 'undecided_F_len': unc_len,
    'undecided_pen': nmix, 'pen_unc': float(pen_unc.upper()),
}))
