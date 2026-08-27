"""Rigorous (ARB) UPPER bound on C_+(A) from a dual multiplier.

WHAT IS PROVED.  Given a step function m with

    m(t) = 1            for |t| <= 1
    1 - A <= m(t) <= 1  for |t| >  1,   m(t) = 0 for |t| > T

the CMS duality (arXiv:1708.04122 sec 4.5; derivation re-checked in UPPER.md)
gives  C_+(A) <= sup_x |mcheck(x)|,  mcheck(x) = int e^{2 pi i x t} m(t) dt.
This script certifies an upper bound on that supremum, so it certifies an
upper bound on C_+(A).  Direction of rounding is the OPPOSITE of the lower-bound
certificates in cplus_certA.py: every quantity is rounded UP.

A max over a finite grid is not a max.  The sup is covered in three pieces:

  [0, delta]   |mcheck(x) - mcheck(0)| <= 4 pi^2 x^2 M2,  M2 = int_0^T |m| t^2 dt.
               (from cos(u) = 1 - 2 sin^2(u/2) and |sin(pi x t)| <= pi x t.)
  [delta, X0]  adaptive bisection with a SECOND-ORDER Taylor test:
               on a cell of centre c and half-width h,
                 |mcheck(x)| <= |mcheck(c)| + |mcheck'(c)| h + (1/2) Mb2 h^2
               with Mb2 = 8 pi^2 M2 >= sup |mcheck''|.  mcheck(c) and mcheck'(c)
               are ARB enclosures, not floats.
  [X0, inf)    |mcheck(x)| <= V/(pi x) with V = total variation of m
               (telescoped form mcheck = (1/(pi x)) sum_j c_j sin(2 pi x tau_j),
               |sin| <= 1).  X0 is chosen so this is already <= the target.

Exact evaluation forms (both cancellation-free where they are used):
  mcheck(x)  = 2 sum_j m_j d_j cos(pi x (tau_j + tau_{j-1})) sinc(x d_j),
               d_j = tau_j - tau_{j-1},  sinc(u) = sin(pi u)/(pi u)
  mcheck'(x) = -(4 pi / w^2) sum_j m_j [ g(w tau_j) - g(w tau_{j-1}) ],
               w = 2 pi x,  g(z) = sin z - z cos z

Usage:  python cplus_dual_cert.py multiplier_A28.txt [target]
"""
import sys, os, time
from fractions import Fraction as Fr
from flint import arb, ctx, fmpq

ctx.dps = 60
PI = arb.pi()
QUIET = os.environ.get('QUIET') == '1'
MINH = Fr(1, 2 ** 34)


def say(*a, **k):
    if not QUIET:
        print(*a, **k)


def Q(f):
    return arb(fmpq(f.numerator, f.denominator))


def ball(lo, hi):
    c, r = (lo + hi) / 2, (hi - lo) / 2
    return arb(fmpq(c.numerator, c.denominator), fmpq(r.numerator, r.denominator))


# --------------------------------------------------------------------- input
def load(path):
    A = T = None
    tau, mval = [], []
    for line in open(path):
        s = line.strip()
        if s.startswith('A ='):
            A = Fr(s.split('=')[1].strip())
        elif s.startswith('T ='):
            T = Fr(s.split('=')[1].strip())
        elif s.startswith('cell'):
            _, lo, hi, v = s.split()
            if not tau:
                tau.append(Fr(lo))
            assert Fr(lo) == tau[-1], 'cells must be contiguous'
            tau.append(Fr(hi))
            mval.append(Fr(v))
    return A, T, tau, mval


path = sys.argv[1] if len(sys.argv) > 1 else 'multiplier_A28.txt'
A, T, tau, mv = load(path)
N = len(mv)

say('== 1. the multiplier, as exact rationals ==')
say('file %s :  A = %s,  T = %s,  %d cells on [1,T]' % (path, A, T, N))
assert tau[0] == 1 and tau[-1] == T
assert all(tau[i] < tau[i + 1] for i in range(N)), 'breakpoints must increase'
lo_clamp = 1 - A
bad = [i for i, v in enumerate(mv) if not (lo_clamp <= v <= 1)]
assert not bad, 'admissibility violated in cells %s' % bad[:5]
say('admissibility  %s <= m_j <= 1  : VERIFIED exactly in rational arithmetic'
    % lo_clamp)
say('min m_j = %s   max m_j = %s' % (min(mv), max(mv)))

# c_j = m_j - m_{j+1} at tau_j, with m_0 = 1 (on [0,1]) and m_{N+1} = 0
mfull = [Fr(1)] + list(mv) + [Fr(0)]
c = [mfull[j] - mfull[j + 1] for j in range(N + 1)]
V = sum(abs(x) for x in c)
d = [tau[j] - tau[j - 1] for j in range(1, N + 1)]

# mcheck(0) = 2 int_0^T m  (two independent exact forms, cross-checked)
m0_cells = 2 * (Fr(1) + sum(mv[j] * d[j] for j in range(N)))
m0_tele = 2 * sum(c[j] * tau[j] for j in range(N + 1))
assert m0_cells == m0_tele, 'the two exact forms of mcheck(0) disagree'
M0 = m0_cells
# M1 = int_0^T |m| t dt,  M2 = int_0^T |m| t^2 dt   (exact)
M1 = Fr(1, 2) + sum(abs(mv[j]) * (tau[j + 1] ** 2 - tau[j] ** 2) for j in range(N)) / 2
M2 = Fr(1, 3) + sum(abs(mv[j]) * (tau[j + 1] ** 3 - tau[j] ** 3) for j in range(N)) / 3
say('mcheck(0) = 2 int m  = %.15f   (exact rational, both forms agree)' % float(M0))
say('V = total variation  = %.6f' % float(V))
say('M1 = int |m| t dt    = %.6f' % float(M1))
say('M2 = int |m| t^2 dt  = %.6f' % float(M2))

TARGET = Fr(sys.argv[2]) if len(sys.argv) > 2 else None

# ------------------------------------------------------------- ARB constants
tauA = [Q(t) for t in tau]
mA = [Q(v) for v in mv]
dA = [Q(x) for x in d]
sumA = [Q(tau[j] + tau[j - 1]) for j in range(1, N + 1)]
MB2 = 8 * PI ** 2 * Q(M2)          # >= sup |mcheck''|
MB1 = 4 * PI * Q(M1)               # >= sup |mcheck'|


def mcheck(x):
    """ARB enclosure of mcheck(x).  Safe at x = 0 (sinc)."""
    s = 2 * (PI * x).cos() * x.sinc_pi()                # cell [0,1]: m=1, d=1
    for j in range(N):
        s += 2 * mA[j] * dA[j] * (PI * x * sumA[j]).cos() * (x * dA[j]).sinc_pi()
    return s


def gz(z):
    return z.sin() - z * z.cos()


def dmcheck(x):
    """ARB enclosure of mcheck'(x).  Used only for x >= delta > 0."""
    w = 2 * PI * x
    s = gz(w * tauA[0]) - gz(w * arb(0))               # cell [0,1], m = 1
    for j in range(N):
        s += mA[j] * (gz(w * tauA[j + 1]) - gz(w * tauA[j]))
    return -4 * PI * s / (w * w)


say('\n== 2. spot checks ==')
say('mcheck(0) by ARB sinc form = %s' % mcheck(arb(0)).str(18))
say('mcheck(0) exact rational   = %s' % Q(M0).str(18))
assert mcheck(arb(0)).overlaps(Q(M0))
h = arb(10) ** -20
say("mcheck'(0.7) closed form   = %s" % dmcheck(arb('0.7')).str(15))
say("mcheck'(0.7) finite diff   = %s"
    % ((mcheck(arb('0.7') + arb('1e-12')) - mcheck(arb('0.7') - arb('1e-12')))
       / arb('2e-12')).str(15))


# ---------------------------------------------------------------- 3. the sup
def certify(target):
    """Try to prove sup_x |mcheck(x)| <= target.  Returns (ok, stats)."""
    tg = Q(target)
    # --- piece 1: [0, delta]
    room = tg - abs(Q(M0))
    if room.lower() <= 0:
        return False, 'target is below |mcheck(0)| itself'
    delta = (room / (4 * PI ** 2 * Q(M2))).sqrt()
    dlo = Fr(int(float(arb(delta.lower())) * 2 ** 30), 2 ** 30)   # rounded DOWN
    assert dlo > 0
    # --- piece 3: [X0, inf)
    X0 = Fr(V) / Fr(3141592653589793, 10 ** 15) / target   # V/(pi*target), rounded UP
    X0 = Fr(int(X0 * 2 ** 20) + 1, 2 ** 20)
    # --- piece 2: [delta, X0] by adaptive second-order bisection
    t0 = time.time()
    stack, ncell, nev, minh = [(dlo, X0)], 0, 0, Fr(1)
    while stack:
        lo, hi = stack.pop()
        cc, hh = (lo + hi) / 2, (hi - lo) / 2
        val = abs(mcheck(Q(cc))) + abs(dmcheck(Q(cc))) * Q(hh) \
            + MB2 * Q(hh) ** 2 / 2
        nev += 1
        if val.upper() <= tg.lower():
            ncell += 1
            minh = min(minh, hh)
            continue
        if hh <= MINH:
            return False, ('undecided cell at x ~ %.9f, width %.3e'
                           % (float(cc), float(2 * hh)))
        stack.append((cc, hi))
        stack.append((lo, cc))
    return True, dict(delta=float(dlo), X0=float(X0), cells=ncell, evals=nev,
                      minh=float(minh), secs=time.time() - t0)


if TARGET is None:
    say('\nno target given; nothing certified')
    sys.exit(0)

say('\n== 3. certifying sup_x |mcheck(x)| <= %s ==' % TARGET)
ok, info = certify(TARGET)
if ok:
    say('[0, delta]   delta = %.6e   |mcheck| <= |mcheck(0)| + 4 pi^2 delta^2 M2 <= target'
        % info['delta'])
    say('[delta, X0]  X0 = %.4f, %d certified cells, %d ARB evaluations, '
        'finest half-width %.3e, %.1f s'
        % (info['X0'], info['cells'], info['evals'], info['minh'], info['secs']))
    say('[X0, inf)    |mcheck| <= V/(pi x) <= target')
    print('\nCERTIFIED   C_+(%s) <= %s' % (A, TARGET))
else:
    print('\nNOT certified at %s : %s' % (TARGET, info))
    sys.exit(2)
