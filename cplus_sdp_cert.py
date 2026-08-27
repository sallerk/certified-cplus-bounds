"""Rigorous certificate for a Gaussian-times-polynomial lower bound on C_+(A).

Every number below is an ARB *ball*: a pair [centre +/- radius] guaranteed to
contain the true value.  Ball arithmetic carries that guarantee through every
operation, so the final number is a proof, not an estimate.

WHAT IS CERTIFIED.  The test function is

    F(x) = sum_k c_k h_k(x),   h_k(x) = (2 pi)^{1/4} psi_{2k}(x sqrt(2 pi))

where psi_n is the n-th Hermite function, normalised so that int psi_n^2 = 1,
and the c_k are the float64 numbers in the input file, read as exact binary
fractions.  F is a polynomial times a Gaussian, so it is even, real, continuous
and in L^1 for ANY coefficients whatever.  Unlike the bandlimited certificate
cplus_certA.py, there is NO admissibility side condition and NO projection step:
the file is the function.  We certify

    C_+(A)  >=  [ F(0) - A int_{|t|>c} (Fhat(t))_+ dt ] / ( c ||F||_1 )

which is a valid lower bound for C_+(A) for every c > 0.  c = 1 is the problem
exactly as Carneiro-Milinovich-Soundararajan state it; c != 1 is the same bound
applied to the dilated function F(x/c).

FOUR EXACT FACTS DO THE WORK.

1.  Fhat(t) = sum_k (-1)^k c_k h_k(t).  Even Hermite functions are
    eigenfunctions of the Fourier transform, eigenvalue (-1)^k.  The transform
    is a sign flip: no error at all.

2.  int_a^b psi_n has a CLOSED FORM.  From psi_n = sqrt((n-1)/n) psi_{n-2}
    - sqrt(2/n) psi_{n-1}',

      J_n(a,b) = sqrt((n-1)/n) J_{n-2}(a,b) - sqrt(2/n)[psi_{n-1}(b)-psi_{n-1}(a)]

    starting from J_0 = pi^{-1/4} sqrt(pi/2)[erf(b/sqrt2) - erf(a/sqrt2)].  So
    there is NO quadrature error anywhere.  All that is left is to decide the
    sign of the integrand, because int|F| = |int F| where the sign is constant.

3.  A GLOBAL bound  |psi_n(y)| <= sqrt(2) (n+1/2)^{1/4}, proved in two lines:
    psi_n(y)^2 = 2 int_{-inf}^y psi psi' <= 2 ||psi_n||_2 ||psi_n'||_2, and
    ||psi_n||_2 = 1, ||psi_n'||_2^2 = -int psi psi'' = (2n+1) - int y^2 psi^2
    = (2n+1) - (n+1/2) = n+1/2.

4.  A DECAYING bound  |psi_n(y)| <= r_n(y) e^{-y^2/2}, where r_n obeys the same
    three-term recurrence with a PLUS sign (induction + triangle inequality) and
    has non-negative coefficients, hence is increasing on y >= 0.

WHY A TAYLOR MODEL AND NOT A PLAIN BALL EVALUATION.  Feeding an interval into
the Hermite recurrence is useless: the recurrence does not track the dependence
on y, so an input of radius 1e-11 comes out with radius 1.8e5, an amplification
of 1.8e16.  That is measured, not assumed.  Worse, in the range where F is
exponentially small the bounds of facts 3 and 4 are both far above |F|, because
F's decay there IS a cancellation between individually O(1) Hermite functions.
So F and its first J-1 derivatives are evaluated at the exact MIDPOINT of each
interval (radius ~1e-90, no amplification), and only the J-th derivative is
bounded crudely -- multiplied by (w/2)^J / J!, which annihilates the crudeness.

Usage:  python cplus_sdp_cert.py <coefficient file> <A> [c]
Env:    PREC (bits, default 300), JORD (Taylor order, default 6),
        NSEED (seed intervals, default 4000), FLOORBITS (default 40),
        XCUT (force the sweep cutoff), TAILTGT, QUIET.
"""

import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

PREC = int(os.environ.get('PREC', 300))
ctx.prec = PREC
JORD = int(os.environ.get('JORD', 6))
NSEED = int(os.environ.get('NSEED', 4000))
FLOOR = 2.0 ** -int(os.environ.get('FLOORBITS', 40))
XCUT = float(os.environ.get('XCUT', 0))
TAILTGT = float(os.environ.get('TAILTGT', 1e-20))


class Herm:
    def __init__(self, coef, J=JORD):
        self.c = [arb(float(v)) for v in coef]
        self.ac = [abs(v) for v in self.c]
        self.nk = len(coef) - 1
        self.N = 2 * self.nk
        self.J = J
        self.M = self.N + J
        self.pi = arb.pi()
        self.rt2 = arb(2).sqrt()
        self.rt2pi = (2 * self.pi).sqrt()
        self.pim14 = self.pi ** arb(-0.25)
        self.sc = (2 * self.pi) ** arb(0.25)
        M = self.M
        self.a = [None, None] + [(arb(2) / n).sqrt() for n in range(2, M + 3)]
        self.b = [None, None] + [(arb(n - 1) / n).sqrt() for n in range(2, M + 3)]
        self.Sglob = [self.rt2 * (arb(n) + arb(0.5)) ** arb(0.25) for n in range(M + 1)]
        self.hf = [(arb(n) / 2).sqrt() for n in range(M + 3)]
        self.fact = [arb(math.factorial(j)) for j in range(J + 2)]
        self.ders = {False: self._ders(False), True: self._ders(True)}

    # ---------------------------------------------------------------- Hermite
    def psis(self, y, upto=None):
        upto = self.M if upto is None else upto
        g = (-y * y / 2).exp() * self.pim14
        out = [g, self.rt2 * y * g]
        for n in range(1, upto):
            out.append(self.a[n + 1] * y * out[n] - self.b[n + 1] * out[n - 1])
        return out[:upto + 1]

    def rpoly(self, y, upto=None):
        upto = self.M if upto is None else upto
        out = [self.pim14, self.rt2 * y * self.pim14]
        for n in range(1, upto):
            out.append(self.a[n + 1] * y * out[n] + self.b[n + 1] * out[n - 1])
        return out[:upto + 1]

    def sup_psi(self, xa, xb):
        A, Bb = xa * self.rt2pi, xb * self.rt2pi
        r = self.rpoly(Bb)
        dec = (-A * A / 2).exp()
        out = []
        for n in range(self.M + 1):
            loc = r[n] * dec
            out.append(loc if float(loc) < float(self.Sglob[n]) else self.Sglob[n])
        return out

    # ------------------------------------------------- derivative coefficients
    def _ders(self, flip):
        # h_k(x) = (-1)^k (2 pi)^{1/4} psi_{2k}(x sqrt(2 pi)):  both sides are
        # L^2-normalised, and h_k(0) > 0 while psi_{2k}(0) has sign (-1)^k
        # because H_{2k}(0) = (-1)^k (2k)!/k!.  Fhat contributes a second (-1)^k.
        a = [arb(0)] * (self.M + 1)
        for k in range(self.nk + 1):
            v = self.c[k] * self.sc
            if (k % 2) and not flip:
                v = -v
            a[2 * k] = v
        out = [a]
        for _ in range(self.J):
            p = out[-1]
            q = [arb(0)] * (self.M + 1)
            for m in range(self.M + 1):
                s = arb(0)
                if m + 1 <= self.M:
                    s += p[m + 1] * self.hf[m + 1]
                if m >= 1:
                    s -= p[m - 1] * self.hf[m]
                q[m] = s
            out.append(q)
        return out

    def F(self, x, flip=False):
        p = self.psis(x * self.rt2pi, self.N)
        d = self.ders[flip][0]
        s = arb(0)
        for m in range(0, self.N + 1, 2):
            s += d[m] * p[m]
        return s

    def encl(self, xc, h, flip, S):
        """Ball containing F (or Fhat) on [xc-h, xc+h]."""
        p = self.psis(xc * self.rt2pi)
        D = self.ders[flip]
        hb = arb(0, h)
        tot = arb(0)
        pw = arb(1)
        sy = arb(1)
        for j in range(self.J):
            v = arb(0)
            Dj = D[j]
            for m in range(self.M + 1):
                v += Dj[m] * p[m]
            tot += v * sy * pw / self.fact[j]
            pw = pw * hb
            sy = sy * self.rt2pi
        rem = arb(0)
        DJ = D[self.J]
        for m in range(self.M + 1):
            rem += abs(DJ[m]) * S[m]
        tot += arb(0, float((rem * sy * arb(h) ** self.J / self.fact[self.J]).upper()))
        return tot

    # ------------------------------------------------------- exact integration
    def integ(self, xa, xb, flip=False):
        ya, yb = xa * self.rt2pi, xb * self.rt2pi
        pa, pb = self.psis(ya, self.N), self.psis(yb, self.N)
        j = [None] * (self.N + 1)
        j[0] = self.pim14 * (self.pi / 2).sqrt() * ((yb / self.rt2).erf() - (ya / self.rt2).erf())
        if self.N >= 1:
            j[1] = -self.pim14 * self.rt2 * ((-yb * yb / 2).exp() - (-ya * ya / 2).exp())
        for n in range(2, self.N + 1):
            j[n] = self.b[n] * j[n - 2] - self.a[n] * (pb[n - 1] - pa[n - 1])
        d = self.ders[flip][0]
        s = arb(0)
        for m in range(0, self.N + 1, 2):
            s += d[m] * j[m]
        return s / self.rt2pi

    def tailbound(self, X):
        Y = X * self.rt2pi
        rr = self.rpoly(Y, self.N)
        RY = arb(0)
        for k in range(self.nk + 1):
            RY += self.ac[k] * rr[2 * k]
        RY *= self.sc
        N = self.N
        g = arb(2) ** (arb(N - 1) / 2) * (Y * Y / 2).gamma_upper(arb(N + 1) / 2)
        return RY * Y ** (-N) * g / self.rt2pi


def sweep(H, lo, hi, flip, positive_part):
    """Adaptive sign sweep.  Definite-sign runs are integrated in closed form.
    Undecided pieces are charged width * (bound on |F| there), an over-estimate.
    Nothing claims that every sign change was found: a missed one shows up as a
    piece that never settles and is charged at the over-estimate rate."""
    leaves = []
    for i in range(NSEED):
        sa = lo + (hi - lo) * i / NSEED
        sb = lo + (hi - lo) * (i + 1) / NSEED
        S = H.sup_psi(arb(sa), arb(sb))
        stack = [(arb(sa), arb(sb))]
        while stack:
            a, b = stack.pop()
            c = (a + b) / 2
            v = H.encl(c, float((b - a) / 2), flip, S)
            if v > 0:
                leaves.append((a, b, 1, None))
            elif v < 0:
                leaves.append((a, b, -1, None))
            elif float(b - a) < FLOOR:
                leaves.append((a, b, 0, (b - a) * arb(abs(v).upper())))
            else:
                stack.append((c, b))
                stack.append((a, c))
    leaves.sort(key=lambda z: float(z[0]))
    total = arb(0); und_len = arb(0); und_chg = arb(0)
    nsign = 0; i = 0; prev = None
    while i < len(leaves):
        a, b, s, chg = leaves[i]
        if s == 0:
            und_len += b - a
            und_chg += chg
            total += chg
            i += 1
            continue
        j = i
        while j + 1 < len(leaves) and leaves[j + 1][2] == s:
            j += 1
        val = H.integ(a, leaves[j][1], flip)
        total += abs(val) if not positive_part else (val if s > 0 else arb(0))
        if prev is not None and prev != s:
            nsign += 1
        prev = s
        i = j + 1
    return total, len(leaves), nsign, und_len, und_chg


def main(path, A, c=1.0):
    t0 = time.time()
    coef = np.loadtxt(path)
    H = Herm(coef)
    if XCUT:
        X = XCUT
    else:
        X = 2.0
        while float(H.tailbound(arb(X))) > TAILTGT and X < 40.0:
            X += 0.25
    Ab = arb(float(A)); cb = arb(float(c))
    F0 = H.F(arb(0))
    tail = H.tailbound(arb(X))
    pen, npc, nsp, ulp, ucp = sweep(H, c, X, True, True)
    pen = 2 * (pen + tail)
    l1, nlc, nsl, ull, ucl = sweep(H, 0.0, X, False, False)
    l1 = 2 * (l1 + tail)
    num = F0 - Ab * pen
    val = num / (cb * l1)
    lo = float(arb(val.lower()))
    if not os.environ.get('QUIET'):
        print(f"file            {path}")
        print(f"A = {A}   c = {c}   Hermite orders 0..{H.N}   Taylor order {H.J}"
              f"   prec {PREC} bits")
        print(f"cutoff X        {X:.4f}   tail bound <= {float(tail):.3e}")
        print(f"|F|  sweep      {nlc} pieces, {nsl} sign changes, undecided length "
              f"{float(ull):.2e}, charge {float(ucl):.2e}")
        print(f"Fhat sweep      {npc} pieces, {nsp} sign changes, undecided length "
              f"{float(ulp):.2e}, charge {float(ucp):.2e}")
        print(f"F(0)            {F0.str(20)}")
        print(f"penalty  <=     {pen.str(20)}")
        print(f"||F||_1  <=     {l1.str(20)}")
        print(f"numerator >=    {num.str(20)}")
        print(f"CERTIFIED >=    {lo:.16f}      ({time.time()-t0:.1f} s)")
    return lo


if __name__ == '__main__':
    main(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 1.0)
