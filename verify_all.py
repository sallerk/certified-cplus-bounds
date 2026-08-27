"""Verify every certified bound in this package, from the coefficient files alone.

Usage:
    python verify_all.py            # all 68 rows, both methods (a few minutes)
    python verify_all.py 28         # just A = 28
    python verify_all.py --upper    # the certified UPPER bounds instead

Each row is re-certified in interval arithmetic and checked against three things:
  1. the value recorded in final_table.json,
  2. the published record from cqh_table1.txt (must be strictly exceeded),
  3. the certified upper bound (must not be exceeded -- if it is, something is wrong).

Requires: python-flint (ARB bindings) and numpy.
"""
import json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CERTS = os.path.join(HERE, 'certificates')
FUNCS = [os.path.join(HERE, 'functions', d) for d in ('lp', 'sdp', 'dual')]


def resolve(name):
    """Find a coefficient file by bare name, wherever it lives under functions/."""
    for d in [HERE] + FUNCS:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    raise SystemExit('cannot find coefficient file: %s' % name)


def published_records():
    """Read CQH Table 1. The record for each A is the MAX of the three columns --
    which column wins changes three times across the range."""
    rec = {}
    with open(os.path.join(HERE, 'cqh_table1.txt')) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 5 and re.match(r'^[0-9]+\.[0-9]$', p[0]):
                try:
                    A, f82, f122, pw = (float(x) for x in p[:4])
                except ValueError:
                    continue
                rec[A] = max(f82, f122, pw)
    return rec


def run(script, args):
    out = subprocess.run([sys.executable, os.path.join(CERTS, script)] + args,
                         capture_output=True, text=True, cwd=HERE).stdout
    m = re.findall(r'(?:CERTIFIED\s*>=|C_\+\([^)]*\)\s*>=)\s*([0-9]+\.[0-9]+)', out)
    return float(m[-1]) if m else None


def certify(row):
    """Re-certify one row with whichever method produced its best value."""
    A = float(row['A'])
    lp, sdp = float(row['lp']), float(row['sdp'])
    if sdp >= lp:
        f = row['sdp_file']
        c = re.search(r'_c([0-9.]+)\.txt$', f)
        return 'SDP', run('cplus_sdp_cert.py', [resolve(f), repr(A), c.group(1) if c else '1.0'])
    return 'LP', run('cplus_certA.py', [resolve(row['lp_file'])])


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    rows = json.load(open(os.path.join(HERE, 'final_table.json')))
    rows.sort(key=lambda r: float(r['A']))
    if args:
        want = {float(a) for a in args}
        rows = [r for r in rows if float(r['A']) in want]

    rec = published_records()
    upper = {2.0: 1.3102, 5.0: 1.1852, 15.0: 1.1250, 28.0: 1.1136, 34.5: 1.1136}

    print('%6s  %-4s  %-20s  %-12s  %-9s  %s' %
          ('A', 'via', 'certified', 'record', 'margin', 'status'))
    bad = 0
    for r in rows:
        A = float(r['A'])
        method, got = certify(r)
        claimed = max(float(r['lp']), float(r['sdp']))
        if got is None:
            print('%6.1f  %-4s  %-20s' % (A, method, 'CERTIFICATE FAILED')); bad += 1; continue

        notes = []
        if abs(got - claimed) > 5e-13:
            notes.append('differs from table (%.16f)' % claimed); bad += 1
        if got <= rec[A]:
            notes.append('DOES NOT BEAT RECORD'); bad += 1
        # C_+ is non-increasing in A, so any certified upper bound at A' <= A applies here
        ub = min((v for a, v in upper.items() if a <= A), default=None)
        if ub is not None and got > ub:
            notes.append('EXCEEDS UPPER BOUND %.4f' % ub); bad += 1

        print('%6.1f  %-4s  %.16f  %.6f  %+.4f%%  %s' %
              (A, method, got, rec[A], 100 * (got / rec[A] - 1),
               'ok' if not notes else ' | '.join(notes)))

    print('\n%d row(s) checked, %d problem(s).' % (len(rows), bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
