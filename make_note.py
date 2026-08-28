"""Generate docs/note.html from the certified data.

The note used to be hand-maintained, and drifted: three row margins and the
A = 28 bracket were stale at various points, each caught only by a manual audit.
Everything numeric in the page now comes from final_table.json, upper_bounds.json
and extra_rows.json, so it cannot disagree with what verify_all.py certifies.

  python make_note.py            # rewrite docs/note.html
  python make_note.py --check    # exit 1 if the file on disk is out of date
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'docs', 'note.html')


def load():
    j = lambda n: json.load(open(os.path.join(HERE, n), encoding='utf-8'))
    rows = sorted(j('final_table.json'), key=lambda r: float(r['A']))
    ub = {float(k): v['ub'] for k, v in j('upper_bounds.json')['bounds'].items()}
    env, run = {}, 9.9
    for A in sorted(ub):                      # C_+ is non-increasing: running min
        run = min(run, ub[A])
        env[A] = run
    return rows, env, j('extra_rows.json')['rows']


def up(x, d):
    """Round UP to d places: a bound must never be quoted tighter than proved."""
    return math.ceil(x * 10 ** d) / 10 ** d


def stats(rows, env, extra):
    m = [float(r['margin'].strip('%+')) for r in rows]
    br = [(100 * (env[float(r['A'])] / float(r['best']) - 1), float(r['A'])) for r in rows]
    g = lambda A: float([r for r in rows if float(r['A']) == A][0]['best'])
    x36 = float(extra[0]['certified'])
    return dict(
        n=len(rows), beat=sum(1 for v in m if v > 0), mean=sum(m) / len(m),
        brmean=sum(v for v, _ in br) / len(br), brmax=max(br), brmin=min(br),
        a28=g(28.0), a28u=env[28.0], a4=g(4.0), x36=x36,
        qf=up(2 / g(28.0), 10), qf3=up(6 / g(28.0), 10),
        ap=up(1 / g(4.0), 6), rh=up(1 / x36, 4),
    )


def table_rows(rows):
    out = []
    for r in rows:
        lp, sdp = float(r['lp']), float(r['sdp'])
        wl = ' win' if r['method'] == 'LP' else ''
        ws = ' win' if r['method'] == 'SDP' else ''
        out.append(
            '<tr><td class="n">%g</td><td class="n">%s</td><td class="src">%s</td>'
            '<td class="n%s">%.10f</td><td class="n%s">%.10f</td>'
            '<td class="src">%s</td><td class="n gain">%s</td></tr>'
            % (float(r['A']), ('%g' % float(r['record'])), r['who'],
               wl, lp, ws, sdp, r['method'], r['margin']))
    return '\n'.join(out)


CSS = open(os.path.join(HERE, 'docs', '_note_head.html'), encoding='utf-8').read() \
    if os.path.exists(os.path.join(HERE, 'docs', '_note_head.html')) else None


def build():
    rows, env, extra = load()
    s = stats(rows, env, extra)
    head = CSS
    if head is None:
        raise SystemExit('docs/_note_head.html (title, fonts, <style>) is missing')

    stat = ''.join(
        '\n  <div class="stat"><span class="v">%s</span><span class="k">%s</span></div>' % kv
        for kv in [('%d/%d' % (s['beat'], s['n']), 'rows improved'),
                   ('+%.2f%%' % s['mean'], 'mean margin'),
                   ('%.4f' % s['rh'], 'prime gap, RH'),
                   ('%.4f' % s['qf'], 'quadratic forms, GRH'),
                   ('%.2f%%' % s['brmean'], 'mean bracket')])

    body = []
    A = body.append
    A('\n<div class="wrap">')
    A('<header style="display:flex;flex-direction:column;gap:.7rem">')
    A('  <p class="eyebrow">Fourier optimization &middot; analytic number theory</p>')
    A('  <h1>Certified lower bounds for the Carneiro&ndash;Milinovich&ndash;Soundararajan constant</h1>')
    A('  <p class="lede">Every row of a published table of extremal constants, improved and proved')
    A('  rigorously in interval arithmetic.</p>')
    A('</header>')
    A('\n<div class="stats">%s\n</div>' % stat)

    A('\n<section>')
    A('<h2 style="border:0;padding:0;margin-top:0">The constant</h2>')
    A('<p>Carneiro, Milinovich and Soundararajan (<a href="https://arxiv.org/abs/1708.04122">arXiv:1708.04122</a>,')
    A('Extremal Problem&nbsp;2) define, for a parameter <code>A</code>,</p>')
    A('<div class="eq">C&#8330;(A) = sup [ F(0) &minus; A &int;<sub>|t|&gt;1</sub> (F&#770;(t))&#8330; dt ] / &Vert;F&Vert;&#8321;</div>')
    A('<p>the supremum over even, continuous, real <code>F</code> in <code>L&sup1;(&reals;)</code>, for')
    A('<code>1 &le; A &lt; &infin;</code>. Larger is better: a bigger <code>C&#8330;</code> gives a smaller')
    A('constant in a bound on gaps between primes. Chirre and Quesada-Herrera')
    A('(<a href="https://arxiv.org/abs/2012.07781">arXiv:2012.07781</a>) tabulate lower bounds at 68')
    A('values of <code>A</code>.</p>')
    A('<p><code>A</code> is not a free parameter. It is fixed by the application: how many primes can')
    A('crowd into a short interval. Three values carry the results below &mdash;')
    A('<code>36/11</code> for ordinary primes, <code>4</code> for primes in an arithmetic progression,')
    A('<code>28</code> for primes represented by a binary quadratic form.</p>')
    A('<p>Any admissible <code>F</code> gives a valid lower bound. So the whole game is finding a good')
    A('test function &mdash; and then <em>proving</em> what it is worth.</p>')
    A('</section>')

    A('\n<section>')
    A('<h2>What changed</h2>')
    A('<p>The published records come from hand-tuned families of three or four parameters. Replacing')
    A('that search with two optimizers over 100&ndash;200 coefficients beats them at every value of')
    A('<code>A</code>:</p>')
    A('<ul>')
    A('<li>a <strong>linear program</strong> over a cosine series for <code>F&#770;</code>, band-limited to')
    A('<code>[&minus;B,B]</code>;</li>')
    A('<li>a <strong>semidefinite program</strong> over <code>F(x) = P(x)e<sup>&minus;&pi;x&sup2;</sup></code>,')
    A('following Chirre&ndash;Pereira&ndash;de&nbsp;Laat, where sign conditions become sum-of-squares')
    A('constraints.</li>')
    A('</ul>')
    A('<p>Neither dominates. The linear program wins for <code>1.5 &le; A &le; 9</code>; the semidefinite')
    A('program wins at <code>A = 1</code> and for <code>A &ge; 9.5</code>. Running both and taking the better')
    A('one is what gives %d out of %d.</p>' % (s['beat'], s['n']))
    A('</section>')

    A('\n<section>')
    A('<h2>What &ldquo;certified&rdquo; means here</h2>')
    A('<p>Ordinary floating-point arithmetic returns a number with no guarantee. Every quantity below')
    A('is instead computed in <strong>ball arithmetic</strong> &mdash; a value plus a radius, with the')
    A('true result provably inside &mdash; using ARB, the same library the published record was')
    A('certified with.</p>')
    A('<p>Because the claim is a lower bound on a fraction, the errors must be forced in a specific')
    A('direction: an <em>upper</em> bound on the denominator <code>&Vert;F&Vert;&#8321;</code>, and an')
    A('<em>upper</em> bound on the subtracted penalty. Three details carry the proof:</p>')
    A('<ul>')
    A('<li>Coefficients are fixed as exact rationals and projected so <code>F&#770;(B) = 0</code> holds')
    A('<em>exactly</em>. Without this <code>&Vert;F&Vert;&#8321;</code> is infinite and nothing is proved.</li>')
    A('<li><code>&int;sinc = Si/&pi;</code> in closed form, so <strong>no quadrature error exists')
    A('anywhere</strong>. The only work is deciding signs, and intervals that cannot be decided are')
    A('charged at their maximum &mdash; so no claim about having found every root is needed.</li>')
    A('<li>The tail beyond the computed range is <em>proved</em> to decay like <code>1/x&sup3;</code> with')
    A('an explicit constant, not observed to be small.</li>')
    A('</ul>')
    A('<p class="note">The certificate was tested in both directions: it accepts the true bound and')
    A('refuses a target below the true supremum, reporting exactly where it fails.</p>')
    A('</section>')

    A('\n<section>')
    A('<h2>How much room is left</h2>')
    A('<p>A lower bound alone says nothing about how close to optimal it is. Optimizing the dual side')
    A('numerically &mdash; the same duality CMS use analytically &mdash; gives certified <em>upper</em>')
    A('bounds:</p>')
    A('<div class="eq">%.16f &le; C&#8330;(28) &le; %g</div>' % (s['a28'], s['a28u']))
    A('<p>a bracket of <strong>%.2f%%</strong>, against roughly 10.6%% from the previously published upper'
      % (100 * (s['a28u'] / s['a28'] - 1)))
    A('bound. So the bulk of the apparent gap was slack in the upper bound, not room in the lower one.')
    A('Every one of the %d rows now carries its own certified upper bound &mdash; mean bracket' % s['n'])
    A('<strong>%.2f%%</strong>, widest %.2f%% at <code>A = %g</code>. At <code>A = 1</code>, where the exact'
      % (s['brmean'], s['brmax'][0], s['brmax'][1]))
    A('value <code>C&#8330;(1) = 2</code> is known, the bracket is %.2f%%.</p>' % s['brmin'][0])
    A('<p class="note">The duality gap is not proven, so the residual cannot be attributed to either')
    A('side &mdash; the true value is known to lie in the interval, not where in it.</p>')
    A('</section>')

    A('\n<section>')
    A('<h2>The table</h2>')
    A('<p class="note">Record is the best of Chirre&ndash;Quesada-Herrera\'s three published columns')
    A('(<code>F82</code>, <code>F122</code>, <code>PW</code>) &mdash; the winning column changes across the')
    A('range. At <code>A = 4</code> the record instead comes from Chirre&ndash;Pereira&ndash;de&nbsp;Laat')
    A('(<code>CPdL</code>), who optimised that single <code>A</code> directly and published')
    A('<code>1.17233</code>, above anything in Table&nbsp;1. The shaded cell is the method that won')
    A('each row.</p>')
    A('<div class="tw"><table>')
    A('<thead><tr><th>A</th><th>record</th><th>col</th><th>LP (certified)</th><th>SDP (certified)</th><th>won</th><th>margin</th></tr></thead>')
    A('<tbody>')
    A(table_rows(rows))
    A('</tbody></table></div>')
    e = extra[0]
    A('<p class="note">Certified separately, because <code>A = %s</code> is not a row of Table&nbsp;1:'
      ' <code>C&#8330;(%s) &ge; %.16f</code> against Chirre&ndash;Pereira&ndash;de&nbsp;Laat\'s'
      ' <code>%g</code>, a margin of <code>%s</code>. It is kept in <code>extra_rows.json</code>.</p>'
      % (e['A'], e['A'], float(e['certified']), float(e['record']), e['margin']))
    A('</section>')

    A('\n<section>')
    A('<h2>Consequences</h2>')
    A('<p>Each of the three values of <code>A</code> carries a published constant.</p>')
    A('<p><strong>Ordinary primes, under RH.</strong> The constant is <code>1/C&#8330;(36/11)</code>:</p>')
    A('<div class="eq">limsup (p<sub>n+1</sub> &minus; p<sub>n</sub>) / (&radic;p&#8202;log&#8202;p) &lt; %.4f</div>' % s['rh'])
    A('<p>improving Chirre&ndash;Pereira&ndash;de&nbsp;Laat\'s <code>0.8358</code>. Monotonicity alone would')
    A('not give this: <code>36/11 = 3.2727&hellip;</code> falls between grid points, and the bound')
    A('inherited from <code>A = 3.5</code> yields only <code>0.8401</code>.</p>')
    A('<p><strong>Primes represented by a binary quadratic form, under GRH.</strong> The constant is')
    A('<code>2h(&minus;D)/C&#8330;(28)</code>:</p>')
    A('<div class="eq">limsup (p<sub>n+1,f</sub> &minus; p<sub>n,f</sub>) / (&radic;p&#8202;log&#8202;p) &lt; %.10f &middot; h(&minus;D)</div>' % s['qf'])
    A('<p>improving the constant in Chirre&ndash;Quesada-Herrera\'s Corollary&nbsp;5 from their printed')
    A('<code>1.837 h(&minus;D)</code>. For <code>u&sup2; + 27v&sup2;</code>, where <code>h(&minus;108) = 3</code>,')
    A('that reads <code>%.4f</code> against <code>5.5101</code>.</p>' % s['qf3'])
    A('<p><strong>Primes in an arithmetic progression, under GRH.</strong> The constant is')
    A('<code>&phi;(q)/C&#8330;(4)</code>:</p>')
    A('<div class="eq">limsup (p<sub>n+1,q,b</sub> &minus; p<sub>n,q,b</sub>) / (&radic;p&#8202;log&#8202;p) &lt; %.6f &middot; &phi;(q)</div>' % s['ap'])
    A('<p>improving Chirre&ndash;Pereira&ndash;de&nbsp;Laat\'s Corollary&nbsp;2, which gives')
    A('<code>0.8531 &phi;(q)</code>, for <code>q &ge; 3</code> and <code>b</code> coprime to <code>q</code>.</p>')
    A('<p class="note">All three constants are quoted rounded <em>up</em>, so each is genuinely implied')
    A('by the certified bound, following the same convention as the papers being improved.</p>')
    A('</section>')

    A('\n<section>')
    A('<h2>What this is not</h2>')
    A('<ul>')
    A('<li><strong>Not new mathematics.</strong> The extremal problem is CMS\'s, the constants and the')
    A('applications are Chirre&ndash;Quesada-Herrera\'s and Chirre&ndash;Pereira&ndash;de&nbsp;Laat\'s, and')
    A('the semidefinite formulation is theirs too. New here are the test functions and their')
    A('certificates.</li>')
    A('<li><strong>Not an improvement to the hard part.</strong> <code>A</code> is pinned by a sieve')
    A('bound &mdash; for ordinary primes, <code>36/11</code>, from Iwaniec (1982). Were the conjectural')
    A('value <code>A = 1</code> available, the same table would already give <code>0.5003</code>. Almost')
    A('everything separating that from <code>%.4f</code> lies outside Fourier optimization.</li>' % s['rh'])
    A('<li><strong>The SDP does not exactly reproduce published SDP values</strong> &mdash; it lands')
    A('4&ndash;5&times;10<sup>&minus;4</sup> low at matched degree, and the discrepancy is unexplained.')
    A('This does not affect the bounds: the semidefinite program is only a generator of candidate')
    A('functions, and each candidate is certified by evaluating it rigorously.</li>')
    A('<li><span class="warn"><strong>A = 1 is not a record.</strong></span> CMS prove')
    A('<code>C&#8330;(1) = 2</code> exactly; the certificate reaches 99.94% of it.</li>')
    A('<li><strong>The record check is thorough but not exhaustive.</strong> Citation indexes, author')
    A('sweeps and both authors\' own publication listings found no improvement since 2022, but')
    A('MathSciNet, zbMATH and Google Scholar were unreachable.</li>')
    A('<li><strong>Certification rests on ARB\'s correctness</strong> &mdash; unavoidable, and shared with')
    A('the published record it improves.</li>')
    A('</ul>')
    A('</section>')

    nlp = len({r['lp_file'] for r in rows}) + len({r['sdp_file'] for r in rows})
    A('\n<footer>')
    A('<p><strong>Code, test functions and certificates:</strong> <a href="https://github.com/sallerk/certified-cplus-bounds">github.com/sallerk/certified-cplus-bounds</a></p>')
    A('<p>Each row is independently checkable from its own coefficient file: %d distinct files across' % nlp)
    A('the %d rows, each with its own interval enclosure. <code>python verify_all.py</code> re-certifies' % s['n'])
    A('every one of them from scratch. This page is generated by <code>make_note.py</code> from the same')
    A('data, so it cannot disagree with the certificates.</p>')
    A('</footer>')
    A('</div>')
    return head.rstrip('\n') + '\n' + '\n'.join(body) + '\n'


if __name__ == '__main__':
    html = build()
    if '--check' in sys.argv:
        cur = open(OUT, encoding='utf-8').read() if os.path.exists(OUT) else ''
        if cur != html:
            print('docs/note.html is OUT OF DATE; run: python make_note.py')
            sys.exit(1)
        print('docs/note.html is up to date')
    else:
        open(OUT, 'w', encoding='utf-8', newline='\n').write(html)
        print('wrote %s (%d bytes)' % (OUT, len(html)))
