"""A USGS steady-state example in the existing-model mode: the stream-capture model of the MODFLOW 6 examples
(ex-gwf-capture: 40 x 20 cells of 250 m, one layer, Newton, seconds as time unit, one steady-state stress period of one
second, a river (RIV), pumping wells (WEL), specified heads and recharge). A steady-state model with one stress period has
no physical time: Raven stretches the period over its run, so that every Raven day is a steady-state solution.

The Liard HRUs of subbasin 52 are laid over the active cells (:OverlapWeights; see usgs_common.py). Raven takes over the
recharge package and the river package carries its water to Raven's reaches (subbasin of the dominant HRU of each cell).

Cases:
  published   the model's own solver (OUTER_DVCLOSE 1e-8 m, 100 outer iterations): days on which MODFLOW does not meet
              that tolerance are flagged; accounts exact on the others
  solver      the same with a copy whose solver stops at 1e-6 m, with up to 500 outer and 300 inner iterations and
              under-relaxation: every day converges, accounts exact
  restart     'solver' split at day 30 (restart from solution.rvc) against the continuous run

Usage: python3 usgs_capture.py PATH      (PATH: mf6examples.zip, its unpacked folder or the ex-gwf-capture folder)
"""
import os, sys, re, shutil, zipfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import testconfig as tc
import usgs_common as U
import pandas as pd

tc.need_mf6()
src = os.path.abspath(sys.argv[1])
if os.path.isdir(os.path.join(src, 'ex-gwf-capture')):   # the unpacked examples folder
    src = os.path.join(src, 'ex-gwf-capture')
W = os.path.dirname(tc.work('usgs', 'x'))
M = os.path.join(W, 'ex-gwf-capture')
if not os.path.exists(os.path.join(M, 'mfsim.nam')):
    if src.endswith('.zip'):
        with zipfile.ZipFile(src) as z:
            for n in z.namelist():
                if n.startswith('ex-gwf-capture/'):
                    z.extract(n, W)
    else:
        shutil.copytree(src, M)
MS = os.path.join(W, 'capture-solver')
if not os.path.exists(os.path.join(MS, 'mfsim.nam')):
    shutil.copytree(M, MS)
    ims = open(os.path.join(MS, 'ex-gwf-capture.ims')).read()
    ims = re.sub(r'OUTER_DVCLOSE\s+\S+', 'OUTER_DVCLOSE  1.0E-6', ims); ims = re.sub(r'OUTER_MAXIMUM\s+\S+', 'OUTER_MAXIMUM  500\n  UNDER_RELAXATION  DBD', ims)
    ims = re.sub(r'INNER_MAXIMUM\s+\S+', 'INNER_MAXIMUM  300', ims)
    open(os.path.join(MS, 'ex-gwf-capture.ims'), 'w').write(ims)
DAYS = 60
RVG = ':RavenTakesOver RCHA_0\n:StreamPackage RIV-1 DOMINANT_HRU\n:SeepageToRaven ON'


def leakance(d, v):
    """seepage leakance of both aquifer profiles: with the Liard value of 1/d the water table of this thin, flat-topped
    aquifer reaches the land surface on the wettest days and MODFLOW's Newton solver does not converge on two of them;
    0.01/d (suited to its conductivity of about 7e-5 m/s) converges every day"""
    p = d + '/Liard.rvp'; t = open(p).read()
    t = re.sub(r'(\n\s*(VALLEY_ALLUVIUM|UPLAND_TILL),\s*[\d.]+,\s*[\d.]+,\s*)[\d.]+', lambda m: m.group(1) + str(v), t)
    open(p, 'w').write(t)


def make(name, model, days, start='1985-10-01', rvc=None):
    d = U.make_case(W, name, model, os.path.join(model, 'ex-gwf-capture.dis'), RVG, days, start, rvc)
    leakance(d, 0.01)
    return d

ok_all = True; res = {}
for name, model in (('published', M), ('solver', MS)):
    d = make('cap_' + name, model, DAYS)
    msg = tc.run(d, timeout=7200)
    good = 'Successful' in msg
    line = '%-10s %s' % (name, msg[:120])
    if good:
        A, b, chk = U.audit(d); res[name] = b
        steady = 'steady-state model' in open(d + '/out/GWModelSummary.txt').read()
        line += '\n      steady-state model recognised: %s' % steady; good &= steady
        conv = b['converged'] > 0
        for n, v, lim in chk:
            if name == 'published' and n == 'unconverged steps':
                line += '\n      %-30s %d (the model\'s own tolerance of 1e-8 m; flagged, as designed)' % (n, v); continue
            if name == 'published' and n == 'MODFLOW cumulative residual':
                bb = b[conv]; e = bb['MF6 balance error [m3/d]'].abs().sum() / max(bb['MF6 recharge [m3/d]'].abs().sum(), 1)
                line += '\n      %-30s %.2e on the converged days %s' % ('MODFLOW residual', e, 'ok' if e < 1e-6 else 'FAIL'); good &= e < 1e-6; continue
            good &= v < lim; line += '\n      %-30s %.2e %s' % (n, v, 'ok' if v < lim else 'FAIL')
        line += '\n      storage release (must be 0 in a steady model): %.1e m3/d; river to Raven %.4g m3/d, Raven recharge %.4g m3/d (means)' % (
            b['storage release [m3/d]'].abs().max(), (b['aquifer discharge to rivers [m3/d]'] - b['river leakage to aquifer [m3/d]']).mean(), b['Raven recharge [m3/d]'].mean())
        good &= b['storage release [m3/d]'].abs().max() == 0
    print(line + ('' if good else '   <== FAIL'), flush=True); ok_all &= good
a = make('cap_restart_a', MS, 30)
ma = tc.run(a)
bdir = make('cap_restart_b', MS, DAYS - 30, start='1985-10-31', rvc=a + '/out/solution.rvc')
mb = tc.run(bdir)
if 'Successful' in ma and 'Successful' in mb and 'solver' in res:
    c = res['solver']; h = pd.read_csv(bdir + '/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]'])
    m = c.merge(h, on='date', suffixes=('_c', '_r'))
    cols = ['MF6 recharge [m3/d]', 'aquifer discharge to rivers [m3/d]', 'seepage to Raven [m3/d]', 'other model packages [m3/d]']
    dev = max(((m[x + '_c'] - m[x + '_r']).abs() / (m[x + '_c'].abs().max() + 1e3)).max() for x in cols)
    print('restart    day 30 -> %d days compared, largest deviation %.1e %s' % (len(m), dev, 'ok' if dev < 1e-4 else 'FAIL'))
    ok_all &= dev < 1e-4
else:
    print('restart    FAILED: %s / %s' % (ma, mb)); ok_all = False
print('capture checks: %s' % ('ALL PASSED' if ok_all else 'FAILED'))
sys.exit(0 if ok_all else 1)
