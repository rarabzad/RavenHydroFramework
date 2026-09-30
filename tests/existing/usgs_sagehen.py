"""A real USGS model in the existing-model mode: the Sagehen Creek watershed model of the MODFLOW 6 examples
(ex-gwf-sagehen: Sagehen Creek, Sierra Nevada, California; 73 x 81 cells of 90 m, one layer, Newton, a steady-state day
and 398 daily periods, UZF with rejected infiltration and groundwater discharge moved to SFR by the MVR package, land
surface drains moved to SFR, three specified heads).

The model is not georeferenced, so the Liard Raven model is laid over it: the coupled HRUs of subbasin 52 are given zones
of the model's active cells (weights 1.0 in an :OverlapWeights table) and .rvh areas equal to their zones. Raven takes over
the UZF package (Raven's recharge replaces the model's infiltration); the model keeps its drains, streams and mover.

Cases (each audited; the model's own packages enter the balance check). The published model solves loosely (OUTER_DVCLOSE
0.03 m, INNER_RCLOSE 1000 m3/d); the other cases use a copy with a tight solver, so that MODFLOW's own budget closes:
  published     :RavenTakesOver UZF-1, :SeepageToRaven OFF, the model's own solver settings
  takeover      the same on the tight copy (the model's drains and streams carry the discharge)
  seepage       the same with :SeepageToRaven ON (Raven's seepage drains beside the model's)
  restart       'takeover' split at day 200 (restart from solution.rvc) against the continuous run
  stream_mover  :StreamPackage DRN-1 while DRN-1 sends its water to SFR through the mover: must be refused
  keep_uzf      UZF-1 kept (:KeepPackages): must run (Raven's recharge is then added beside the model's infiltration)

Usage: python3 usgs_sagehen.py PATH      (PATH: mf6examples.zip, its unpacked folder or the ex-gwf-sagehen folder;
       https://github.com/MODFLOW-ORG/modflow6-examples/releases/download/current/mf6examples.zip)
"""
import os, sys, re, shutil, subprocess, zipfile, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import testconfig as tc
import numpy as np, pandas as pd

tc.need_mf6()
src = os.path.abspath(sys.argv[1])
if os.path.isdir(os.path.join(src, 'ex-gwf-sagehen')):   # the unpacked examples folder
    src = os.path.join(src, 'ex-gwf-sagehen')
W = tc.work('usgs', 'x'); W = os.path.dirname(W)
M = os.path.join(W, 'ex-gwf-sagehen')
if not os.path.exists(os.path.join(M, 'mfsim.nam')):
    if src.endswith('.zip'):
        with zipfile.ZipFile(src) as z:
            for n in z.namelist():
                if n.startswith('ex-gwf-sagehen/'):
                    z.extract(n, W)
    else:
        shutil.copytree(src, M)
DAYS = 399
# the published model solves loosely (OUTER_DVCLOSE 0.03 m, INNER_RCLOSE 1000 m3/d), so MODFLOW's own budget can be off by up
# to about 1% on a day; a copy with a tight solver shows that the coupling itself is exact
MT = os.path.join(W, 'sagehen-tight')
if not os.path.exists(os.path.join(MT, 'mfsim.nam')):
    shutil.copytree(M, MT)
    ims = open(os.path.join(MT, 'ex-gwf-sagehen.ims')).read()
    ims = re.sub(r'OUTER_DVCLOSE\s+\S+', 'OUTER_DVCLOSE  1.0E-5', ims); ims = re.sub(r'INNER_DVCLOSE\s+\S+', 'INNER_DVCLOSE  1.0E-6', ims)
    ims = re.sub(r'inner_rclose\s+\S+', 'inner_rclose  1.0E-3', ims, flags=re.I)
    open(os.path.join(MT, 'ex-gwf-sagehen.ims'), 'w').write(ims)

import usgs_common as U
nrow, ncol, delr, delc, act = U.active_cells(os.path.join(M, 'ex-gwf-sagehen.dis'))
print('Sagehen: %d x %d cells of %.0f m, %d active (%.2f km2)' % (nrow, ncol, delr, len(act), len(act) * delr * delc / 1e6))


def case(name, rvg_extra, days=DAYS, start='1985-10-01', rvc=None, model=None):
    m = model or MT
    return U.make_case(W, name, m, os.path.join(m, 'ex-gwf-sagehen.dis'), rvg_extra, days, start, rvc)


audit = U.audit


ok_all = True
res = {}
for name, extra, expect in [
        ('published', ':RavenTakesOver UZF-1\n:SeepageToRaven OFF', 'ok'),
        ('takeover', ':RavenTakesOver UZF-1\n:SeepageToRaven OFF', 'ok'),
        ('seepage', ':RavenTakesOver UZF-1\n:SeepageToRaven ON', 'ok'),
        ('stream_mover', ':RavenTakesOver UZF-1\n:StreamPackage DRN-1 DOMINANT_HRU\n:SeepageToRaven OFF', 'error'),
        ('keep_uzf', ':KeepPackages UZF-1\n:SeepageToRaven OFF', 'ok')]:
    d = case(name, extra, days=(DAYS if expect == 'ok' else 30), model=(M if name == 'published' else MT))
    msg = tc.run(d, timeout=7200)
    good = ('Successful' in msg) == (expect == 'ok')
    line = '%-13s expected %-5s: %s' % (name, expect, msg[:150])
    if good and expect == 'ok':
        # the published solver settings: MODFLOW's own discrepancy (up to ~1% a day) limits the cumulative residual
        A, b, chk = audit(d, resid=(2e-3 if name == 'published' else 1e-8)); res[name] = b
        for n, v, lim in chk:
            if name == 'keep_uzf' and n == 'unconverged steps':
                # recharge counted twice (the model's infiltration and Raven's): in the snowmelt peak a day may miss the
                # tight tolerance within 300 iterations; it is flagged and its exchange limited, as designed
                line += '\n      %-30s %d (flagged; accounts still exact)' % (n, v); good &= v <= 3; continue
            good &= v < lim; line += '\n      %-30s %.2e %s' % (n, v, 'ok' if v < lim else 'FAIL')
        line += '\n      Raven recharge %.4g m3/d, model packages %.4g m3/d, storage %.4g m3/d, seepage to Raven %.4g m3/d (means)' % (
            b['Raven recharge [m3/d]'].mean(), b['other model packages [m3/d]'].mean(), b['storage release [m3/d]'].mean(),
            b['seepage to Raven [m3/d]'].mean())
    print(line + ('' if good else '   <== FAIL'), flush=True); ok_all &= good

# restart inside the daily periods: 200 days, then 199 days from solution.rvc
a = case('restart_a', ':RavenTakesOver UZF-1\n:SeepageToRaven OFF', days=200)
ma = tc.run(a, timeout=7200)
bdir = case('restart_b', ':RavenTakesOver UZF-1\n:SeepageToRaven OFF', days=DAYS - 200, start='1986-04-19', rvc=a + '/out/solution.rvc')
mb = tc.run(bdir, timeout=7200)
if 'Successful' in ma and 'Successful' in mb and 'takeover' in res:
    c = res['takeover']; h = pd.read_csv(bdir + '/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]'])
    m = c.merge(h, on='date', suffixes=('_c', '_r'))
    cols = ['MF6 recharge [m3/d]', 'storage release [m3/d]', 'other model packages [m3/d]']
    dev = max(((m[x + '_c'] - m[x + '_r']).abs() / (m[x + '_c'].abs().max() + 1e3)).max() for x in cols)
    print('restart       day 200 -> %d days compared, largest deviation %.1e %s' % (len(m), dev, 'ok' if dev < 1e-4 else 'FAIL'))
    ok_all &= dev < 1e-4
else:
    print('restart       FAILED: %s / %s' % (ma, mb)); ok_all = False
print('Sagehen checks: %s' % ('ALL PASSED' if ok_all else 'FAILED'))
sys.exit(0 if ok_all else 1)
