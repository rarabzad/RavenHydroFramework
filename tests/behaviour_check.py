"""Behaviour tests of the manual's table "Behaviour tests" that involve reservoirs and ensembles, on the full Liard example:

  reservoir  a 50 km2 lake with seepage on the largest coupled HRU of subbasin 52 (:ReservoirExchange), two years, plus a
             pumping well and a monitoring well: the seepage MODFLOW receives each step must equal the volume that left
             the lake the step before, all of it delivered; exact water accounts on both sides
  ensemble   the same set-up as a two-member Monte Carlo ensemble with a fixed parameter: both members must equal each
             other and the single run (groundwater budget, heads, hydrographs, reservoir stages)

Usage: python3 behaviour_check.py      (needs RAVEN_MF6_LIB)
"""
import os, sys, re, json, shutil, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import testconfig as tc
import numpy as np, pandas as pd
from shapely.geometry import shape

tc.need_mf6()
B = tc.LIARD
DAYS = 731


def setup(d):
    tc.copy_case(B, d, duration=DAYS)
    rvh = open(os.path.join(d, 'Liard.rvh')).read()
    hdr = [l for l in rvh.split('\n') if ':Attributes' in l and 'AQUIFER_PROFILE' in l][0].replace(',', ' ').split()[1:]
    ib, ia, ie, iar = (hdr.index(k) + 1 for k in ('BASIN_ID', 'AQUIFER_PROFILE', 'ELEVATION', 'AREA'))
    rows = [l.replace(',', ' ').split() for l in rvh.split(':HRUs')[1].split(':EndHRUs')[0].split('\n')]
    cand = sorted([r for r in rows if r and r[0].isdigit() and r[ib] == '52' and r[ia] != '[NONE]'], key=lambda r: -float(r[iar]))
    lid, elev = int(cand[0][0]), float(cand[0][ie])        # lake on the largest coupled HRU
    wid, wel = int(cand[1][0]), float(cand[1][ie])         # well in the second largest
    open(os.path.join(d, 'Liard.rvh'), 'a').write(
        '\n:Reservoir TestLake\n  :SubBasinID 52\n  :HRUID %d\n  :Type RESROUTE_STANDARD\n  :WeirCoefficient 0.6\n'
        '  :CrestWidth 20.0\n  :MaxDepth 8.0\n  :LakeArea 5.0E7\n  :AbsoluteCrestHeight %.2f\n  :SeepageParameters 0.05 %.2f\n'
        ':EndReservoir\n' % (lid, elev, elev - 5))
    geo = [f for f in json.load(open(os.path.join(d, 'Liard_HRUs_standin.geojson')))['features']
           if int(f['properties']['HRU_ID']) == wid][0]
    pt = shape(geo['geometry']).representative_point()
    open(os.path.join(d, 'Liard.rvg'), 'a').write(
        '\n:ReservoirExchange\n:Wells\n  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_TOP, SCREEN_BOT\n'
        '  1, PW_1, %.5f, %.5f, %.1f, %.1f\n:EndWells\n:ObservationWells\n  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_ELEV\n'
        '  101, OW_NEAR, %.5f, %.5f, %.1f\n:EndObservationWells\n' % (pt.y, pt.x, wel - 5, wel - 25, pt.y, pt.x + 0.0044, wel - 15))
    n = 7306
    open(os.path.join(d, 'Liard.rvt'), 'a').write('\n:WellRate 1 m3/d\n  1985-10-01 00:00:00 1.0 %d\n' % n + '  -20000.0\n' * n + ':EndWellRate\n')
    return lid


ok_all = True
d = tc.work('behaviour', 'reservoir')
lid = setup(d)
msg = tc.run(d)
print('reservoir: lake on HRU %d, 50 km2, %d days: %s' % (lid, DAYS, msg))
if 'Successful' not in msg:
    sys.exit(1)
A = tc.audit.audit(d + '/out')
b = pd.read_csv(d + '/out/GWBudget.csv')
sent, deliv = b['reservoir seepage sent [m3/d]'], b['reservoir seepage delivered [m3/d]']
checks = [
    ('recharge identity', abs(A['recharge: (Raven - MF6 - dStore - toCHD)/Raven']), 1e-9),
    ('MODFLOW cumulative residual', abs(A['MF6 cumulative residual / gross flow']), 1e-8),
    ('river bookkeeping', abs(A['river: applied - net MF6 - carry change (rel)']), 1e-8),
    ('reservoir: sent(t) - awaiting(t-1)', abs(A['reservoir: sent(t) - awaiting(t-1) (rel)']), 1e-8),
    ('reservoir: (delivered + not supplied) / sent - 1', abs(A['reservoir: (delivered + not supplied) / sent'] - 1), 1e-9),
    ('unconverged steps', A['unconverged steps'], 0.5),
    ('Raven balance error [mm]', A['Raven MB error max [mm]'], 1.0),
]
for name, v, lim in checks:
    ok = v < lim; ok_all &= ok
    print('  %-50s %.2e  %s' % (name, v, 'ok' if ok else 'FAIL'))
print('  seepage sent: mean %.0f m3/d, total %.3g m3; delivered/sent %.12f; wells delivered/specified %.3f' % (
    sent.mean(), sent.sum(), deliv.sum() / max(sent.sum(), 1e-30), b['wells actual [m3/d]'].sum() / min(b['wells specified [m3/d]'].sum(), -1e-30)))

# ensemble: two identical members against the single run
e = tc.work('behaviour', 'ensemble')
shutil.rmtree(e, ignore_errors=True)
shutil.copytree(d, e, ignore=shutil.ignore_patterns('out', '*.gwcache'))
t = open(e + '/Liard.rvi').read()
open(e + '/Liard.rvi', 'w').write(t + '\n:EnsembleMode ENSEMBLE_MONTECARLO 2\n')
open(e + '/Liard.rve', 'w').write(':OutputDirectoryFormat ./ens_*/\n:ParameterDistributions\n'
                                   '  BASEFLOW_COEFF SOIL SLOW_RES 2.258343E-02 DIST_UNIFORM 2.258343E-02 2.258343E-02\n'
                                   ':EndParameterDistributions\n')
r = subprocess.run([tc.RAVEN, 'Liard', '-o', 'out/'], cwd=e, capture_output=True, text=True, env=tc.ENV)
emsg = (re.findall(r'Exiting Gracefully: (.*)', r.stdout + r.stderr) or ['(no exit line)'])[-1]
print('ensemble: 2 Monte Carlo members: %s' % emsg)
worst = 0.0
for m in ('ens_1', 'ens_2'):
    for f in ('GWBudget.csv', 'GWHeads.csv', 'Hydrographs.csv', 'ReservoirStages.csv', 'GWHRUState.csv'):
        pa, pb = os.path.join(d, 'out', f), os.path.join(e, m, f)
        if not os.path.exists(pa):
            continue
        if not os.path.exists(pb):
            print('  %s/%s missing' % (m, f)); ok_all = False; continue
        a, c = pd.read_csv(pa), pd.read_csv(pb)
        num = [k for k in a.columns if k in c.columns and pd.api.types.is_numeric_dtype(a[k]) and pd.api.types.is_numeric_dtype(c[k])]
        dev = max(float(((a[k] - c[k]).abs() / (a[k].abs().max() + 1e-9)).max()) for k in num)
        worst = max(worst, dev)
        if len(a) != len(c):
            print('  %s/%s: %d rows instead of %d' % (m, f, len(c), len(a))); ok_all = False
print('  largest relative difference of a member from the single run: %.1e  %s' % (worst, 'ok' if worst < 1e-12 else 'FAIL'))
ok_all &= worst < 1e-12 and 'Successful' in emsg
print('behaviour checks: %s' % ('ALL PASSED' if ok_all else 'FAILED'))
sys.exit(0 if ok_all else 1)
