"""Runs the whole verification suite and writes one report.

    python3 run_all.py            # everything (needs RAVEN_MF6_LIB; about an hour on one core)
    python3 run_all.py --nolib    # only the checks that do not need the MODFLOW 6 library (a few minutes)
    python3 run_all.py --quick    # everything, with smaller randomized rounds
    python3 run_all.py --only thiem,edge    # selected steps (names as in the report)

The USGS example-model steps need the MODFLOW 6 examples (mf6examples.zip from
https://github.com/MODFLOW-ORG/modflow6-examples/releases): export RAVEN_USGS_EXAMPLES=/path/to/mf6examples.zip

Settings come from testconfig.py (RAVEN_EXE, RAVEN_MF6_LIB, RAVEN_TEST_WORK). The report is written to
<scratch>/run_all_report.md, the full output of every step to <scratch>/run_all_logs/. The exit code is the number of
failed steps.

Each step is judged by what its script prints: a step passes when the script finishes without error and its own pass
criterion holds (all cases pass, all random set-ups pass, residuals and restart deviations within the manual's limits).
Steps that only report numbers for the manual (Thiem slopes, Liard grid table) pass when every run finishes and is
audited; their numbers are copied into the report for comparison with the manual.
"""
import os, sys, re, json, time, subprocess, glob

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import testconfig as tc   # noqa: E402

NOLIB = '--nolib' in sys.argv
QUICK = '--quick' in sys.argv
ONLY = None
for i, a in enumerate(sys.argv):
    if a == '--only' and i + 1 < len(sys.argv):
        ONLY = set(sys.argv[i + 1].split(','))

LOGS = os.path.join(tc.WORK, 'run_all_logs')
os.makedirs(LOGS, exist_ok=True)
PY = sys.executable
ROOT = tc.SRC   # Raven's source files (the unit tests compile some of them)


def sh(cmd, cwd=HERE, timeout=None):
    """runs a command; returns (return code, combined output)"""
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=tc.ENV, timeout=timeout,
                           shell=isinstance(cmd, str))
        return r.returncode, r.stdout + r.stderr
    except subprocess.TimeoutExpired as ex:
        out = (ex.stdout or b'')
        out = out.decode() if isinstance(out, bytes) else out
        return -9, out + '\n(timeout after %d s)' % timeout


# ---- pass criteria: each takes the output text and returns (passed, one-line summary) ----
def last_line(pat):
    def judge(out):
        m = re.findall(pat, out)
        return (bool(m), m[-1] if m else 'expected line not found: ' + pat)
    return judge


def n_of_n(pat):
    """'<k> of <n> ...' with k == n"""
    def judge(out):
        m = re.findall(pat, out)
        if not m:
            return False, 'summary line not found'
        k, n = map(int, m[-1][:2])
        return k == n, '%d of %d' % (k, n)
    return judge


def zero_failed(out):
    m = re.findall(r'(\d+) cases, (\d+) checks, (\d+) failed', out)
    if not m:
        return False, 'summary line not found'
    c, k, f = map(int, m[-1])
    return f == 0, '%d cases, %d checks, %d failed' % (c, k, f)


def no_word(word, summary_pat=None):
    def judge(out):
        bad = [l for l in out.split('\n') if word in l]
        s = re.findall(summary_pat, out)[-1] if summary_pat and re.findall(summary_pat, out) else ''
        return (not bad, ('%d lines with %s' % (len(bad), word)) if bad else ('no %s %s' % (word, s)).strip())
    return judge


def small_numbers(pat, limit, what):
    """every number captured by pat must be below limit"""
    def judge(out):
        v = [abs(float(x)) for x in re.findall(pat, out)]
        if not v:
            return False, 'no %s found' % what
        return max(v) < limit, '%d %s, largest %.1e (limit %.0e)' % (len(v), what, max(v), limit)
    return judge


def fuzz_json(paths_fn):
    """all randomized set-ups of the given results files must pass (restart comparisons included)"""
    def judge(out):
        P = paths_fn(); R = []
        for p in P:
            if not os.path.exists(p):
                return False, 'no results file ' + p
            R += json.load(open(p))
        k = sum(1 for r in R if str(r.get('pass')) == 'True')   # stored as JSON true or as the text 'True'
        nr = sum(1 for r in R if 'hotstart_dev' in r or r.get('restart'))
        return k == len(R) and len(R) > 0, '%d of %d set-ups pass, %d with a restart (%s)' % (
            k, len(R), nr, ', '.join(os.path.basename(p) for p in P))
    return judge


def thiem_judge(out):
    # every run must be audited (no failed steps); slopes are reported for the manual's Thiem table. Grids whose cells
    # are larger than the fitting window (regular 1000 m: four cells > 3 km) have no cells to fit and are not judged.
    fails = re.findall(r'failed (\d+)', out)
    fits = [l for l in out.split('\n') if 'slope error' in l]
    fitted = [l for l in fits if 'cells used 0' not in l]
    ok = bool(fitted) and all(f == '0' for f in fails)
    return ok, '%d grids fitted (%d with no cells in the window), %d runs with failed steps' % (
        len(fitted), len(fits) - len(fitted), sum(1 for f in fails if f != '0'))


# ---- the steps: (name, needs MODFLOW, command, pass criterion, timeout [s]) ----
# the USGS example models (mf6examples.zip of the MODFLOW 6 examples, or its unpacked folder): set RAVEN_USGS_EXAMPLES
USGS = os.environ.get('RAVEN_USGS_EXAMPLES', '')
FZ = [('3', '8')] if QUICK else [('11', '24'), ('21', '24'), ('22', '24')]     # (seed, set-ups) of the grid rounds
FZL = [('12', '4')] if QUICK else [('12', '12'), ('13', '12')]                  # rounds with a reservoir and lake refinement
FE = ('5', '6') if QUICK else ('5', '20')
cxx = 'g++ -std=c++11 -O2 -I%s' % ROOT
STEPS = [
    ('gridtest', False, '%s gridtest.cpp %s/GWGrid.cpp %s/GWGeometry.cpp -o %s/gridtest && %s/gridtest'
     % (cxx, ROOT, ROOT, LOGS, LOGS), last_line(r'ALL PASSED: 0 failures'), 600),
    ('geomtest', False, '%s geomtest.cpp %s/GWGeometry.cpp -o geomtest && %s geomtest_check.py'
     % (cxx, ROOT, PY), last_line(r'ALL PASSED'), 600),
    ('crs', False, 'cd existing && %s crstest.cpp %s/GWGeometry.cpp -o crstest && %s crs_check.py'
     % (cxx, ROOT, PY), small_numbers(r'worst ([0-9.e+-]+) m', 1e-3, 'worst position difference [m]'), 600),
    ('mf6files', False, [PY, 'mf6files_check.py'], zero_failed, 3600),
    ('mf6lib', True, [PY, 'mf6lib_check.py'], last_line(r'library checks: ALL PASSED'), 1800),
    ('edge_nolib', False, [PY, 'edge.py', '--nolib'], n_of_n(r'(\d+) of (\d+) cases pass'), 3600),
    ('edge', True, [PY, 'edge.py'], n_of_n(r'(\d+) of (\d+) cases pass'), 4 * 3600),
    ('thiem', True, '%s thiem/run_thiem.py && %s thiem/fit_thiem.py' % (PY, PY), thiem_judge, 6 * 3600),
    ('restart', True, [PY, 'restart_test.py'], small_numbers(r'restart deviation ([0-9.e+-]+)', 1e-4, 'restart deviations'), 4 * 3600),
    ('lake', True, [PY, 'lake_test.py'], no_word('FAIL'), 4 * 3600),
    ('liard_grids', True, [PY, 'liard_grids.py'], no_word('FAILED', r'results in .*'), 6 * 3600),
    ('behaviour', True, [PY, 'behaviour_check.py'], last_line(r'behaviour checks: ALL PASSED'), 3 * 3600),
    ('budget_xcheck', True, lambda: [PY, 'budget_xcheck.py'] + sorted(glob.glob(os.path.join(tc.WORK, 'liard_grids', '*', 'out'))
                                    + glob.glob(os.path.join(tc.WORK, 'lake', '*', 'out')) + glob.glob(os.path.join(tc.WORK, 'behaviour', 'reservoir', 'out'))),
     small_numbers(r'WORST relative difference .*: ([0-9.e+-]+)', 1e-6, 'worst differences'), 1800),
    ('fuzz_grids', True, ' && '.join('%s fuzz_grids.py %s %s' % (PY, s, n) for s, n in FZ),
     fuzz_json(lambda: [os.path.join(tc.WORK, 'fuzz', 'results_%s.json' % s) for s, _ in FZ]), 12 * 3600),
    ('fuzz_lakes', True, ' && '.join('%s fuzz_grids.py %s %s lakes' % (PY, s, n) for s, n in FZL),
     fuzz_json(lambda: [os.path.join(tc.WORK, 'fuzz', 'results_%s_lakes.json' % s) for s, _ in FZL]), 12 * 3600),
    ('ext_guards', True, [PY, 'existing/guards.py'], n_of_n(r'(\d+) of (\d+) behave as required'), 3 * 3600),
    ('ext_roundtrip', True, [PY, 'existing/roundtrip.py'],
     small_numbers(r'max relative difference ([0-9.e+-]+)', 1e-4, 'round-trip differences'), 3 * 3600),
    ('ext_restart', True, [PY, 'existing/restart_ext.py'],
     small_numbers(r'max deviation ([0-9.e+-]+)', 1e-4, 'restart deviations'), 3 * 3600),
    ('usgs_sagehen', True, lambda: [PY, 'existing/usgs_sagehen.py', USGS], last_line(r'Sagehen checks: ALL PASSED'), 6 * 3600),
    ('usgs_capture', True, lambda: [PY, 'existing/usgs_capture.py', USGS], last_line(r'capture checks: ALL PASSED'), 3 * 3600),
    ('ext_fuzz', True, [PY, 'existing/fuzz_ext.py', FE[0], FE[1]],
     fuzz_json(lambda: [os.path.join(tc.WORK, 'extmodel', 'fz_results_%s.json' % FE[0])]), 12 * 3600),
]

have_lib = bool(tc.MF6LIB) and os.path.exists(tc.MF6LIB)
if not NOLIB and not have_lib:
    print('RAVEN_MF6_LIB is not set to the MODFLOW 6 library: running only the checks that do not need it (--nolib).')
    NOLIB = True
if not os.path.exists(tc.RAVEN):
    sys.exit('Raven not found at %s: build it (make) or set RAVEN_EXE' % tc.RAVEN)

rows = []
t_all = time.time()
for name, needs_lib, cmd, judge, timeout in STEPS:
    if ONLY and name not in ONLY:
        continue
    if needs_lib and NOLIB:
        rows.append((name, 'skipped', 0, 'needs the MODFLOW 6 library', ''))
        continue
    if name.startswith('usgs') and not USGS:
        rows.append((name, 'skipped', 0, 'set RAVEN_USGS_EXAMPLES to mf6examples.zip (MODFLOW 6 examples)', ''))
        continue
    print('%-14s ...' % name, end=' ', flush=True)
    t0 = time.time()
    rc, out = sh(cmd() if callable(cmd) else cmd, timeout=timeout)
    dt = time.time() - t0
    open(os.path.join(LOGS, name + '.txt'), 'w').write(out)
    ok, summary = judge(out)
    if rc != 0:
        ok, summary = False, 'exit code %d; %s' % (rc, summary)
    status = 'PASS' if ok else 'FAIL'
    print('%s  %s  (%.0f s)' % (status, summary, dt), flush=True)
    tail = '\n'.join([l for l in out.strip().split('\n') if 'Warning' not in l and 'warn(' not in l][-25:])
    rows.append((name, status, dt, summary, tail))

# with --only, the report keeps the earlier results of the steps not run this time (each row shows when it ran)
RES = os.path.join(tc.WORK, 'run_all_results.json')
now = time.strftime('%Y-%m-%d %H:%M')
rows = [tuple(r) + (now,) for r in rows]                      # this run's steps
if ONLY and os.path.exists(RES):
    old = {r[0]: tuple(r) for r in json.load(open(RES))}
    new = {r[0]: r for r in rows}
    rows = [new.get(n, old.get(n)) for n in [s[0] for s in STEPS] if n in new or n in old]
json.dump([list(r) for r in rows], open(RES, 'w'), indent=1)
rep = os.path.join(tc.WORK, 'run_all_report.md')
with open(rep, 'w') as f:
    f.write('# Verification report\n\n')
    f.write('- Raven: `%s`\n- MODFLOW 6 library: `%s`\n- scratch folder: `%s`\n- date: %s\n- mode: %s\n\n' % (
        tc.RAVEN, tc.MF6LIB if have_lib else '(not used)', tc.WORK, time.strftime('%Y-%m-%d %H:%M'),
        'no library' if NOLIB else ('quick' if QUICK else 'full')))
    f.write('| Step | Result | Time [s] | Summary | Run |\n|---|---|---|---|---|\n')
    for name, status, dt, summary, _, when in rows:
        f.write('| %s | %s | %.0f | %s | %s |\n' % (name, status, dt, summary.replace('|', '/'), when))
    f.write('\nTotal time %.0f s. Full output of each step: `%s/<step>.txt`.\n' % (time.time() - t_all, LOGS))
    f.write('\nNumbers reported in the manual come from: `thiem` (Table "Thiem test"), `liard_grids` (Table "The Liard '
            'model on each grid"), `restart`, `lake` and the randomized steps (chapter "Verification").\n')
    for name, status, dt, summary, tail, when in rows:
        if tail:
            f.write('\n## %s (%s)\n\n```\n%s\n```\n' % (name, status, tail))
print('report:', rep)
sys.exit(sum(1 for r in rows if r[1] == 'FAIL'))
