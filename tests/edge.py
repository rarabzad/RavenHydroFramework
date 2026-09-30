# Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
# SPDX-License-Identifier: Artistic-2.0
"""Input errors and edge cases: every faulty set-up must stop with a clear message (no crash); every valid one must run.

Usage:  python3 edge.py            full runs (needs RAVEN_MF6_LIB)
        python3 edge.py --nolib    without the MODFLOW library: Raven checks every input and builds the model (for an
                                   existing model: reads, copies and rewrites it) before loading the library, so a valid
                                   case passes when it stops at the library load, and a faulty one must stop earlier.
Cases starting with E use the existing-model example; the others use the Liard example.
Results: <scratch>/edge/results.json
"""
import os, re, shutil, subprocess, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import testconfig as tc

NOLIB = '--nolib' in sys.argv
RV = tc.RAVEN
if not NOLIB:
    tc.need_mf6()
ENV = dict(tc.ENV, RAVEN_MF6_LIB='/nonexistent/libmf6.so') if NOLIB else tc.ENV
W = os.path.dirname(tc.work('edge', 'x'))
LIB_MSG = 'cannot load MODFLOW 6 library'
# the existing-model example reads its forcing from ../Liard_groundwater: provide that sibling in the scratch folder
if not os.path.exists(os.path.join(W, 'Liard_groundwater')):
    try: os.symlink(tc.LIARD, os.path.join(W, 'Liard_groundwater'))
    except OSError: shutil.copytree(tc.LIARD, os.path.join(W, 'Liard_groundwater'))


def base(name, days=10, src=None):
    d = os.path.join(W, name)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(src or tc.LIARD, d, ignore=shutil.ignore_patterns('out', 'output', '__pycache__', '*.gwcache'))
    x = open(f'{d}/Liard.rvi').read(); x = re.sub(r':Duration\s+\d+', f':Duration {days}', x); open(f'{d}/Liard.rvi', 'w').write(x)
    return d
def edit(d, fn, old, new, count=1, regex=False):
    p = f'{d}/{fn}'; s = open(p).read()
    if regex:
        s2, n = re.subn(old, new, s, count=count); assert n > 0, (fn, old); s = s2
    else:
        assert old in s, (fn, old[:50]); s = s.replace(old, new, count)
    open(p, 'w').write(s)
def app(d, fn, txt): open(f'{d}/{fn}', 'a').write(txt)
cases = {}
def case(name, expect, fn, days=10, src=None, verify=None):
    d = base(name, days, src); fn(d); cases[name] = (d, expect, verify)

EXT = tc.EXISTING
RES = ("\n:Reservoir TestLake\n  :SubBasinID 52\n  :HRUID {hru}\n  :Type RESROUTE_STANDARD\n  :WeirCoefficient 0.6\n  :CrestWidth 20.0\n"
       "  :MaxDepth 8.0\n  :LakeArea 2.0E7\n{crest}  :SeepageParameters 0.05 {seep}\n:EndReservoir\n")
def lake_hru(d):
    rvh = open(f'{d}/Liard.rvh').read()
    hdr = [l for l in rvh.split('\n') if ':Attributes' in l and 'AQUIFER_PROFILE' in l][0].replace(',', ' ').split()[1:]
    ib, ia, ie, iar = hdr.index('BASIN_ID') + 1, hdr.index('AQUIFER_PROFILE') + 1, hdr.index('ELEVATION') + 1, hdr.index('AREA') + 1
    rows = [l.replace(',', ' ').split() for l in rvh.split(':HRUs')[1].split(':EndHRUs')[0].split('\n')]
    r = max([r for r in rows if r and r[0].isdigit() and r[ib] == '52' and r[ia] != '[NONE]'], key=lambda r: float(r[iar]))
    return int(r[0]), float(r[ie])

# ---- generated models: faulty inputs
case('N01_missing_geojson', 'error', lambda d: edit(d, 'Liard.rvg', 'Liard_HRUs_standin.geojson', 'no_such_file.geojson'))
case('N02_bad_id_field', 'error', lambda d: edit(d, 'Liard.rvg', 'Liard_HRUs_standin.geojson  HRU_ID', 'Liard_HRUs_standin.geojson  NOT_A_FIELD'))
case('N03_unknown_profile', 'error', lambda d: edit(d, 'Liard.rvh', 'UPLAND_TILL', 'UNKNOWN_PROFILE'))
case('N04_aquiclude_top', 'error', lambda d: edit(d, 'Liard.rvp', 'UPLAND_TILL,     1, UPLAND_TILL, TO_BEDROCK, AQUIFER', 'UPLAND_TILL,     1, UPLAND_TILL, TO_BEDROCK, AQUICLUDE'))
case('N05_undefined_class', 'error', lambda d: edit(d, 'Liard.rvp', 'UPLAND_TILL,     1, UPLAND_TILL,', 'UPLAND_TILL,     1, NO_CLASS,'))
case('N06_row_length', 'error', lambda d: edit(d, 'Liard.rvp', '  SAND,        10.0,    1.0,    1.0E-5,       0.20,       0.35', '  SAND,        10.0,    1.0,    1.0E-5,       0.20'))
case('N07_missing_Kh', 'error', lambda d: edit(d, 'Liard.rvp', ':Attributes, K_HORIZ, K_VERT,', ':Attributes, K_H_TYPO, K_VERT,'))
case('N08_drain_GW_in_coupled', 'error', lambda d: edit(d, 'Liard.rvi', ':HydrologicProcesses', ':HydrologicProcesses\n  :Baseflow BASE_LINEAR GROUNDWATER SURFACE_WATER'))
case('N09_no_recharge_process', 'error', lambda d: edit(d, 'Liard.rvi', 'FAST_RESERVOIR \tGROUNDWATER', 'FAST_RESERVOIR \tSLOW_RESERVOIR'))
case('N10_bad_library', 'badlib', lambda d: app(d, 'Liard.rvg', ':MF6Library /no/such/libmf6.so\n'))
case('N11_well_outside', 'error', lambda d: app(d, 'Liard.rvg', ':Wells\n  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_TOP, SCREEN_BOT\n  1, W1, 40.0, -100.0, 10, 0\n:EndWells\n'))
case('N13_swx_bad_group', 'error', lambda d: app(d, 'Liard.rvg', ':SurfaceWaterExchange DEPRESSION NoSuchGroup LEAKANCE 0.01\n'))
case('N14_swx_bad_sv', 'error', lambda d: app(d, 'Liard.rvg', ':SurfaceWaterExchange LAKE_STORAGE AquiferHRUs LEAKANCE 0.01\n'))
case('N15_ghb_zero_cond', 'error', lambda d: app(d, 'Liard.rvg', ':GeneralHeadBoundary G1 Liard_rivers.geojson HEAD 1000 CONDUCTANCE 0\n'))
case('N16_hotstart_dims', 'error', lambda d: app(d, 'Liard.rvc', '\n:GWHeads 1 2 2\n 1 2 3 4\n:EndGWHeads\n'))
case('N21_bad_raster', 'error', lambda d: edit(d, 'Liard.rvg', 'Liard_dem_standin.asc', 'missing_dem.asc'))
case('N23_usg_command', 'error', lambda d: app(d, 'Liard.rvg', ':NameFile model.nam\n'))
case('N24_refine_lakes_no_reservoir', 'error', lambda d: app(d, 'Liard.rvg', ':GridType QUADTREE\n:GridRefinement LAKES 500\n'))
case('N25_refinement_too_large', 'error', lambda d: app(d, 'Liard.rvg', ':GridType QUADTREE\n:GridRefinement RIVERS 5000\n'))
case('N26_rotation_mesh', 'error', lambda d: app(d, 'Liard.rvg', ':GridType HRU_MESH\n:GridRotation 20\n'))
case('N27_gridfile_missing', 'error', lambda d: app(d, 'Liard.rvg', ':GridType FILE\n'))
case('N28_negative_smoothing', 'error', lambda d: edit(d, 'Liard.rvp', 'SEEPAGE_LEAKANCE, SOIL_ZONE_DEPTH', 'SEEPAGE_LEAKANCE, SOIL_ZONE_DEPTH, SEEPAGE_SMOOTHING_DEPTH') or
     edit(d, 'Liard.rvp', '1/d,              m\n', '1/d,              m, m\n') or edit(d, 'Liard.rvp', '1.0,              0.0\n', '1.0,              0.0, -1.0\n', count=2))
case('N29_reservoir_relative_stages', 'error', lambda d: app(d, 'Liard.rvh', RES.format(hru=lake_hru(d)[0], crest='', seep=-5.0)) or app(d, 'Liard.rvg', ':ReservoirExchange\n'))
# ---- generated models: valid edge cases
case('N12_obswell_outside', 'ok', lambda d: app(d, 'Liard.rvg', ':ObservationWells\n  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_ELEV\n  9, OW, 40.0, -100.0, 10\n:EndObservationWells\n'))
case('N17_huge_cells', 'either', lambda d: edit(d, 'Liard.rvg', r':GridCellSize\s+\S+', ':GridCellSize 200000.0', regex=True))
case('N18_gw_off', 'uncoupled', lambda d: edit(d, 'Liard.rvi', ':GroundwaterModel MODFLOW6', ':GroundwaterModel NONE'))
case('N19_tobedrock_middle', 'ok', lambda d: edit(d, 'Liard.rvp', 'GRAVEL, 15.0, AQUIFER, CLAY_TILL, 5.0, AQUITARD', 'GRAVEL, TO_BEDROCK, AQUIFER, CLAY_TILL, 5.0, AQUITARD'))
case('N20_unknown_rvg_cmd', 'ok', lambda d: app(d, 'Liard.rvg', ':NotACommand 1 2 3\n'))
case('N22_no_rasters', 'ok', lambda d: (edit(d, 'Liard.rvg', ':LandSurface', '#:LandSurface'), edit(d, 'Liard.rvg', ':BedrockSurface', '#:BedrockSurface')))
case('N30_reservoir_absolute_stages', 'ok', lambda d: (lambda h: app(d, 'Liard.rvh', RES.format(hru=h[0], crest='  :AbsoluteCrestHeight %.2f\n' % h[1], seep=h[1] - 5)))(lake_hru(d)) or app(d, 'Liard.rvg', ':ReservoirExchange\n:GridType HRU_MESH\n:GridRefinement LAKES 500\n'))
case('N31_mesh_obswell_outside_refined', 'ok', lambda d: app(d, 'Liard.rvg', ':GridType HRU_MESH\n:GridRefinement WELLS 400\n:ObservationWells\n  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_ELEV\n  9, OW, 58.2, -127.9, 10\n:EndObservationWells\n'))
case('N32_reservoir_hru_of_other_subbasin', 'error', lambda d: app(d, 'Liard.rvh', RES.format(hru=2745, crest='  :AbsoluteCrestHeight 1115.0\n', seep=1110.0)))   # HRU 2745 lies in subbasin 43 (crashed Raven 4.1.1)
case('P01_cell1000', 'ok', lambda d: edit(d, 'Liard.rvg', r':GridCellSize\s+\S+', ':GridCellSize 1000.0', regex=True), days=365)
case('P02_cell2000', 'ok', lambda d: None, days=365)
case('P03_cell4000', 'ok', lambda d: edit(d, 'Liard.rvg', r':GridCellSize\s+\S+', ':GridCellSize 4000.0', regex=True), days=365)
case('P04_halfday', 'ok', lambda d: edit(d, 'Liard.rvi', ':TimeStep', ':TimeStep 0.5\n#:TimeStep'), days=60)

# ---- existing models: faulty inputs (checked before MODFLOW is loaded)
case('E01_two_link_methods', 'error', lambda d: app(d, 'Liard.rvg', ':MF6GridFile cells.geojson CELL_ID\n'), src=EXT)
case('E02_generated_setting', 'error', lambda d: app(d, 'Liard.rvg', ':GridCellSize 1000\n'), src=EXT)
case('E03_grid_type', 'error', lambda d: app(d, 'Liard.rvg', ':GridType QUADTREE\n'), src=EXT)
case('E04_netcdf_with_weights', 'error', lambda d: (edit(d, 'Liard.rvg', ':MF6CRS', '# :MF6CRS'), edit(d, 'Liard.rvg', ':HRUGeometry', '# :HRUGeometry'),
     app(d, 'Liard.rvg', open(f'{d}/overlap_weights_alternative.rvg').read() + ':WriteNetCDFHeads\n')), src=EXT)
case('E05_start_before_model', 'error', lambda d: edit(d, 'Liard.rvi', ':StartDate          1985-10-01', ':StartDate          1985-09-01'), src=EXT)
case('E06_two_solutions', 'error', lambda d: edit(d, 'model/mfsim.nam', '  ims6  ext.ims  gwf_liard\n', '  ims6  ext.ims  gwf_liard\n  ims6  other.ims  other_model\n'), src=EXT)
case('E07_adaptive_time_steps', 'error', lambda d: edit(d, 'model/ext.tdis', '  TIME_UNITS  days\n', '  TIME_UNITS  days\n  ATS6  FILEIN  ext.ats\n'), src=EXT)
case('E08_noleap_calendar', 'error', lambda d: edit(d, 'Liard.rvi', ':StartDate', ':Calendar NOLEAP\n:StartDate'), src=EXT)
case('E09_model_in_working_folder', 'error', lambda d: (shutil.copytree(f'{d}/model', f'{d}/out/mf6'), edit(d, 'Liard.rvg', 'model/mfsim.nam', 'out/mf6/mfsim.nam')), src=EXT)
case('E10_restart_with_time_series', 'error', lambda d: (edit(d, 'Liard.rvi', ':StartDate          1985-10-01', ':StartDate          1985-11-04'),
     edit(d, 'model/gwf_liard.wel', 'BEGIN options\n', 'BEGIN options\n  TS6  FILEIN  wel.ts\n')), src=EXT)
case('E11_restart_with_tvk', 'error', lambda d: (edit(d, 'Liard.rvi', ':StartDate          1985-10-01', ':StartDate          1985-11-04'),
     edit(d, 'model/gwf_liard.npf', '  SAVE_FLOWS\n', '  SAVE_FLOWS\n  TVK6  FILEIN  npf.tvk\n')), src=EXT)
case('E12_usg_command', 'error', lambda d: app(d, 'Liard.rvg', ':GWRiverConnection\n'), src=EXT)
case('E13_period_not_whole_steps', 'error', lambda d: edit(d, 'model/ext.tdis', '       1.00000000  1       1.00000000', '       1.50000000  1       1.00000000'), src=EXT)

# ---- existing models: valid set-ups; the working copy is checked
def v_plain(d):
    w = f'{d}/out/mf6'; out = []
    sim = open(f'{w}/mfsim.nam').read()
    out.append(('mfsim.nam has CONTINUE', 'CONTINUE' in sim.upper()))
    t = open(f'{w}/ext.tdis').read()
    per = [l.split() for l in t.split('BEGIN PERIODDATA')[1].split('END PERIODDATA')[0].strip().split('\n')]
    out.append(('TDIS: periods split into daily steps, ending with the 100-day run', all(int(p[1]) == round(float(p[0])) for p in per) and abs(sum(float(p[0]) for p in per) - 100) < 1e-9))
    nam = open(f'{w}/gwf_liard.nam').read()
    out.append(('name file: RCH_MODEL and EVT_MODEL switched off', bool(re.search(r'(?mi)^#.*RCH_MODEL', nam)) and bool(re.search(r'(?mi)^#.*EVT_MODEL', nam))))
    out.append(('original model folder unchanged', open(f'{d}/model/mfsim.nam').read() == open(f'{EXT}/model/mfsim.nam').read()))
    return out
def v_restart(d, start='1985-11-04'):
    """expected result computed from the original model: drop the periods that end on or before Raven's start, shorten the
    one that holds it, renumber the package periods (new period 1 = the data in force at the start)"""
    import datetime
    w = f'{d}/out/mf6'; out = []
    per0 = [float(l.split()[0]) for l in open(f'{EXT}/model/ext.tdis').read().split('BEGIN perioddata')[1].split('END perioddata')[0].strip().split('\n')]
    off = (datetime.date.fromisoformat(start) - datetime.date(1985, 10, 1)).days
    drop, cum = 0, 0.0
    while cum + per0[drop] <= off: cum += per0[drop]; drop += 1
    first = per0[drop] - (off - cum)
    t = open(f'{w}/ext.tdis').read()
    per = [l.split() for l in t.split('BEGIN PERIODDATA')[1].split('END PERIODDATA')[0].strip().split('\n')]
    out.append(('TDIS: %d elapsed periods dropped, the current one shortened to %g days' % (drop, first),
                abs(float(per[0][0]) - first) < 1e-9 and int(per[0][1]) == round(first) and abs(float(per[1][0]) - per0[drop + 1]) < 1e-9))
    out.append(("TDIS: START_DATE_TIME set to Raven's start", start in t))
    blk = lambda s: {int(m.group(1)): m.group(2).strip() for m in re.finditer(r'(?s)BEGIN period\s+(\d+)[^\n]*\n(.*?)END period', s, re.I)}
    bo = blk(open(f'{EXT}/model/gwf_liard.wel').read()); bn = blk(open(f'{w}/gwf_liard.wel').read())
    out.append(('WEL: period k+%d of the model is period k of the copy' % drop, all(bn.get(k) == bo.get(k + drop) for k in range(2, 12) if k + drop in bo)))
    out.append(('WEL: new period 1 holds the data in force at the start', bn.get(1) == bo.get(max(k for k in bo if k <= drop + 1))))
    return out
case('E20_valid', 'ok', lambda d: None, days=100, src=EXT, verify=v_plain)
case('E21_valid_restart', 'ok', lambda d: edit(d, 'Liard.rvi', ':StartDate          1985-10-01', ':StartDate          1985-11-04'), days=100, src=EXT, verify=v_restart)
case('E22_valid_sim_file_other_name', 'ok', lambda d: (os.rename(f'{d}/model/mfsim.nam', f'{d}/model/run1.nam'), edit(d, 'Liard.rvg', 'model/mfsim.nam', 'model/run1.nam')), src=EXT,
     verify=lambda d: [('working copy holds mfsim.nam with CONTINUE', os.path.exists(f'{d}/out/mf6/mfsim.nam') and 'CONTINUE' in open(f'{d}/out/mf6/mfsim.nam').read().upper())])
case('E23_valid_output_inside_model_folder', 'ok', lambda d: None, src=EXT)   # run with -o model/out/ (below)

# ------------------------------------------------------------------------------------------------ run and judge
res = {}; npass = 0
for n, (d, exp, verify) in cases.items():
    out = 'model/out/' if n == 'E23_valid_output_inside_model_folder' else 'out/'
    os.makedirs(os.path.join(d, out), exist_ok=True)
    try:
        r = subprocess.run([RV, 'Liard', '-o', out], cwd=d, capture_output=True, text=True, timeout=1800, env=ENV); rc = r.returncode; so = r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        rc = 'timeout'; so = ''
    m = re.findall(r'Exiting Gracefully: (.*)', so) or re.findall(r'Error Statement: (.*)', so); msg = m[-1] if m else ''
    crash = (rc in (-11, -6, 134, 139)) or ('Segmentation' in so) or (rc == 'timeout')
    ok_run = 'Successful Simulation' in so
    reached = ok_run or (NOLIB and LIB_MSG in msg)
    if exp == 'error': good = (not crash) and (not ok_run) and msg != '' and LIB_MSG not in msg
    elif exp == 'ok': good = (not crash) and reached
    elif exp == 'uncoupled': good = ok_run
    elif exp == 'badlib': good = (not crash) and (not ok_run) and '/no/such/libmf6.so' in msg
    else: good = not crash
    checks = []
    if good and verify is not None:
        try: checks = verify(d)
        except Exception as e: checks = [('verification ran', False, repr(e)[:100])]
        good = good and all(c[1] for c in checks)
    if n == 'E23_valid_output_inside_model_folder' and good:
        nest = os.path.exists(os.path.join(d, 'model/out/mf6/out'))
        checks.append(('working copy did not copy itself', not nest)); good = good and not nest
    res[n] = dict(expect=exp, passed=bool(good), rc=str(rc), crash=bool(crash), msg=msg[:200], checks=[(c[0], bool(c[1])) for c in checks])
    npass += good
    print('%-40s %-9s %s  %s' % (n, exp, 'PASS' if good else 'FAIL', msg[:110]), flush=True)
    for c in checks:
        if not c[1]: print('      FAIL %s' % c[0])
json.dump(res, open(os.path.join(W, 'results.json'), 'w'), indent=1)
print('%d of %d cases pass%s; results in %s' % (npass, len(res), ' (without the MODFLOW library)' if NOLIB else '', os.path.join(W, 'results.json')))
