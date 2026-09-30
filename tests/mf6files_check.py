"""Checks the MODFLOW 6 files Raven writes for every grid type - no MODFLOW 6 library needed.

Raven builds the whole groundwater model (grid, layers, links, boundaries, all MODFLOW input files and
grid_cells.geojson) before it loads the MODFLOW library. Run with a library path that does not exist, Raven stops right
there, and the files can be checked on their own:

  - every file loads in FloPy;
  - cells: valid polygons that tile the active area (grid_cells.geojson);
  - rotated regular grids: the DIS origin and rotation put every cell where grid_cells.geojson puts it;
  - polygon boundaries select exactly the cells whose centre lies in the polygon;
  - DISU: connections symmetric; each face direction (ANGLDEGX) agrees with the cell vertices within 1 degree (the check
    MODFLOW itself applies); vertical connections bridge pass-through cells only, never an aquiclude;
  - drains: the smoothing interval stays inside the top cell (MODFLOW refuses it otherwise);
  - HRU meshes: share of each active cell's area in its dominant HRU, computed independently with shapely.

Usage:  python3 mf6files_check.py [OLD_RAVEN_EXE]
With OLD_RAVEN_EXE (an earlier build) the written files of both builds are compared case by case, to show which
results can have changed between the two builds. Results: <scratch>/mf6files/results.json
"""
import os, sys, re, json, shutil, subprocess, math, filecmp
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import testconfig as tc
import flopy, pyproj
from shapely.geometry import shape, Polygon, Point
from shapely.ops import unary_union

NOLIB = '/nonexistent/libmf6.so'
OUT = os.path.dirname(tc.work('mf6files', 'x'))
OLD = sys.argv[1] if len(sys.argv) > 1 else None
EX = tc.LIARD


def rvg_uncomment(g, key):
    """uncomment the example's line(s) starting with '# key'"""
    return re.sub(r'(?m)^# (' + re.escape(key) + ')', r'\1', g)


def build(name, extra='', edits=(), base=EX, exe=None, tag='new', model='Liard'):
    """copies the base case, applies edits (file, old, new) and extra .rvg lines, runs Raven up to the library load"""
    d = os.path.join(OUT, tag, name)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(base, d, ignore=shutil.ignore_patterns('out', 'output', '__pycache__', '*.gwcache'))
    for f, a, b in edits:
        p = os.path.join(d, f); s = open(p).read()
        if callable(a):
            s = a(s)
        else:
            assert a in s, (name, f, a[:60]); s = s.replace(a, b, 1)
        open(p, 'w').write(s)
    if extra:
        open(os.path.join(d, model + '.rvg'), 'a').write('\n' + extra)
    os.makedirs(os.path.join(d, 'out'), exist_ok=True)
    env = dict(os.environ, RAVEN_MF6_LIB=NOLIB)
    r = subprocess.run([exe or tc.RAVEN, model, '-o', 'out/'], cwd=d, capture_output=True, text=True, env=env, timeout=1800)
    m = re.findall(r'Exiting Gracefully: (.*)', r.stdout + r.stderr)
    msg = m[-1] if m else '(no exit message; rc %d)' % r.returncode
    return d, msg, r.returncode


def laea_of(d):
    head = open([os.path.join(d, 'out', 'mf6', f) for f in os.listdir(os.path.join(d, 'out', 'mf6')) if re.match(r'gwf\.dis', f)][0]).readline()
    lat0, lon0 = map(float, re.findall(r'lat=([-0-9.]+) lon=([-0-9.]+)', head)[0])
    return pyproj.Proj(proj='laea', lat_0=lat0, lon_0=lon0, R=6371007.181)


def cells_xy(d, P):
    g = json.load(open(os.path.join(d, 'out', 'mf6', 'grid_cells.geojson')))['features']
    polys = {}
    for f in g:
        ring = np.array(f['geometry']['coordinates'][0])[:, :2]
        x, y = P(ring[:, 0], ring[:, 1]); polys[f['properties']['cell']] = (Polygon(np.c_[x, y]), f['properties'])
    return polys


def check_case(name, d, msg, rc, want=('cannot load MODFLOW 6 library',), aquiclude=False, chd_polygon=None, mesh_fidelity=False):
    R = dict(name=name, msg=msg[:160], checks=[])
    def chk(label, ok, value=''):
        R['checks'].append((label, bool(ok), str(value)))
    reached = any(w in msg for w in want)
    chk('build reached the MODFLOW library load', reached, msg[:90])
    if not reached or rc in (-11, -6, 134, 139):
        return R
    ws = os.path.join(d, 'out', 'mf6')
    try:
        sim = flopy.mf6.MFSimulation.load(sim_ws=ws, verbosity_level=0)
        gwf = sim.get_model(sim.model_names[0]); chk('FloPy loads every file', True)
    except Exception as e:
        chk('FloPy loads every file', False, repr(e)[:120]); return R
    P = laea_of(d); polys = cells_xy(d, P)
    bad = [c for c, (p, pr) in polys.items() if not p.is_valid or p.area <= 0 or len(p.exterior.coords) < 4]
    chk('grid_cells.geojson: every cell a valid polygon', not bad, '%d bad of %d' % (len(bad), len(polys)))
    act = [p for p, pr in polys.values() if pr.get('active')]
    tot = sum(p.area for p in act); uni = unary_union(act).area
    chk('active cells do not overlap (area of union / sum)', abs(uni / tot - 1) < 1e-6, '%.2e' % (uni / tot - 1))
    dis = gwf.get_package('dis') if 'dis' in gwf.package_type_dict else None
    disu = gwf.get_package('disu') if 'disu' in gwf.package_type_dict else None
    disv = gwf.get_package('disv') if 'disv' in gwf.package_type_dict else None
    # ---- rotated/refined regular grids: DIS origin and rotation against the cell polygons
    if dis is not None:
        delr = dis.delr.array; delc = dis.delc.array; nrow, ncol = len(delc), len(delr)
        x0 = dis.xorigin.get_data() or 0.0; y0 = dis.yorigin.get_data() or 0.0; a = math.radians(dis.angrot.get_data() or 0.0)
        xe = np.r_[0, np.cumsum(delr)]; ye_from_bottom = np.r_[0, np.cumsum(delc[::-1])]
        worst = 0.0
        for c_id, (p, pr) in polys.items():
            r, c = divmod(c_id - 1, ncol)
            xl = 0.5 * (xe[c] + xe[c + 1]); yl = 0.5 * (ye_from_bottom[nrow - 1 - r] + ye_from_bottom[nrow - r])
            X = x0 + xl * math.cos(a) - yl * math.sin(a); Y = y0 + xl * math.sin(a) + yl * math.cos(a)
            worst = max(worst, math.hypot(X - p.centroid.x, Y - p.centroid.y))
        chk('DIS origin/rotation place every cell where grid_cells.geojson does (worst, m)', worst < 0.5, '%.3f' % worst)
    # ---- DISV/DISU: vertices, origin and rotation against the cell polygons
    vp = disv if disv is not None else disu
    if vp is not None and vp.vertices.has_data():
        x0 = vp.xorigin.get_data() or 0.0; y0 = vp.yorigin.get_data() or 0.0; a = math.radians(vp.angrot.get_data() or 0.0)
        V = {int(v[0]): (v[1], v[2]) for v in vp.vertices.array}
        worst = 0.0; ncpl = len(polys)
        for c in vp.cell2d.array:
            icell = int(c[0])
            if icell >= ncpl: break
            pts = [V[int(k)] for k in list(c)[4:4 + int(c[3])]]
            if len(pts) < 3: worst = 1e9; continue
            pg = Polygon([(x0 + px * math.cos(a) - py * math.sin(a), y0 + px * math.sin(a) + py * math.cos(a)) for px, py in pts])
            q = polys[icell + 1][0]
            worst = max(worst, math.hypot(pg.centroid.x - q.centroid.x, pg.centroid.y - q.centroid.y))
        chk('vertices, origin and rotation place every cell where grid_cells.geojson does (worst, m)', worst < 0.5, '%.3f' % worst)
    # ---- polygon boundary: cells whose centre lies inside
    if chd_polygon is not None:
        chd = gwf.get_package('chd')
        rec = chd.stress_period_data.get_data()[0]
        got = set()
        for cid in rec['cellid']:
            if dis is not None: got.add(cid[1] * dis.delr.array.size + cid[2] + 1)
            elif disv is not None: got.add(cid[1] + 1)
            else:
                ncpl = len(polys); got.add(cid[0] % ncpl + 1)
        poly = unary_union([shape(f['geometry']) for f in json.load(open(os.path.join(d, chd_polygon)))['features']])
        ext = np.array(poly.exterior.coords); x, y = P(ext[:, 0], ext[:, 1]); pl = Polygon(np.c_[x, y])
        exp = set(c for c, (p, pr) in polys.items() if pr.get('active') and pl.contains(p.centroid))
        chk('polygon boundary takes exactly the cells whose centre lies inside', got == exp and len(exp) > 0,
            'selected %d, expected %d, differing %d' % (len(got), len(exp), len(got ^ exp)))
    # ---- DISU: symmetry, face directions, aquicludes
    if disu is not None:
        iac = disu.iac.array; ja = disu.ja.array; ihc = disu.ihc.array; cl12 = disu.cl12.array; hwva = disu.hwva.array
        ang = disu.angldegx.array if disu.angldegx.has_data() else None
        idom = disu.idomain.array if disu.idomain.has_data() else np.ones(len(iac), int)
        nodes = len(iac); ia = np.r_[0, np.cumsum(iac)]
        conn = {}
        for n in range(nodes):
            for q in range(ia[n] + 1, ia[n + 1]):
                conn[(n, int(ja[q]))] = (ihc[q], hwva[q], q)   # (FloPy gives JA zero-based)
        asym = sum(1 for (n, m), (h, w, q) in conn.items() if (m, n) not in conn or conn[(m, n)][0] != h or abs(conn[(m, n)][1] - w) > 1e-6 * max(abs(w), 1))
        chk('DISU connections symmetric (IHC, HWVA)', asym == 0, '%d asymmetric of %d' % (asym, len(conn)))
        ncpl = len(polys)
        if ang is not None:
            verts = {int(v[0]): (v[1], v[2]) for v in disu.vertices.array}
            c2d = {int(c[0]): c for c in disu.cell2d.array}
            worst = 0.0; nbad = 0; ntest = 0
            for (n, m), (h, w, q) in conn.items():
                if h == 0 or n >= m: continue
                cn, cm = c2d[n], c2d[m]
                vn = set(int(v) for v in list(cn)[4:4 + int(cn[3])]); vm = set(int(v) for v in list(cm)[4:4 + int(cm[3])])
                sh = list(vn & vm)
                if len(sh) < 2: continue
                (xa, ya), (xb, yb) = verts[sh[0]], verts[sh[1]]
                nx, ny = yb - ya, -(xb - xa)
                if (cm[1] - cn[1]) * nx + (cm[2] - cn[2]) * ny < 0: nx, ny = -nx, -ny
                a_geo = math.degrees(math.atan2(ny, nx)) % 360
                dlt = abs((ang[q] - a_geo + 180) % 360 - 180); ntest += 1
                worst = max(worst, dlt); nbad += dlt > 1.0
            chk('DISU face directions agree with the vertices (within 1 degree)', nbad == 0,
                'worst %.3f deg over %d faces' % (worst, ntest))
        bridged_aquiclude = 0; nvert = 0
        for (n, m), (h, w, q) in conn.items():
            if h != 0 or m <= n: continue
            nvert += 1
            for k in range(n + ncpl, m, ncpl):
                if idom[k] == 0: bridged_aquiclude += 1
        chk('vertical connections never cross an aquiclude (IDOMAIN 0)', bridged_aquiclude == 0,
            '%d crossing of %d vertical connections' % (bridged_aquiclude, nvert))
        if aquiclude:
            chk('the case has aquiclude cells', int((idom == 0).sum()) > 0, int((idom == 0).sum()))
    # ---- drains inside the top cell
    drn = gwf.get_package('drn') if 'drn' in gwf.package_type_dict else None
    if drn is not None:
        rec = drn.stress_period_data.get_data()[0]
        if dis is not None: bot = dis.botm.array[0].ravel(); idx = [c[1] * dis.delr.array.size + c[2] for c in rec['cellid']]
        elif disv is not None: bot = disv.botm.array[0]; idx = [c[1] for c in rec['cellid']]
        else: bot = disu.bot.array; idx = [c[0] for c in rec['cellid']]
        aux = rec['ddrn'] if 'ddrn' in rec.dtype.names else np.zeros(len(rec))
        below = sum(1 for e, dd, i in zip(rec['elev'], aux, idx) if e + dd < bot[i] - 1e-9)
        chk('drain smoothing interval inside the top cell', below == 0, '%d of %d below the cell bottom' % (below, len(rec)))
    # ---- HRU-mesh fidelity (independent)
    if mesh_fidelity:
        H = json.load(open(os.path.join(d, 'Liard_HRUs_standin.geojson')))['features']
        hp = []
        for f in H:
            gg = shape(f['geometry'])
            from shapely.ops import transform
            hp.append(transform(lambda x, y, z=None: P(x, y), gg).buffer(0))
        from shapely.strtree import STRtree
        tree = STRtree(hp); dom = []
        for p, pr in polys.values():
            if not pr.get('active'): continue
            best = 0.0
            for j in tree.query(p):
                best = max(best, p.intersection(hp[j]).area)
            dom.append(best / p.area)
        R['mesh_dominant_share'] = float(np.average(dom, weights=[p.area for p, pr in polys.values() if pr.get('active')]))
        chk('HRU mesh: area-weighted share of each active cell in its dominant HRU', R['mesh_dominant_share'] > 0.95, '%.4f' % R['mesh_dominant_share'])
    return R


# ------------------------------------------------------------------------------------------------ cases
aq = ('Liard.rvp', 'CLAY_TILL, 5.0, AQUITARD', 'CLAY_TILL, 5.0, AQUICLUDE')
res_blk = None
def lake_edits():
    """a reservoir with seepage on the largest coupled HRU of subbasin 52 (as in lake_test.py)"""
    rvh = open(os.path.join(EX, 'Liard.rvh')).read()
    hdr = [l for l in rvh.split('\n') if ':Attributes' in l and 'AQUIFER_PROFILE' in l][0].replace(',', ' ').split()[1:]
    ib, ia, ie, iar = hdr.index('BASIN_ID') + 1, hdr.index('AQUIFER_PROFILE') + 1, hdr.index('ELEVATION') + 1, hdr.index('AREA') + 1
    rows = [l.replace(',', ' ').split() for l in rvh.split(':HRUs')[1].split(':EndHRUs')[0].split('\n')]
    cand = [r for r in rows if r and r[0].isdigit() and r[ib] == '52' and r[ia] != '[NONE]']
    lake = max(cand, key=lambda r: float(r[iar])); lid = int(lake[0]); el = float(lake[ie])
    blk = ("\n:Reservoir TestLake\n  :SubBasinID 52\n  :HRUID %d\n  :Type RESROUTE_STANDARD\n  :WeirCoefficient 0.6\n  :CrestWidth 20.0\n"
           "  :MaxDepth 8.0\n  :LakeArea 2.0E7\n  :AbsoluteCrestHeight %.2f\n  :SeepageParameters 0.05 %.2f\n:EndReservoir\n" % (lid, el, el - 5))
    return [('Liard.rvh', lambda s: s + blk, None)]

boundaries = [('Liard.rvg', lambda s: rvg_uncomment(rvg_uncomment(s, ':SpecifiedHeadBoundary'), ':GeneralHeadBoundary'), None)]
obs_out = ':ObservationWells\n  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_ELEV\n  9, OW_FAR, 58.2, -127.9, 1000\n:EndObservationWells\n'
thin = ('Liard.rvp', 'GRAVEL, 15.0, AQUIFER', 'GRAVEL, 0.6, AQUIFER')
smooth = ('Liard.rvp', lambda s: s.replace('SEEPAGE_LEAKANCE, SOIL_ZONE_DEPTH', 'SEEPAGE_LEAKANCE, SOIL_ZONE_DEPTH, SEEPAGE_SMOOTHING_DEPTH')
          .replace('1/d,              m\n', '1/d,              m, m\n').replace('1.0,              0.0\n', '1.0,              0.0, 1.0\n'), None)

CASES = [
    ('regular (default)', '', [], {}),
    ('rotated 20, rivers 500', ':GridRotation 20\n:GridRefinement RIVERS 500\n', [], {}),
    ('quadtree, rivers 500', ':GridType QUADTREE\n:GridRefinement RIVERS 500\n', [], {}),
    ('HRU mesh, rivers 1000', ':GridType HRU_MESH\n:GridRefinement RIVERS 1000\n', [], dict(mesh_fidelity=True)),
    ('layers by elevation', ':LayerConnection ELEVATION\n', [], {}),
    ('rotated 20 uniform, boundaries', ':GridRotation 20\n', boundaries, dict(chd_polygon='lake.geojson')),
    ('regular, boundaries', '', boundaries, dict(chd_polygon='lake.geojson')),
    ('HRU mesh, aquiclude', ':GridType HRU_MESH\n', [aq], dict(aquiclude=True, mesh_fidelity=True)),
    ('elevation, aquiclude', ':LayerConnection ELEVATION\n', [aq], dict(aquiclude=True)),
    ('rotated 30, elevation', ':GridRotation 30\n:LayerConnection ELEVATION\n', [], {}),
    ('thin top layer, smoothing 1 m', '', [thin, smooth], {}),
    ('HRU mesh, wells, monitoring well outside', ':GridType HRU_MESH\n:GridRefinement WELLS 400\n' + obs_out, [], {}),
    ('HRU mesh, lakes 500', ':ReservoirExchange\n:GridType HRU_MESH\n:GridRefinement LAKES 500\n', lake_edits(), dict(mesh_fidelity=True)),
    ('quadtree, lakes 500', ':ReservoirExchange\n:GridType QUADTREE\n:GridRefinement LAKES 500\n', lake_edits(), {}),
    ('rotated 20, lakes 750', ':ReservoirExchange\n:GridRotation 20\n:GridRefinement LAKES 750\n', lake_edits(), {}),
]


def main():
    results = []
    qt = None
    for name, extra, edits, kw in CASES:
        d, msg, rc = build(re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_'), extra, edits)
        R = check_case(name, d, msg, rc, **kw)
        if name.startswith('quadtree, rivers'):
            qt = os.path.join(d, 'out', 'mf6', 'grid_cells.geojson')
        if OLD:
            do, mo, ro = build(re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_'), extra, edits, exe=OLD, tag='old')
            fs = sorted(f for f in os.listdir(os.path.join(d, 'out', 'mf6')) if f != 'linkage.cache')
            diff = [f for f in fs if not os.path.exists(os.path.join(do, 'out', 'mf6', f)) or not filecmp.cmp(os.path.join(d, 'out', 'mf6', f), os.path.join(do, 'out', 'mf6', f), shallow=False)]
            R['files_differ_from_old_build'] = diff
            if kw.get('mesh_fidelity') and 'cannot load' in mo:
                Ro = check_case(name + ' (old build)', do, mo, ro, mesh_fidelity=True)
                R['old_mesh_dominant_share'] = Ro.get('mesh_dominant_share')
        results.append(R)
        print('%-44s %s' % (name, 'PASS' if all(c[1] for c in R['checks']) else 'FAIL'),
              ('files differ from old build: %s' % (R['files_differ_from_old_build'] or 'none')) if OLD else '', flush=True)
        for c in R['checks']:
            if not c[1]: print('      FAIL %s: %s' % (c[0], c[2]))
    if qt:   # imported cells: the quadtree case's cells
        name = 'imported cells (FILE)'
        d, msg, rc = build('imported_cells', ':GridType FILE\n:GridFile %s cell\n' % qt, [])
        R = check_case(name, d, msg, rc); results.append(R)
        print('%-44s %s' % (name, 'PASS' if all(c[1] for c in R['checks']) else 'FAIL'))
        for c in R['checks']:
            if not c[1]: print('      FAIL %s: %s' % (c[0], c[2]))
    json.dump(results, open(os.path.join(OUT, 'results.json'), 'w'), indent=1)
    n = sum(1 for R in results for c in R['checks']); nf = sum(1 for R in results for c in R['checks'] if not c[1])
    print('%d cases, %d checks, %d failed; results in %s' % (len(results), n, nf, os.path.join(OUT, 'results.json')))


if __name__ == '__main__':
    main()
