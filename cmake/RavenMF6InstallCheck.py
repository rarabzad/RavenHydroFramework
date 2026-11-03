"""Install check of Raven with the MODFLOW 6 library: runs the installed Raven on a tiny coupled model and checks that it
used the library beside it, with the Fortran runtime installed beside it (cmake/RavenMF6Bundle.cmake), not the one of the
compiler that built it.

    python RavenMF6InstallCheck.py <install>/bin [expected MODFLOW version, e.g. 6.8.1]

The model (written to a temporary folder): one 200 m x 200 m HRU over a flat 4-cell aquifer, two daily steps, rain into
the soil and percolation to GROUNDWATER (:SoilZone CONCEPTUAL). Raven runs with RAVEN_MF6_LIB unset, so it takes the
library next to itself. Checks:
  - the run ends with "Successful Simulation" and GWBudget.csv has two steps;
  - GWModelSummary.txt names <bin>/libmf6.* (and the expected version, if given);
  - every Fortran runtime library the library loads comes from <bin> (none for the official release libraries; a
    source-built one must have its runtime installed beside it, not taken from the compiler's folder):
      Linux    ldd of <bin>/libmf6.so: every libgfortran/libquadmath it resolves lies in <bin> (libgcc_s is a system
               library there);
      macOS    DYLD_PRINT_LIBRARIES of the Raven run: every libgfortran/libquadmath/libgcc_s image loaded lies in <bin>;
      Windows  the run's PATH has every folder holding libgfortran-*.dll removed except <bin> (Raven loads an absolute
               library path with LOAD_WITH_ALTERED_SEARCH_PATH: the DLLs beside it are found first).
Prints one line per check and "install check: PASS" or "install check: FAIL"; exit code 0 or 1.
"""
import json, math, os, platform, re, shutil, subprocess, sys, tempfile

RUNTIME = re.compile(r'(lib)?(gfortran|quadmath|gcc_s)[^/\\\s]*\.(so[.0-9]*|[0-9.]*dylib|dll)', re.I)


def ll(x, y, lat0=50.0, lon0=-100.0):
    """lon, lat of x, y [m] (inverse Lambert azimuthal equal area about lat0, lon0, as Raven projects HRUs)"""
    R = 6371007.181; p0, l0 = math.radians(lat0), math.radians(lon0); rho = math.hypot(x, y)
    if rho < 1e-9:
        return lon0, lat0
    c = 2.0 * math.asin(rho / (2.0 * R))
    lat = math.asin(math.cos(c) * math.sin(p0) + y * math.sin(c) * math.cos(p0) / rho)
    lon = l0 + math.atan2(x * math.sin(c), rho * math.cos(p0) * math.cos(c) - y * math.sin(p0) * math.sin(c))
    return math.degrees(lon), math.degrees(lat)


def write_model(d):
    ring = [ll(-100, -100), ll(100, -100), ll(100, 100), ll(-100, 100), ll(-100, -100)]
    json.dump({"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {"HRU_ID": 1},
               "geometry": {"type": "Polygon", "coordinates": [ring]}}]}, open(os.path.join(d, 'hru.geojson'), 'w'))
    lonA, latA = ll(-300, -300); lonB, latB = ll(300, 300); n = 20; cs = (lonB - lonA) / n
    nr = int(math.ceil((latB - latA) / cs)) + 1
    with open(os.path.join(d, 'dem.asc'), 'w') as f:
        f.write('ncols %d\nnrows %d\nxllcorner %.10f\nyllcorner %.10f\ncellsize %.12f\nNODATA_value -9999\n' % (n, nr, lonA, latA, cs))
        f.write(('100.0 ' * n + '\n') * nr)
    lon, lat = ll(0, 0)
    open(os.path.join(d, 'm.rvi'), 'w').write(
        ':StartDate 2000-01-01 00:00:00\n:Duration 2\n:TimeStep 1.0\n:Method ORDERED_SERIES\n:SoilModel SOIL_MULTILAYER 2\n'
        ':Routing ROUTE_NONE\n:CatchmentRoute ROUTE_DUMP\n:Evaporation PET_CONSTANT\n:RainSnowFraction RAINSNOW_DATA\n'
        ':OroTempCorrect OROCORR_NONE\n:OroPrecipCorrect OROCORR_NONE\n:HydrologicProcesses\n'
        '  :Precipitation PRECIP_RAVEN ATMOS_PRECIP MULTIPLE\n  :Infiltration INF_RATIONAL PONDED_WATER MULTIPLE\n'
        '  :Percolation PERC_CONSTANT SOIL[0] SOIL[1]\n  :Percolation PERC_CONSTANT SOIL[1] GROUNDWATER\n'
        ':EndHydrologicProcesses\n:GroundwaterModel MODFLOW6\n:SilentMode\n')
    open(os.path.join(d, 'm.rvh'), 'w').write(
        ':SubBasins\n  :Attributes, NAME, DOWNSTREAM_ID, PROFILE, REACH_LENGTH, GAUGED\n  :Units, none, none, none, km, none\n'
        '  1, box, -1, NONE, ZERO-, 1\n:EndSubBasins\n:HRUs\n  :Attributes, AREA, ELEVATION, LATITUDE, LONGITUDE, BASIN_ID, '
        'LAND_USE_CLASS, VEG_CLASS, SOIL_PROFILE, AQUIFER_PROFILE, TERRAIN_CLASS, SLOPE, ASPECT\n'
        '  :Units, km2, m, deg, deg, none, none, none, none, none, none, deg, deg\n'
        '  1, 0.04, 100.0, %.6f, %.6f, 1, OPEN, GRASS, SOILP, AQ, [NONE], 0, 0\n:EndHRUs\n' % (lat, lon))
    open(os.path.join(d, 'm.rvp'), 'w').write(
        ':SoilClasses\n  :Attributes\n  :Units\n  SOIL\n:EndSoilClasses\n:SoilProfiles\n  SOILP, 2, SOIL, 0.5, SOIL, 0.5\n'
        ':EndSoilProfiles\n:VegetationClasses\n  :Attributes, MAX_HT, MAX_LAI, MAX_LEAF_COND\n  :Units, m, none, mm_per_s\n'
        '  GRASS, 1, 1, 5\n:EndVegetationClasses\n:LandUseClasses\n  :Attributes, IMPERM, FOREST_COV\n  :Units, frac, frac\n'
        '  OPEN, 0, 0\n:EndLandUseClasses\n:SoilParameterList\n  :Parameters, POROSITY, FIELD_CAPACITY, SAT_WILT, HYDRAUL_COND, '
        'MAX_PERC_RATE\n  :Units, none, none, none, mm_per_d, mm_per_d\n  [DEFAULT], 0.4, 0.5, 0.0, 100, 5\n'
        ':EndSoilParameterList\n:LandUseParameterList\n  :Parameters, PARTITION_COEFF\n  :Units, none\n  [DEFAULT], 0.0\n'
        ':EndLandUseParameterList\n:AquiferClasses\n  :Attributes, K_HORIZ, K_VERT, SPEC_STORAGE, SPEC_YIELD, POROSITY\n'
        '  :Units, m/d, m/d, 1/m, none, none\n  SAND, 10, 1, 1.0E-5, 0.2, 0.3\n:EndAquiferClasses\n:AquiferProfiles\n'
        '  AQ, 1, SAND, TO_BEDROCK, AQUIFER\n:EndAquiferProfiles\n:AquiferProfileParameters\n'
        '  :Attributes, INITIAL_HEAD_DEPTH, DEFAULT_BEDROCK_DEPTH\n  :Units, m, m\n  AQ, 2.0, 20.0\n:EndAquiferProfileParameters\n')
    data = lambda v: '    2000-01-01 00:00:00 1.0 2\n    %s\n    %s\n  :EndData\n' % (v, v)
    open(os.path.join(d, 'm.rvt'), 'w').write(
        ':Gauge G\n  :Latitude %.6f\n  :Longitude %.6f\n  :Elevation 100\n  :Data RAINFALL mm/d\n%s  :Data SNOWFALL mm/d\n%s'
        '  :Data TEMP_AVE C\n%s:EndGauge\n' % (lat, lon, data(10), data(0), data(10)))
    open(os.path.join(d, 'm.rvc'), 'w').write(':UniformInitialConditions SOIL[0] 100\n:UniformInitialConditions SOIL[1] 100\n')
    open(os.path.join(d, 'm.rvg'), 'w').write(':HRUGeometry hru.geojson HRU_ID\n:LandSurface dem.asc\n:GridCellSize 100\n')


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    bindir = os.path.abspath(sys.argv[1]); want = sys.argv[2] if len(sys.argv) > 2 else ''
    sysname = platform.system()
    exe = os.path.join(bindir, 'Raven.exe' if sysname == 'Windows' else 'Raven')
    lib = os.path.join(bindir, {'Windows': 'libmf6.dll', 'Darwin': 'libmf6.dylib'}.get(sysname, 'libmf6.so'))
    ok = True

    def check(label, good, detail=''):
        nonlocal ok
        ok &= bool(good)
        print('%-4s %s%s' % ('ok' if good else 'FAIL', label, ('  (%s)' % detail) if detail else ''), flush=True)

    check('Raven and the library in %s' % bindir, os.path.isfile(exe) and os.path.isfile(lib))
    if not ok:
        print('install check: FAIL'); return 1
    bundled = sorted(f for f in os.listdir(bindir) if RUNTIME.fullmatch(f))
    print('     runtime libraries beside the library: %s' % (', '.join(bundled) or 'none'))
    env = {k: v for k, v in os.environ.items() if k not in ('RAVEN_MF6_LIB', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH',
                                                              'DYLD_FALLBACK_LIBRARY_PATH')}
    if sysname == 'Darwin':
        env['DYLD_PRINT_LIBRARIES'] = '1'
    elif sysname == 'Windows':            # hide the compiler's runtime: only <bin> may supply it
        keep, hidden = [], []
        for p in env.get('PATH', '').split(os.pathsep):
            if p and os.path.isdir(p) and os.path.normcase(os.path.abspath(p)) != os.path.normcase(bindir) and \
                    any(re.match(r'libgfortran-\d+\.dll$', f, re.I) for f in os.listdir(p)):
                hidden.append(p)
            else:
                keep.append(p)
        env['PATH'] = os.pathsep.join(keep)
        print('     PATH folders hidden from the run (they hold libgfortran): %s' % ('; '.join(hidden) or 'none'))
    d = tempfile.mkdtemp(prefix='raven_mf6_check_')
    try:
        write_model(d); os.makedirs(os.path.join(d, 'out'))
        r = subprocess.run([exe, 'm', '-o', 'out' + os.sep], cwd=d, env=env, capture_output=True, text=True, timeout=600)
        out = r.stdout + r.stderr
        msg = re.findall(r'Exiting Gracefully: (.*)', out)
        check('coupled run finishes', 'Successful' in ''.join(msg), (msg[-1] if msg else 'exit code %d: %s' % (r.returncode, out[-400:])))
        b = os.path.join(d, 'out', 'GWBudget.csv')
        rows = len(open(b).read().strip().split('\n')) - 1 if os.path.exists(b) else 0
        check('GWBudget.csv has the two steps', rows == 2, '%d rows' % rows)
        s = os.path.join(d, 'out', 'GWModelSummary.txt')
        line = ([l for l in open(s).read().split('\n') if 'MODFLOW 6 library' in l] or [''])[0] if os.path.exists(s) else ''
        m = re.search(r':\s*(.*?)(\s+\(version ([^)]*)\))?\s*$', line)
        used = m.group(1) if m else ''
        same = bool(used) and os.path.exists(used) and os.path.samefile(used, lib)
        check('the library beside Raven was used', same, line.strip())
        if want:
            v = (m.group(3) or '') if m else ''
            check('MODFLOW %s' % want, v.startswith(want), v or 'no version')
        if sysname == 'Linux':
            ld = subprocess.run(['ldd', lib], capture_output=True, text=True).stdout
            res = re.findall(r'^\s*(\S+)\s+=>\s+(\S+)', ld, re.M)
            rt = [(n, p) for n, p in res if RUNTIME.fullmatch(n) and not n.startswith('libgcc_s')]
            bad = ['%s => %s' % (n, p) for n, p in rt if os.path.dirname(os.path.abspath(p)) != bindir]
            check('the Fortran runtime resolves to %s' % bindir, not bad, '; '.join(bad) or ('%d libraries' % len(rt)))
        elif sysname == 'Darwin':
            imgs = re.findall(r'dyld(?:\[\d+\])?: <[^>]*> (\S+)|dyld: loaded: (?:<[^>]*> )?(\S+)', out)
            imgs = [a or b for a, b in imgs]
            rt = [p for p in imgs if RUNTIME.search(os.path.basename(p))]
            bad = [p for p in rt if os.path.realpath(os.path.dirname(p)) != os.path.realpath(bindir)]
            check('the Fortran runtime loaded from %s' % bindir, not bad, '; '.join(bad) or ('%d images' % len(rt)))
    finally:
        shutil.rmtree(d, ignore_errors=True)
    print('install check: %s' % ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
