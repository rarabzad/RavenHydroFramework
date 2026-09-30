"""Checks the geometry library (GWGeometry) against independent implementations.

Build and run from this folder:
    g++ -std=c++11 -O2 -I.. geomtest.cpp ../GWGeometry.cpp -o geomtest
    python3 geomtest_check.py            (needs numpy, shapely, pyproj)

geomtest writes clipping areas, projected points and raster samples; this script writes the test raster first,
runs geomtest, then compares:
  - clipping: area of 400 random (often concave) polygons inside random squares, against shapely;
  - projection: Lambert azimuthal equal-area (sphere R=6371007.181 m) against PROJ, and the inverse round trip;
  - raster: bilinear sampling of an ESRI ASCII grid holding a plane, which bilinear sampling must reproduce exactly.
"""
import os, subprocess, sys
import numpy as np
from shapely.geometry import Polygon, box
import pyproj

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)

# test raster: a plane v = 3*lon + 2*lat + 1 on a 0.5-degree grid (lon -3..3, lat 9..14); cell centres are sampled
ncols, nrows, x0, y0, cs = 12, 10, -3.0, 9.0, 0.5
with open('test.asc', 'w') as f:
    f.write(f'ncols {ncols}\nnrows {nrows}\nxllcorner {x0}\nyllcorner {y0}\ncellsize {cs}\nNODATA_value -9999\n')
    for r in range(nrows):
        lat = y0 + (nrows - 1 - r + 0.5) * cs
        f.write(' '.join('%.6f' % (3 * (x0 + (c + 0.5) * cs) + 2 * lat + 1) for c in range(ncols)) + '\n')

exe = os.path.join(here, 'geomtest')
if not os.path.exists(exe):
    sys.exit('build geomtest first (see the docstring)')
out = subprocess.run([exe], capture_output=True, text=True)
if 'done' not in out.stdout:
    sys.exit('geomtest failed: ' + out.stdout + out.stderr)

fails = 0
def report(name, ok, value):
    global fails
    print('%-60s %s  (%.3g)' % (name, 'ok' if ok else 'FAIL', value))
    fails += 0 if ok else 1

# clipping against shapely
worst = 0.0
for line in open('clip.txt'):
    v = line.split(); n = int(v[0]); xa, ya, xb, yb, A = map(float, v[1:6])
    pts = np.array(v[6:6 + 2 * n], float).reshape(-1, 2)
    P = Polygon(pts)
    if not P.is_valid:
        P = P.buffer(0)
    ref = P.intersection(box(xa, ya, xb, yb)).area
    worst = max(worst, abs(A - ref) / max(ref, 1e-9) if ref > 1e-9 else abs(A))
report('clipping: 400 polygons against shapely (worst relative error)', worst < 1e-9, worst)

# projection against PROJ
d = np.loadtxt('proj.txt')
laea = pyproj.Proj(proj='laea', lat_0=59.0, lon_0=-125.3, R=6371007.181)
x, y = laea(d[:, 0], d[:, 1])
e = np.max(np.hypot(x - d[:, 2], y - d[:, 3]))
report('projection: forward against PROJ (worst position error, m)', e < 1e-3, e)
rt = np.max(np.abs(np.r_[d[:, 4] - d[:, 0], d[:, 5] - d[:, 1]]))
report('projection: inverse round trip (worst error, degrees)', rt < 1e-9, rt)

# raster: a plane is reproduced exactly by bilinear sampling (points inside the outer cell centres)
r = np.loadtxt('raster.txt')
inside = (r[:, 0] >= x0 + cs / 2) & (r[:, 0] <= x0 + (ncols - 0.5) * cs) & (r[:, 1] >= y0 + cs / 2) & (r[:, 1] <= y0 + (nrows - 0.5) * cs)
ok = r[:, 2] == 1
e = np.max(np.abs(r[inside & ok, 3] - (3 * r[inside & ok, 0] + 2 * r[inside & ok, 1] + 1))) if np.any(inside & ok) else 1e9
report('raster: bilinear sampling of a plane (worst error)', e < 1e-6 and np.all(ok[inside]), e)

for f in ('clip.txt', 'proj.txt', 'raster.txt', 'test.asc'):
    if os.path.exists(f):   # another run in the same checkout may have removed it already
        os.remove(f)
print('ALL PASSED' if fails == 0 else '%d FAILED' % fails)
sys.exit(1 if fails else 0)
