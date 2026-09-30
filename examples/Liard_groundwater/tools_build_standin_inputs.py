# Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
# SPDX-License-Identifier: Artistic-2.0
"""Builds the Liard groundwater test case for the Raven-MODFLOW 6 coupling.
The Liard package has no HRU polygons, DEM or bedrock data, so stand-ins are generated
here (FOR TESTING ONLY): Voronoi HRU polygons from .rvh centroids clipped to subbasins,
a DEM interpolated from HRU elevations, and bedrock = DEM - depth by land type."""
import json, re, sys, numpy as np
# Run this in a fresh copy of the ORIGINAL Liard package (after converting Windows paths and renaming
# DEP_THRESHHOLD -> DEP_THRESHOLD for Raven 4.x). It edits Liard.rvh/.rvi/.rvp and writes Liard.rvg.
if ':AquiferClasses' in open('Liard.rvp').read():
    sys.exit('Liard.rvp already contains groundwater definitions: run this script on the original package files only')
from scipy.spatial import Voronoi
import shapely.geometry as sg
from shapely.ops import unary_union

SUBBASINS=[43,48,52]
VALLEY=('FOREST','SHRUBLAND','GRASSLAND','WETLAND','CROPLAND')
rvh=open('Liard.rvh').read()
head,rest=rvh.split(':HRUs',1); body,tail=rest.split(':EndHRUs',1)
lines=body.split('\n'); hrus=[]
for i,l in enumerate(lines):
    c=[x.strip() for x in l.split(',')]
    if len(c)>=13 and not l.strip().startswith(':') and c[0].isdigit():
        hrus.append(dict(i=i,id=int(c[0]),area=float(c[1]),elev=float(c[2]),lat=float(c[3]),lon=float(c[4]),sb=int(c[5]),lu=c[6],soil=c[8],cols=c))
sel=[h for h in hrus if h['sb'] in SUBBASINS]
def profile(h):
    if h['lu']=='GLACIER' or h['soil']=='GLACIER': return '[NONE]'
    return 'VALLEY_ALLUVIUM' if h['lu'] in VALLEY else 'UPLAND_TILL'
for h in sel: h['prof']=profile(h)

# ---- 1. stand-in HRU polygons: Voronoi of centroids (in locally scaled lon/lat) clipped to subbasin
sbgeo={f['properties']['SubId']:sg.shape(f['geometry']) for f in json.load(open('Liard_subbasins.geojson'))['features']}
lat0=np.mean([h['lat'] for h in sel]); kx=np.cos(np.radians(lat0))
feats=[]; rng=np.random.default_rng(1)
for sb in SUBBASINS:
    poly=sbgeo[sb]; hs=[h for h in hrus if h['sb']==sb]
    pts=np.array([[h['lon']*kx,h['lat']] for h in hs])
    pts+=rng.normal(0,1e-4,pts.shape)                    # break ties between identical centroids
    minx,miny,maxx,maxy=poly.bounds; far=10*max(maxx-minx,maxy-miny)
    cx,cy=(minx+maxx)/2*kx,(miny+maxy)/2
    allpts=np.vstack([pts,[[cx-far,cy-far],[cx+far,cy-far],[cx-far,cy+far],[cx+far,cy+far]]])
    vor=Voronoi(allpts)
    for j,h in enumerate(hs):
        reg=vor.regions[vor.point_region[j]]
        cell=sg.Polygon([(vor.vertices[v][0]/kx,vor.vertices[v][1]) for v in reg])
        g=cell.intersection(poly)
        if g.is_empty: continue
        feats.append(dict(type='Feature',properties=dict(HRU_ID=h['id'],SubId=sb,AQUIFER_PROFILE=h['prof']),geometry=sg.mapping(g)))
json.dump(dict(type='FeatureCollection',features=feats),open('Liard_HRUs_standin.geojson','w'))

# ---- 2. DEM and 3. bedrock rasters (ESRI ASCII, lon/lat) by inverse-distance weighting of HRU values
dom=unary_union([sbgeo[s] for s in SUBBASINS]).buffer(0.3)
x0,y0,x1,y1=dom.bounds; cs=0.02
nc,nr=int(np.ceil((x1-x0)/cs)),int(np.ceil((y1-y0)/cs))
gx=x0+(np.arange(nc)+0.5)*cs; gy=y1-(np.arange(nr)+0.5)*cs
GX,GY=np.meshgrid(gx,gy)
src=[h for h in hrus if h['sb'] in SUBBASINS or dom.contains(sg.Point(h['lon'],h['lat']))]
sx=np.array([h['lon'] for h in src])*kx; sy=np.array([h['lat'] for h in src])
def idw(vals,p=2.0):
    d2=(GX[...,None]*kx-sx)**2+(GY[...,None]-sy)**2+1e-6
    w=1/d2**(p/2); return (w*vals).sum(-1)/w.sum(-1)
dem=idw(np.array([h['elev'] for h in src]))
depth=idw(np.array([60.0 if h['lu'] in VALLEY else (5.0 if h['lu']=='GLACIER' else 15.0) for h in src]),p=1.0)
def asc(fn,a):
    with open(fn,'w') as f:
        f.write(f"ncols {nc}\nnrows {nr}\nxllcorner {x0}\nyllcorner {y1-nr*cs}\ncellsize {cs}\nNODATA_value -9999\n")
        np.savetxt(f,a,fmt='%.2f')
asc('Liard_dem_standin.asc',dem); asc('Liard_bedrock_standin.asc',dem-depth)

# ---- 4. .rvh: AQUIFER_PROFILE column + AquiferHRUs group
aq=[h for h in sel if h['prof']!='[NONE]']
for h in aq:
    c=list(h['cols']); c[9]=h['prof']; lines[h['i']]=' , '.join(c)
grp='\n:HRUGroup AquiferHRUs\n  '+', '.join(str(h['id']) for h in aq)+'\n:EndHRUGroup\n'
open('Liard.rvh','w').write(head+':HRUs'+'\n'.join(lines)+':EndHRUs'+tail+grp)

# ---- 5. .rvi: switch on groundwater, route FAST_RESERVOIR percolation to GROUNDWATER in aquifer HRUs
rvi=open('Liard.rvi').read()
rvi=rvi.replace(':DefineHRUGroups WetlandHRUs',':GroundwaterModel MODFLOW6\n:DefineHRUGroups AquiferHRUs\n:DefineHRUGroups WetlandHRUs',1)
old=re.search(r':Percolation\s+PERC_CONSTANT\s+FAST_RESERVOIR\s+\tSLOW_RESERVOIR\r?\n\s*:-->Conditional HRU_GROUP IS_NOT NahanniHRUs\r?\n',rvi)
assert old, 'percolation line not found'
rvi=rvi.replace(old.group(0),old.group(0)+'    :-->Conditional HRU_GROUP IS_NOT AquiferHRUs\n'
 '  :Percolation           PERC_CONSTANT      FAST_RESERVOIR \tGROUNDWATER\n    :-->Conditional HRU_GROUP IS AquiferHRUs\n',1)
open('Liard.rvi','w').write(rvi)

# ---- 6. .rvp: hydrogeology
open('Liard.rvp','a').write('''
#-----------------------------------------------------------------
# Groundwater (Raven-MODFLOW 6): hypothetical aquifer for testing
#-----------------------------------------------------------------
:AquiferClasses
  :Attributes, K_HORIZ, K_VERT, SPEC_STORAGE, SPEC_YIELD, POROSITY
  :Units,      m/d,     m/d,    1/m,          none,       none
  GRAVEL,      50.0,    5.0,    1.0E-5,       0.25,       0.30
  CLAY_TILL,   0.001,   0.0001, 1.0E-4,       0.03,       0.40
  SAND,        10.0,    1.0,    1.0E-5,       0.20,       0.35
  UPLAND_TILL, 0.2,     0.02,   1.0E-4,       0.08,       0.30
:EndAquiferClasses

:AquiferProfiles
  VALLEY_ALLUVIUM, 3, GRAVEL, 15.0, AQUIFER, CLAY_TILL, 5.0, AQUITARD, SAND, TO_BEDROCK, CONFINED_AQUIFER
  UPLAND_TILL,     1, UPLAND_TILL, TO_BEDROCK, AQUIFER
:EndAquiferProfiles

:AquiferProfileParameters
  :Attributes,     INITIAL_HEAD_DEPTH, DEFAULT_BEDROCK_DEPTH, SEEPAGE_LEAKANCE, SOIL_ZONE_DEPTH
  :Units,          m,                  m,                     1/d,              m
  VALLEY_ALLUVIUM, 3.0,                60.0,                  1.0,              0.0
  UPLAND_TILL,     8.0,                15.0,                  1.0,              0.0
:EndAquiferProfileParameters
''')
# ---- 7. .rvg
open('Liard.rvg','w').write('''# Raven-MODFLOW 6 groundwater file: Liard test (subbasins 43, 48, 52)
:HRUGeometry       Liard_HRUs_standin.geojson  HRU_ID
:LandSurface       Liard_dem_standin.asc
:BedrockSurface    Liard_bedrock_standin.asc
:GridCellSize      2000.0
:MinCellCoverage   0.25
:HeadSaveFrequency 365
''')
print(f"HRU polygons: {len(feats)}  aquifer HRUs: {len(aq)} (valley {sum(h['prof']=='VALLEY_ALLUVIUM' for h in aq)}, upland {sum(h['prof']=='UPLAND_TILL' for h in aq)})  raster {nc}x{nr}")
