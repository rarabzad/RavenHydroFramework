"""Leland-format :OverlapWeights (w = A_overlap / A_cell) for the synthetic model, computed independently with shapely (dev-time only)."""
import os,sys; sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
import json,sys,pyproj
from shapely.geometry import shape
from shapely.ops import transform
from shapely.strtree import STRtree
md=sys.argv[1]; out=sys.argv[2]
tr=pyproj.Transformer.from_crs(4326,32610,always_xy=True); P=lambda g: transform(lambda x,y,z=None: tr.transform(x,y),g)
cells=[(f['properties']['CELL_ID'],P(shape(f['geometry']))) for f in json.load(open(md+'/ext_cells.geojson'))['features']]
tree=STRtree([c[1] for c in cells])
hru=json.load(open(tc.LIARD+'/Liard_HRUs_standin.geojson'))
n=0
with open(out,'w') as F:
    F.write(':OverlapWeights\n  :Attributes, HRU_ID, CELL, WEIGHT\n')
    for f in hru['features']:
        h=P(shape(f['geometry'])); hid=int(f['properties']['HRU_ID'])
        for i in tree.query(h):
            cid,c=cells[int(i)]; a=h.intersection(c).area
            if a>1e-6*c.area: F.write('  %d %d %.10f\n'%(hid,cid,a/c.area)); n+=1
    F.write(':EndOverlapWeights\n')
print(n,'weights written to',out)
