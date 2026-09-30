"""Lake refinement: a reservoir on a coupled HRU of subbasin 52; LAKES refinement on each refinable grid type,
POLYGONS on the HRU mesh; audited; cell sizes inside/around the lake against elsewhere."""
import os,shutil,subprocess,re,sys,json,math
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import testconfig as tc; audit=tc.audit
from shapely.geometry import shape
tc.need_mf6(); B=tc.grid_base(); ENV=tc.ENV; RV=tc.RAVEN; LAKEPOLY=tc.work('lake','lake_poly.geojson')
rvh=open(B+'/Liard.rvh').read(); hdr=[l for l in rvh.split('\n') if ':Attributes' in l and 'AQUIFER_PROFILE' in l][0].replace(',',' ').split()[1:]
ib,ia,ie,iar=hdr.index('BASIN_ID')+1,hdr.index('AQUIFER_PROFILE')+1,hdr.index('ELEVATION')+1,hdr.index('AREA')+1
rows=[l.replace(',',' ').split() for l in rvh.split(':HRUs')[1].split(':EndHRUs')[0].split('\n')]
cand=[r for r in rows if r and r[0].isdigit() and r[ib]=='52' and r[ia]!='[NONE]']
lake=max(cand,key=lambda r:float(r[iar])); lid=int(lake[0]); elev=float(lake[ie])
H=json.load(open(B+'/Liard_HRUs_standin.geojson'))
poly=[f for f in H['features'] if int(f['properties']['HRU_ID'])==lid][0]
json.dump({'type':'FeatureCollection','features':[poly]},open(LAKEPOLY,'w'))
LP=shape(poly['geometry']); print('lake HRU',lid,'area %.1f km2 (rvh), elevation %.0f m'%(float(lake[iar]),elev))
res=f"\n:Reservoir TestLake\n  :SubBasinID 52\n  :HRUID {lid}\n  :Type RESROUTE_STANDARD\n  :WeirCoefficient 0.6\n  :CrestWidth 20.0\n  :MaxDepth 8.0\n  :LakeArea 2.0E7\n  :AbsoluteCrestHeight {elev:.2f}\n  :SeepageParameters 0.05 {elev-5:.2f}\n:EndReservoir\n"
V=[('regular, rotated, LAKES 750',':GridRotation 20\n:GridRefinement LAKES 750\n'),('quadtree, LAKES 375',':GridType QUADTREE\n:GridRefinement LAKES 375\n'),
   ('HRU mesh, LAKES 500',':GridType HRU_MESH\n:GridRefinement LAKES 500\n'),('HRU mesh, POLYGONS 500',':GridType HRU_MESH\n:GridRefinement POLYGONS lake_poly.geojson 500\n'),
   ('HRU mesh, no refinement',':GridType HRU_MESH\n')]
for name,extra in V:
    d=tc.work('lake','lk_'+re.sub('[^a-z0-9]+','_',name.lower()).strip('_')); shutil.rmtree(d,ignore_errors=True); shutil.copytree(B,d,ignore=shutil.ignore_patterns('out')); os.makedirs(d+'/out')
    shutil.copy(LAKEPOLY,d); open(d+'/Liard.rvh','a').write(res)
    g=open(d+'/Liard.rvg').read(); open(d+'/Liard.rvg','w').write(g+'\n:ReservoirExchange\n'+extra)
    t=open(d+'/Liard.rvi').read(); t=re.sub(r':Duration\s+\S+',':Duration 60',t); open(d+'/Liard.rvi','w').write(t)
    r=subprocess.run([RV,'Liard','-o','out/'],cwd=d,capture_output=True,text=True,env=ENV)
    ex=(re.findall(r'Exiting Gracefully: (.*)',r.stdout + r.stderr) or ['(no exit line)'])[-1]
    if 'Successful' not in ex: print('%-28s FAILED: %s'%(name,ex)); continue
    cells=json.load(open(d+'/out/mf6/grid_cells.geojson'))['features']
    near=[];far=[]
    for f in cells:
        if not f['properties']['active']: continue
        c=shape(f['geometry']); s=math.sqrt(f['properties']['area_m2'])
        (near if c.intersects(LP.buffer(0.01)) else far).append(s)
    A=audit.audit(d+'/out')
    print('%-28s cells near the lake: %3d, median %4.0f m | elsewhere: %4d, median %4.0f m | recharge %.0e residual %.1e river %.0e failed %s'%(name,len(near),sorted(near)[len(near)//2],len(far),sorted(far)[len(far)//2],
          abs(A['recharge: (Raven - MF6 - dStore - toCHD)/Raven']),A['MF6 cumulative residual / gross flow'],abs(A['river: applied - net MF6 - carry change (rel)']),A['unconverged steps']))
    b=[k for k in A if 'reservoir' in k.lower()]
    print('      reservoir checks:',{k:A[k] for k in b})
