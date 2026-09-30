"""Synthetic 'existing' MODFLOW 6 model for the Liard coupled area (UTM 10N, rotated DIS grid, monthly periods).
Usage: python3 make_ext_model.py OUTDIR [meters_days|feet_seconds]"""
import os,sys; sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
import sys,os,json,math,datetime,numpy as np,flopy,pyproj
from shapely.geometry import shape, Polygon, Point
from shapely.ops import unary_union
out=sys.argv[1]; units=sys.argv[2] if len(sys.argv)>2 else 'meters_days'
L=3.280839895 if units=='feet_seconds' else 1.0          # model length units per metre
T=86400.0 if units=='feet_seconds' else 1.0              # model time units per day
B=tc.LIARD
hru=json.load(open(B+'/Liard_HRUs_standin.geojson'))
rvh=open(B+'/Liard.rvh').read()
# coupled HRUs = those with an AQUIFER_PROFILE other than [NONE]
import re
hdr=[l for l in rvh.split('\n') if ':Attributes' in l and 'AQUIFER_PROFILE' in l][0].replace(',',' ').split()[1:]
ia=hdr.index('AQUIFER_PROFILE'); coupled=set()
blk=rvh.split(':HRUs')[1].split(':EndHRUs')[0]
for l in blk.split('\n'):
    p=l.replace(',',' ').split()
    if p and p[0].isdigit() and len(p)>ia and p[ia+1] not in ('[NONE]',):   # +1: first token is the ID
        coupled.add(int(p[0]))
geoms=[shape(f['geometry']) for f in hru['features'] if int(f['properties']['HRU_ID']) in coupled]
dom=unary_union(geoms)
tr=pyproj.Transformer.from_crs(4326,32610,always_xy=True); inv=pyproj.Transformer.from_crs(32610,4326,always_xy=True)
from shapely.ops import transform
domu=transform(lambda x,y,z=None: tr.transform(x,y),dom).buffer(500)
ang=15.0; cs=1500.0
cx,cy=domu.centroid.x,domu.centroid.y; a=math.radians(ang)
# grid extent in the rotated frame
pts=np.array(domu.exterior.coords) if domu.geom_type=='Polygon' else np.vstack([np.array(g.exterior.coords) for g in domu.geoms])
lx=(pts[:,0]-cx)*math.cos(a)+(pts[:,1]-cy)*math.sin(a); ly=-(pts[:,0]-cx)*math.sin(a)+(pts[:,1]-cy)*math.cos(a)
x0,y0=lx.min()-cs,ly.min()-cs; ncol=int(math.ceil((lx.max()+cs-x0)/cs)); nrow=int(math.ceil((ly.max()+cs-y0)/cs))
xoff=cx+x0*math.cos(a)-y0*math.sin(a); yoff=cy+x0*math.sin(a)+y0*math.cos(a)
mg=flopy.discretization.StructuredGrid(delr=np.full(ncol,cs),delc=np.full(nrow,cs),xoff=xoff,yoff=yoff,angrot=ang)
XC,YC=mg.xcellcenters,mg.ycellcenters
# DEM (stand-in raster, lon/lat)
def read_asc(f):
    h={}; lines=open(f).read().split('\n')
    for i in range(6): k,v=lines[i].split(); h[k.lower()]=float(v)
    arr=np.loadtxt(f,skiprows=6); return h,arr
def sample(h,arr,lon,lat):
    c=((lon-h['xllcorner'])/h['cellsize']-0.5); r=h['nrows']-1-((lat-h['yllcorner'])/h['cellsize']-0.5)
    c=np.clip(np.round(c).astype(int),0,int(h['ncols'])-1); r=np.clip(np.round(r).astype(int),0,int(h['nrows'])-1); v=arr[r,c]
    return np.where(v==h.get('nodata_value',-9999),np.nan,v)
hd,dem=read_asc(B+'/Liard_dem_standin.asc'); hb,bed=read_asc(B+'/Liard_bedrock_standin.asc')
LON,LAT=inv.transform(XC,YC); top=sample(hd,dem,LON,LAT); br=sample(hb,bed,LON,LAT)
active=np.array([[domu.contains(Point(XC[r,c],YC[r,c])) and np.isfinite(top[r,c]) for c in range(ncol)] for r in range(nrow)])
top=np.where(np.isfinite(top),top,np.nanmean(top)); br=np.where(np.isfinite(br),br,top-150)
th=[20.0,10.0,40.0]; bot=np.zeros((4,nrow,ncol)); z=top.copy()
for k in range(3): z=z-th[k]; bot[k]=z
bot[3]=np.minimum(br,bot[2]-20.0)
idom=np.repeat(active[None,:,:].astype(int),4,axis=0)
rng=np.random.default_rng(11)
def field(mean,sd): f=rng.normal(0,1,(nrow,ncol)); 
f=None
def smooth(m,sd):
    g=rng.normal(0,1,(nrow,ncol))
    for _ in range(4): g=(g+np.roll(g,1,0)+np.roll(g,-1,0)+np.roll(g,1,1)+np.roll(g,-1,1))/5
    g=g/g.std(); return 10**(np.log10(m)+sd*g)
K=np.stack([smooth(40,0.3),smooth(0.05,0.3),smooth(8,0.3),smooth(2,0.2)])
K33=K/10.0
# periods: steady 1 day, then 24 calendar months from 1985-10-01
start=datetime.date(1985,10,1); per=[(1.0,1)]
d=start
for m in range(24):
    n=(datetime.date(d.year+(d.month//12),d.month%12+1,1)-d).days; per.append((float(n),n)); d=d+datetime.timedelta(days=n)
sim=flopy.mf6.MFSimulation(sim_name='ext',sim_ws=out,exe_name='mf6')
flopy.mf6.ModflowTdis(sim,time_units='seconds' if T!=1 else 'days',nper=len(per),perioddata=[(p*T,n,1.0) for p,n in per],start_date_time='1985-10-01')
flopy.mf6.ModflowIms(sim,complexity='complex',outer_dvclose=1e-5*L,outer_maximum=500,inner_dvclose=1e-6*L,rcloserecord=[1.0*L**3/T,'strict'],linear_acceleration='bicgstab')
gwf=flopy.mf6.ModflowGwf(sim,modelname='gwf_liard',newtonoptions='newton',save_flows=True)   # plain NEWTON: its UNDER_RELAXATION keyword stalls MODFLOW 6.6 on drying cells
flopy.mf6.ModflowGwfdis(gwf,length_units='feet' if L!=1 else 'meters',nlay=4,nrow=nrow,ncol=ncol,delr=cs*L,delc=cs*L,top=top*L,botm=bot*L,idomain=idom,xorigin=xoff*L,yorigin=yoff*L,angrot=ang)
flopy.mf6.ModflowGwfnpf(gwf,icelltype=[1,1,1,0],k=K*L/T,k33=K33*L/T,save_flows=True)
flopy.mf6.ModflowGwfsto(gwf,iconvert=[1,1,1,0],ss=1e-5/L,sy=[0.2,0.05,0.15,0.1],steady_state={0:True},transient={1:True})
flopy.mf6.ModflowGwfic(gwf,strt=np.repeat(((top-3.0)*L)[None,:,:],4,axis=0))
# rivers along the river lines (subbasin ID as auxiliary)
riv=json.load(open(B+'/Liard_rivers.geojson')); gi=flopy.utils.GridIntersect(mg)
from shapely.geometry import LineString
rivcells={}
for f in riv['features']:
    g=shape(f['geometry']); sb=f['properties'].get('SubId',-1)
    gu=transform(lambda x,y,z=None: tr.transform(x,y),g)
    if not gu.intersects(domu): continue
    res=gi.intersect(gu)
    for rec in res:
        r,c=rec['cellids']; 
        if not active[r,c]: continue
        Lr=rec['lengths']; key=(r,c)
        if key in rivcells: rivcells[key][0]+=Lr
        else: rivcells[key]=[Lr,int(sb)]
sub=json.load(open(B+'/Liard_subbasins.geojson')); subg=[(shape(f['geometry']),int(str(f['properties']['Sub_B']).split('_')[-1])) for f in sub['features']]
def subbasin_at(lon,lat):
    pt=Point(lon,lat)
    for g,i in subg:
        if g.contains(pt): return i
    return -1
rivdata=[]
for (r,c),(Lr,sb) in rivcells.items():
    lo,la=inv.transform(XC[r,c],YC[r,c]); sb=subbasin_at(lo,la)
    if sb<0: continue
    stage=top[r,c]-1.5; rb=stage-2.0; cond=0.5*20.0*Lr/1.0   # Kb=0.5 m/d, width 20 m, thickness 1 m
    rivdata.append([(0,r,c),stage*L,cond*L**2/T,rb*L,float(sb)])
flopy.mf6.ModflowGwfriv(gwf,pname='RIV_STREAMS',auxiliary=['SUBBASIN'],stress_period_data={0:rivdata},save_flows=True,maxbound=len(rivdata))
# own recharge and evaporation (to be taken over by Raven)
flopy.mf6.ModflowGwfrcha(gwf,pname='RCH_MODEL',recharge=0.0004*L/T,save_flows=True)
flopy.mf6.ModflowGwfevta(gwf,pname='EVT_MODEL',surface=top*L,rate=0.001*L/T,depth=2.0*L,save_flows=True)
# wells: two pumping wells with seasonal rates; general-head boundary at the lowest active cells (outlet)
ac=np.argwhere(active); tops=[top[r,c] for r,c in ac]; low=ac[np.argsort(tops)[:6]]
flopy.mf6.ModflowGwfghb(gwf,pname='GHB_OUTLET',stress_period_data={0:[[(3,int(r),int(c)),(top[r,c]-5)*L,2000.0*L**2/T] for r,c in low]},save_flows=True)
mid=ac[len(ac)//2]; w2=ac[len(ac)//3]
wel={}
for p in range(len(per)):
    q=-8000.0 if (p==0 or (p-1)%12 in (5,6,7,8)) else -3000.0     # summer pumping higher (months Mar-Jun from Oct)
    wel[p]=[[(2,int(mid[0]),int(mid[1])),q*L**3/T],[(2,int(w2[0]),int(w2[1])),0.5*q*L**3/T]]
flopy.mf6.ModflowGwfwel(gwf,pname='WEL_PUMPS',stress_period_data=wel,save_flows=True)
flopy.mf6.ModflowGwfoc(gwf,head_filerecord='gwf_liard.hds',budget_filerecord='gwf_liard.cbc',saverecord=[('HEAD','LAST'),('BUDGET','LAST')])
sim.write_simulation(silent=True)
# cell polygons in lon/lat, CELL_ID = row*ncol+col+1 (the cell number of a layer)
feats=[]
for r in range(nrow):
    for c in range(ncol):
        v=mg.get_cell_vertices(r,c); lonlat=[list(inv.transform(x,y)) for x,y in v]; lonlat.append(lonlat[0])
        feats.append({'type':'Feature','properties':{'CELL_ID':r*ncol+c+1},'geometry':{'type':'Polygon','coordinates':[lonlat]}})
json.dump({'type':'FeatureCollection','features':feats},open(out+'/ext_cells.geojson','w'))
json.dump({'nrow':nrow,'ncol':ncol,'active':int(active.sum()),'riv':len(rivdata),'units':units,'periods':len(per)},open(out+'/ext_info.json','w'))
print('model written:',nrow,'x',ncol,'x 4 layers; active columns',int(active.sum()),'; river cells',len(rivdata),'; periods',len(per),'; units',units)
