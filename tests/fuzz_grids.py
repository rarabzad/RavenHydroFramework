import os, re, sys, shutil, subprocess, json, random, math
import pandas as pd, numpy as np, flopy
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import testconfig as tc
# usage: fuzz_grids.py SEED N [lakes]   ('lakes': every set-up has a reservoir, a refinable grid and lake refinement)
tc.need_mf6(); SEED=sys.argv[1] if len(sys.argv)>1 else '7'; LAKES=(len(sys.argv)>3) and (sys.argv[3]=='lakes')
B=tc.LIARD; RV=tc.RAVEN; ENV=tc.ENV; W=tc.work('fuzz','runs_s'+SEED); os.makedirs(W,exist_ok=True)
RVP_TEMPLATE=tc.liard_rvp_template(); QT_CELLS=tc.quadtree_cells()
rvh=open(f'{B}/Liard.rvh').read(); body=rvh.split(':HRUs')[1].split(':EndHRUs')[0]
hrus=[[x.strip() for x in l.split(',')] for l in body.split('\n')]
hrus=[c for c in hrus if len(c)>12 and c[0].isdigit() and c[9] not in ('[NONE]','NONE','') and c[5]=='52']   # the subbasin-52 model (testconfig.sb52)
def _dem():
    """the example's stand-in DEM (ESRI ASCII, lon/lat)"""
    L=open(f'{B}/Liard_dem_standin.asc').read().split('\n'); h={l.split()[0].lower():float(l.split()[1]) for l in L[:6]}
    z=np.array([[float(v) for v in l.split()] for l in L[6:] if l.strip()]); return h,z
DEMH,DEMZ=_dem()
from shapely.geometry import shape as _shape
HPOLY={str(f['properties']['HRU_ID']):_shape(f['geometry']) for f in json.load(open(f'{B}/Liard_HRUs_standin.geojson'))['features']}
def land_at_lake(hid):
    """highest land surface (stand-in DEM) inside the lake HRU's polygon: the reservoir crest is set there, so the lake
    stands on the aquifer's datum above every cell under it (the .rvh ELEVATION of the stand-in HRUs can lie below the
    DEM-based aquifer, which Raven rightly refuses)"""
    from shapely.geometry import Point
    P=HPOLY[str(hid)]; x0,y0,x1,y1=P.bounds; cs=DEMH['cellsize']; nr=int(DEMH['nrows']); zmax=-1e30
    for r in range(nr):
        y=DEMH['yllcorner']+(nr-r-0.5)*cs
        if y<y0-cs or y>y1+cs: continue
        for c in range(int(DEMH['ncols'])):
            x=DEMH['xllcorner']+(c+0.5)*cs
            if x0-cs<=x<=x1+cs and DEMZ[r,c]!=DEMH['nodata_value'] and P.buffer(cs).contains(Point(x,y)): zmax=max(zmax,float(DEMZ[r,c]))
    return zmax
def setup(name,cfg,days,start='1985-10-01'):
    d=f'{W}/{name}'
    if os.path.exists(d): shutil.rmtree(d)
    os.makedirs(d)
    for f in os.listdir(B):
        if f.endswith(('.rvi','.rvh','.rvt','.rvc','.rvg','.rvp','.geojson','.asc')) and f!='Liard.rvp': shutil.copy(f'{B}/{f}',d)
    for l in ('data_forcing','data_obs'):
        try: os.symlink(f'{B}/{l}',f'{d}/{l}')
        except OSError: shutil.copytree(f'{B}/{l}',f'{d}/{l}')
    tc.sb52(d)   # groundwater under subbasin 52 only; wells, boundaries and the reservoir are placed in its HRUs
    rvp=RVP_TEMPLATE.replace('par_K_gravel',str(cfg['Kg'])).replace('par_Sy_gravel','0.25').replace('par_K_till',str(cfg['Kt']))
    # profile parameters: ET, delay, smoothing
    rvp=re.sub(r':AquiferProfileParameters.*?:EndAquiferProfileParameters',
      ':AquiferProfileParameters\n  :Attributes, INITIAL_HEAD_DEPTH, SEEPAGE_LEAKANCE, EXTINCTION_DEPTH, RECHARGE_DELAY, SEEPAGE_SMOOTHING_DEPTH\n  :Units, m, 1/d, m, d, m\n'
      f"  VALLEY_ALLUVIUM, {cfg['h0']}, {cfg['leak']}, {cfg['etV']}, {cfg['delay']}, {cfg['smooth']}\n  UPLAND_TILL, {cfg['h0']}, {cfg['leak']}, {cfg['etU']}, 0.0, {cfg['smooth']}\n:EndAquiferProfileParameters",rvp,flags=re.S)
    open(f'{d}/Liard.rvp','w').write(rvp)
    x=open(f'{d}/Liard.rvi').read(); x=re.sub(r':Duration\s+\S+',f':Duration {days}',x); x=re.sub(r':StartDate\s+\S+\s+\S+',f':StartDate {start} 00:00:00',x)
    if cfg['dt']!=1: x=x.replace(':TimeStep',f":TimeStep {cfg['dt']}\n#:TimeStep",1)
    open(f'{d}/Liard.rvi','w').write(x)
    g=open(f'{d}/Liard.rvg').read(); g=re.sub(r':GridCellSize\s+\S+',f":GridCellSize {cfg['cs']}",g); g=re.sub(r':MinCellCoverage\s+\S+',f":MinCellCoverage {cfg['cov']}",g)
    if ':MinCellCoverage' not in g: g+=f":MinCellCoverage {cfg['cov']}\n"
    add=''
    if cfg['ss']: add+=':GWInitialization STEADY_STATE 0.3\n'
    if cfg['swx']: add+=f":SurfaceWaterExchange SOIL[0] AquiferHRUs LEAKANCE {cfg['swx']}\n"
    if cfg['wells']:
        add+=':Wells\n  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_TOP, SCREEN_BOT\n'
        for i,(h,q) in enumerate(cfg['wells']): add+=f"  {i+1}, W{i+1}, {h[3]}, {h[4]}, {float(h[2])-4}, {float(h[2])-35}\n"
        add+=':EndWells\n'
    if cfg['obs']:
        add+=':ObservationWells\n  :Attributes, NAME, LATITUDE, LONGITUDE, SCREEN_ELEV\n'
        for i,h in enumerate(cfg['obs']): add+=f"  {101+i}, OW{i}, {h[3]}, {h[4]}, {float(h[2])-10}\n"
        add+=':EndObservationWells\n'
    if cfg['ghb']: add+=f":GeneralHeadBoundary G1 Liard_rivers.geojson HEAD {cfg['ghb']} CONDUCTANCE 500 LAYERS TOP\n"
    if cfg['chd']:
        h=cfg['chd']; la,lo=float(h[3]),float(h[4]); dd=0.03
        json.dump({"type":"FeatureCollection","features":[{"type":"Feature","properties":{},"geometry":{"type":"Polygon","coordinates":[[[lo-dd,la-dd],[lo+dd,la-dd],[lo+dd,la+dd],[lo-dd,la+dd],[lo-dd,la-dd]]]}}]},open(f'{d}/chd.geojson','w'))
        add+=f":SpecifiedHeadBoundary L1 chd.geojson HEAD {float(h[2])-1.0} LAYERS TOP\n"
    if cfg['res']: add+=':ReservoirExchange\n'
    G=cfg['grid']; cs=cfg['cs']
    if G=='ROT': add+=f":GridRotation {cfg['rot'] or 30}\n"
    if G=='QT': add+=':GridType QUADTREE\n'+(f":GridRotation {cfg['rot']}\n" if cfg['rot'] else '')
    if G=='MESH': add+=':GridType HRU_MESH\n'
    if G=='ELEV': add+=':LayerConnection ELEVATION\n'
    if G=='FILE': shutil.copy(QT_CELLS,d); add+=':GridType FILE\n:GridFile qt_cells.geojson cell\n'
    if G in ('ROT','QT','MESH'):
        if cfg['refr'] or G=='QT': add+=f":GridRefinement RIVERS {cs*(cfg['refr'] or 0.25)}\n"
        if cfg['refw'] and (cfg['wells'] or cfg['obs']): add+=f":GridRefinement WELLS {cs*cfg['refw']}\n"
    if cfg['fc']!='AUTO': add+=f":FlowCorrection {cfg['fc']}\n"
    if cfg['res'] and G in ('ROT','QT','MESH') and cfg['lake']:
        if cfg['lake']=='poly' and G=='MESH':
            h=cfg['res']; H=json.load(open(B+'/Liard_HRUs_standin.geojson'))
            json.dump({'type':'FeatureCollection','features':[f for f in H['features'] if str(f['properties']['HRU_ID'])==str(h[0])]},open(d+'/lake_poly.geojson','w'))
            add+=f":GridRefinement POLYGONS lake_poly.geojson {cs*0.25}\n"
        else: add+=f":GridRefinement LAKES {cs*0.25}\n"
    open(f'{d}/Liard.rvg','w').write(g+add)
    if cfg['wells']:
        n=int(days/cfg['dt'])+2; t=f"\n:WellRate 1 m3/d\n  {start} 00:00:00 {cfg['dt']} {n}\n"+''.join(f"  {cfg['wells'][0][1]}\n" for _ in range(n))+":EndWellRate\n"
        for i,(h,q) in enumerate(cfg['wells'][1:],start=2): t+=f":WellRate {i} m3/d\n  {start} 00:00:00 {cfg['dt']} {n}\n"+''.join(f"  {q}\n" for _ in range(n))+":EndWellRate\n"
        open(f'{d}/Liard.rvt','a').write(t)
    if cfg['res']:
        h=cfg['res']; zc=land_at_lake(h[0])
        open(f'{d}/Liard.rvh','a').write(f"\n:Reservoir FuzzLake\n  :SubBasinID 52\n  :HRUID {h[0]}\n  :Type RESROUTE_STANDARD\n  :WeirCoefficient 0.6\n  :CrestWidth 20.0\n  :MaxDepth 8.0\n  :LakeArea 2.0E7\n  :AbsoluteCrestHeight {zc:.2f}\n  :SeepageParameters {cfg['resk']} {zc-5:.2f}\n:EndReservoir\n")
    return d
def run(d):
    os.makedirs(d+'/out',exist_ok=True)
    try: r=subprocess.run([RV,'Liard','-o','out/'],cwd=d,capture_output=True,text=True,timeout=600,env=ENV); so=r.stdout+r.stderr; rc=r.returncode
    except subprocess.TimeoutExpired: return False,'timeout'
    m=re.findall(r'Exiting Gracefully: (.*)',so); msg=m[-1] if m else ''
    return ('Successful Simulation' in so), (msg if 'Successful' not in msg else '') + ('' if rc==0 else f' rc={rc}')
def check(d,cfg):
    out=d+'/out/'; b=pd.read_csv(out+'GWBudget.csv').rename(columns=lambda c:c.strip()).dropna(subset=['Raven recharge [m3/d]'])
    dt=cfg['dt']; P=[]
    tot=b['Raven recharge [m3/d]'].sum()*dt
    e=(tot-b['MF6 recharge [m3/d]'].sum()*dt-b['recharge delay storage [m3]'].iloc[-1]-b['recharge onto CHD cells [m3/d]'].sum()*dt)/max(tot,1); P.append(('recharge',abs(e)<1e-8,e))
    conv=b['converged'].astype(int)>=1
    e=((b['returned to Raven [m3]']-b['seepage to Raven [m3/d]']*dt)[conv]).abs().sum()/max((b['seepage to Raven [m3/d]']*dt).sum(),1); P.append(('seepage',e<1e-9,e))
    net=(b['aquifer discharge to rivers [m3/d]']-b['river leakage to aquifer [m3/d]'])*dt
    e=(b['river exchange applied to Raven [m3]'].sum()-net.sum()-b['river loss carried to next step [m3]'].iloc[-1])/max(net.abs().sum(),1); P.append(('river',abs(e)<1e-8,e))
    if cfg['swx']:
        sw=(b['aquifer discharge to surface stores [m3/d]']-b['surface-store leakage to aquifer [m3/d]'])*dt
        e=(b['surface-store exchange applied to Raven [m3]']-sw)[conv].abs().sum()/max(sw.abs().sum(),1); P.append(('stores',e<1e-8,e))
    if cfg['res']:
        s=b['reservoir seepage sent [m3/d]'].values*dt; pnd=b['reservoir seepage awaiting delivery [m3]'].values
        e=np.abs(s[1:]-pnd[:-1]).max()/max(np.abs(s).max(),1); P.append(('res lag',e<1e-9,e))
    err=b['MF6 balance error [%]'][conv].abs().max(); P.append(('MF6 step',err<0.1,err))
    gc=[c for c in ['MF6 recharge [m3/d]','seepage to Raven [m3/d]','storage release [m3/d]','river leakage to aquifer [m3/d]','aquifer discharge to rivers [m3/d]','water-table ET [m3/d]','wells actual [m3/d]','GHB [m3/d]','CHD [m3/d]','surface-store leakage to aquifer [m3/d]','aquifer discharge to surface stores [m3/d]','reservoir seepage delivered [m3/d]'] if c in b]
    gross=sum(b[c].abs().sum() for c in gc)
    ce=b['MF6 balance error [m3/d]'].sum()/max(gross,1); P.append(('MF6 cum',abs(ce)<1e-4,ce))
    w=pd.read_csv(out+'WatershedStorage.csv').rename(columns=lambda c:c.strip()); mb=w['MB Error [mm]']
    drift=mb.abs().max(); P.append(('Raven MB',drift<1.0,drift))
    inc,_=flopy.utils.Mf6ListBudget(out+'mf6/gwf.lst').get_dataframes(start_datetime=None); last=inc.iloc[-1]; bl=b.iloc[-1]
    g=lambda c: bl[c] if c in bl else 0.0
    ty={'RCH':g('MF6 recharge [m3/d]'),'DRN':-g('seepage to Raven [m3/d]'),'EVT':-g('water-table ET [m3/d]'),'GHB':g('GHB [m3/d]'),'CHD':g('CHD [m3/d]'),
        'RIV':g('river leakage to aquifer [m3/d]')-g('aquifer discharge to rivers [m3/d]')+g('surface-store leakage to aquifer [m3/d]')-g('aquifer discharge to surface stores [m3/d]'),
        'WEL':g('wells actual [m3/d]')+g('reservoir seepage delivered [m3/d]'),'STO':g('storage release [m3/d]')}
    worst=0
    for t,v in ty.items():
        ci=[c for c in last.index if c.startswith(t) and c.endswith('_IN')]; co=[c for c in last.index if c.startswith(t) and c.endswith('_OUT')]
        if not ci: continue
        mf=sum(last[c] for c in ci)-sum(last[c] for c in co); worst=max(worst,abs(mf-v)/max(abs(mf),abs(v),1000.0))
    P.append(('list-file',worst<1e-5,worst))
    _au=tc.audit
    A=_au.audit(d+'/out/')
    if 'reservoir: sent(t) - awaiting(t-1) (rel)' in A:
        v=A['reservoir: sent(t) - awaiting(t-1) (rel)']; P.append(('reservoir lag',v<1e-8,v))
    return P
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 7); N=int(sys.argv[2]) if len(sys.argv)>2 else 24
res=[]
for i in range(N):
    pick=lambda: random.choice(hrus)
    cfg=dict(cs=random.choice([1000,2000,3000]),cov=random.choice([0.1,0.25,0.5]),dt=random.choice([1,1,1,0.5]),Kg=round(10**random.uniform(0.5,2.5),2),Kt=round(10**random.uniform(-2,0.5),3),
             h0=random.choice([0.5,3,10]),leak=random.choice([0.1,1,5]),etV=random.choice([0,0,2]),etU=random.choice([0,1.5]),delay=random.choice([0,0,30]),smooth=random.choice([0,0.5]),
             ss=random.random()<0.3,swx=random.choice([0,0,0.02]),wells=[(pick(),-random.choice([5000,20000,60000])) for _ in range(random.choice([0,0,1,2]))],
             obs=[pick() for _ in range(random.choice([0,2]))],ghb=random.choice([0,0,1250]),chd=random.choice([None,None,pick()]),res=pick() if LAKES else random.choice([None,pick(),pick()]),resk=random.choice([0.02,0.2]),
             grid=random.choice(['ROT','QT','MESH','MESH'] if LAKES else ['REG','ROT','QT','QT','MESH','MESH','ELEV','FILE']),rot=random.choice([0,17,30,-45]),refr=random.choice([0,0.25,0.5]),refw=random.choice([0,0.2]),fc=random.choice(['AUTO','AUTO','NONE']),lake=random.choice(['lakes','lakes','poly'] if LAKES else [None,'lakes','lakes','poly']))
    if cfg['grid'] in ('ROT','QT','MESH') and (cfg['refr'] or cfg['lake']): cfg['cs']=max(cfg['cs'],2000)
    days=random.choice([45,90])
    name=f'f{i:02d}'; d=setup(name,cfg,days); ok,msg=run(d)
    rec=dict(name=name,cfg={k:(v if k not in ('wells','obs','chd','res') else bool(v)) for k,v in cfg.items()},ok=ok,msg=msg)
    if ok:
        try: P=check(d,cfg); rec['checks']=[(n,bool(p),float(v)) for n,p,v in P]; rec['pass']=bool(all(bool(p) for _,p,_ in P))
        except Exception as ex: rec['pass']=False; rec['exc']=repr(ex)[:200]
        if rec.get('pass') and random.random()<0.35 and cfg['dt']==1:   # hotstart split test
            half=days//2; d1=setup(name+'a',cfg,half); ok1,_=run(d1)
            import datetime; st=(datetime.date(1985,10,1)+datetime.timedelta(days=half)).isoformat()
            d2=setup(name+'b',dict(cfg,ss=False),days-half,start=st); shutil.copy(d1+'/out/solution.rvc',d2+'/Liard.rvc'); ok2,m2=run(d2)
            if ok1 and ok2:
                c=pd.read_csv(d+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]']); h=pd.read_csv(d2+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]']); m=c.merge(h,on='date',suffixes=('_c','_h'))
                cols=['seepage to Raven [m3/d]','MF6 recharge [m3/d]','storage release [m3/d]','river exchange applied to Raven [m3]']
                # deviation relative to each flow's largest value in the run (as in existing/fuzz_ext.py): a small net flow such
                # as the storage release can differ by a fraction of a cubic metre, which is not a restart error
                dev=max(((m[x+'_c']-m[x+'_h']).abs()/(m[x+'_c'].abs().max()+1e3)).max() for x in cols)
                rec['hotstart_dev']=float(dev); rec['pass']=bool(rec['pass'] and dev<1e-4)
            else: rec['pass']=False; rec['hot_msg']=m2
    else: rec['pass']=False
    res.append(rec); json.dump(res,open(tc.work('fuzz','results_'+SEED+('_lakes' if LAKES else '')+'.json'),'w'),indent=1,default=lambda o: o.item() if hasattr(o,'item') else str(o))
print('results in',tc.work('fuzz','results_'+SEED+('_lakes' if LAKES else '')+'.json'))
