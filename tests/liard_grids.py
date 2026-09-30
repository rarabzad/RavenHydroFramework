# The Liard example on each grid type for two years (the manual's grid comparison, Table "The Liard model on each grid").
# Usage: python3 liard_grids.py            (runs every grid, then prints the table)
#        python3 liard_grids.py --table    (prints the table from earlier runs)
import os,shutil,subprocess,re,time,sys,json
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import testconfig as tc
tc.need_mf6(); B=tc.LIARD; ENV=tc.ENV; RV=tc.RAVEN; R=tc.work('liard_grids','results.json')
V=[('regular 2 km (default)',''),('regular 2 km, rotated 20 deg, 1 km along rivers',':GridRotation 20\n:GridRefinement RIVERS 1000\n'),
   ('quadtree 2 km to 500 m along rivers',':GridType QUADTREE\n:GridRefinement RIVERS 500\n'),
   ('HRU mesh 2 km, 1 km along rivers',':GridType HRU_MESH\n:GridRefinement RIVERS 1000\n'),
   ('regular 2 km, layers by elevation',':LayerConnection ELEVATION\n')]
def table(res):
    """the manual's table: active cells, run time, NSE at the three gauges, mean seepage to Raven and mean net aquifer
    discharge to rivers [million m3/d], plus the audit's recharge identity, MODFLOW residual and failed steps"""
    import pandas as pd
    print('%-50s %7s %6s  %-23s %8s %7s | %8s %9s %s'%('grid','cells','time','NSE 10BE004/005/010','seepage','rivers','recharge','residual','failed'))
    for r in res:
        o=r['dir']+'/out'
        if not r['ok']: print('%-50s FAILED'%r['name']); continue
        s=open(o+'/GWModelSummary.txt').read(); cells=int(re.search(r'active cells\s*:\s*(\d+)',s).group(1))
        d=pd.read_csv(o+'/Diagnostics.csv'); nse=[float(d[d['filename'].str.contains(g)]['DIAG_NASH_SUTCLIFFE'].iloc[0]) for g in ('10BE004','10BE005','10BE010')]
        b=pd.read_csv(o+'/GWBudget.csv'); seep=b['seepage to Raven [m3/d]'].mean()/1e6
        riv=(b['aquifer discharge to rivers [m3/d]']-b['river leakage to aquifer [m3/d]']).mean()/1e6
        A=tc.audit.audit(o)
        print('%-50s %7d %6.0f  %-23s %8.2f %7.2f | %8.0e %9.1e %s'%(r['name'],cells,r['sec'],' / '.join('%.3f'%x for x in nse),seep,riv,
              abs(A['recharge: (Raven - MF6 - dStore - toCHD)/Raven']),A['MF6 cumulative residual / gross flow'],A['unconverged steps']))
if '--table' in sys.argv:
    table(json.load(open(R))); sys.exit()
res=[]
for name,extra in V:
    d=tc.work('liard_grids',re.sub('[^a-z0-9]+','_',name.lower()).strip('_'))
    shutil.rmtree(d,ignore_errors=True); os.makedirs(d)
    for f in os.listdir(B):
        p=B+'/'+f
        if os.path.isfile(p) and not f.startswith('out'): shutil.copy(p,d)
        elif os.path.isdir(p) and f in ('data_obs','data_forcing'):
            try: os.symlink(p,d+'/'+f)
            except OSError: shutil.copytree(p,d+'/'+f)
    for f in os.listdir(B):
        if os.path.islink(B+'/'+f) and not os.path.exists(d+'/'+f): os.symlink(os.readlink(B+'/'+f),d+'/'+f)
    t=open(d+'/Liard.rvi').read(); t=re.sub(r':Duration\s+\S+(\s+\S+)?',':Duration 731',t); open(d+'/Liard.rvi','w').write(t)
    open(d+'/Liard.rvg','a').write('\n'+extra); os.makedirs(d+'/out',exist_ok=True)
    t0=time.time(); r=subprocess.run([RV,'Liard','-o','out/'],cwd=d,capture_output=True,text=True,env=ENV); dt=time.time()-t0
    ok='Successful' in (r.stdout + r.stderr)
    res.append(dict(name=name,dir=d,ok=ok,sec=dt)); json.dump(res,open(R,'w'),indent=1)
    print(name,'ok' if ok else 'FAILED',round(dt),'s',flush=True)
table(res)
print('results in',R)
