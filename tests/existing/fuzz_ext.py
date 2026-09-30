"""Randomized tests of the existing-model mode. Usage: fuzz_ext.py SEED N"""
import os,sys; sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
import os,shutil,subprocess,re,sys,random,json,datetime,pandas as pd
audit=tc.audit; import flopy
tc.need_mf6(); X=tc.ext_models(); ENV=tc.ENV; RV=tc.RAVEN
seed=int(sys.argv[1]); N=int(sys.argv[2]); random.seed(seed); W=X+'/fz%d'%seed; os.makedirs(W,exist_ok=True); res=[]
def prep(d,cfg,start,dur,rvc=None):
    md=X+('/md' if cfg['units']=='md' else '/fs')
    subprocess.run([sys.executable,tc.SETUP_EXT,d,md],capture_output=True)
    t=open(d+'/Liard.rvi').read(); t=re.sub(r':StartDate\s+\S+\s+\S+',':StartDate %s 00:00:00'%start,t); t=re.sub(r':Duration\s+\S+',':Duration %d'%dur,t); open(d+'/Liard.rvi','w').write(t)
    g=open(d+'/Liard.rvg').read()
    if cfg['link']=='crs': g=re.sub(r':MF6GridFile.*',':MF6CRS            EPSG:32610',g)
    if cfg['link']=='weights': g=re.sub(r':MF6GridFile.*\n','',g); g=re.sub(r':HRUGeometry.*\n','',g); g+=open(X+'/md/weights.rvg').read()
    g=g.replace(':MinCellCoverage   0.25',':MinCellCoverage   %s'%cfg['cov'])
    if cfg['stream']=='none': g=re.sub(r':StreamPackage.*\n','',g)
    elif cfg['stream']=='dominant': g=g.replace('RIV_STREAMS AUX SUBBASIN','RIV_STREAMS DOMINANT_HRU')
    if not cfg['seep']: g=g.replace(':SeepageToRaven    ON',':SeepageToRaven    OFF')
    if cfg['evt']=='keep': g=g.replace(':RavenTakesOver    RCH_MODEL EVT_MODEL',':RavenTakesOver    RCH_MODEL\n:KeepPackages      EVT_MODEL')
    open(d+'/Liard.rvg','w').write(g)
    if rvc: shutil.copy(rvc,d+'/Liard.rvc')
def run(d):
    try: r=subprocess.run([RV,'Liard','-o','out/'],cwd=d,capture_output=True,text=True,env=ENV,timeout=600)
    except subprocess.TimeoutExpired: return 'timeout'
    return (re.findall(r'Exiting Gracefully: (.*)',r.stdout + r.stderr) or ['(crash rc %d)'%r.returncode])[-1]
def listcheck(d):
    """MODFLOW's own saved flows of RCH_RAVEN and DRN_RAVEN (budget file, at each saved time) against the coupler's budget [max relative]"""
    cbf=[f for f in os.listdir(d+'/out/mf6') if f.endswith('.cbc')][0]
    cb=flopy.utils.CellBudgetFile(d+'/out/mf6/'+cbf,precision='double')
    b=pd.read_csv(d+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]'])
    fs='feet' in open(d+'/out/mf6/gwf_liard.dis').read().lower(); fl=3.280839895 if fs else 1.0; ft=86400.0 if fs else 1.0
    worst=0.0; n=0
    for idx,h in enumerate(cb.recordarray):
        pk=h['paknam2'].decode().strip().upper() if 'paknam2' in h.dtype.names else ''
        if pk not in ('RCH_RAVEN','DRN_RAVEN'): continue
        q=cb.get_record(idx)['q'].sum()*ft/fl**3   # [m3/d], positive into the aquifer
        t=h['totim']/ft
        row=b[(b['time [d]']-t).abs()<1e-6]
        if row.empty: continue
        ours=row['MF6 recharge [m3/d]'].iloc[0] if pk=='RCH_RAVEN' else -row['seepage to Raven [m3/d]'].iloc[0]
        worst=max(worst,abs(q-ours)/max(abs(ours),1.0)); n+=1
    return worst if n>0 else float('nan')
for i in range(N):
    cfg=dict(units=random.choice(['md','md','fs']),link=random.choice(['crs','crs','file','weights']),stream=random.choice(['aux','aux','dominant','none']),
             seep=random.random()<0.8,evt=random.choice(['takeover','takeover','keep']),cov=random.choice([0.1,0.25,0.5]),days=random.choice([45,90]))
    name='f%02d'%i; d=W+'/'+name; prep(d,cfg,'1985-10-01',cfg['days']); msg=run(d); rec=dict(name=name,cfg=cfg,msg=msg[:120])
    if 'Successful' in msg:
        A=audit.audit(d+'/out')
        chk=[('recharge',abs(A['recharge: (Raven - MF6 - dStore - toCHD)/Raven'])<1e-9),('residual',abs(A['MF6 cumulative residual / gross flow'])<1e-7),
             ('river',abs(A['river: applied - net MF6 - carry change (rel)'])<1e-9),('failed steps',A['unconverged steps']==0)]
        if cfg['seep']: chk.append(('seepage',abs(A.get('seepage returned vs MF6 seepage (rel)',0))<1e-9))
        lc=listcheck(d); chk.append(('MF6 budget file (RCH/DRN_RAVEN)',lc<1e-6)); rec['budgetfile']=lc
        rec['checks']=chk; rec['pass']=bool(all(c[1] for c in chk))
        if rec['pass'] and random.random()<0.4:
            k=random.randint(5,cfg['days']-5); a=d+'_a'; b=d+'_b'
            prep(a,cfg,'1985-10-01',k); ra=run(a); start=(datetime.date(1985,10,1)+datetime.timedelta(days=k)).isoformat()
            prep(b,cfg,start,cfg['days']-k,a+'/out/solution.rvc'); rb=run(b)
            if ('Successful' in ra) and ('Successful' in rb):
                C=pd.read_csv(d+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]']); B=pd.read_csv(b+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]'])
                m=C.merge(B,on='date',suffixes=('_c','_r'))
                dev=max(((m[x+'_c']-m[x+'_r']).abs()/(m[x+'_c'].abs().max()+1e3)).max() for x in ['seepage to Raven [m3/d]','storage release [m3/d]','MF6 recharge [m3/d]','other model packages [m3/d]'])
                rec['restart']=(k,float(dev)); rec['pass']=bool(rec['pass'] and dev<1e-4)
            else: rec['restart']=(k,'FAILED '+ra[:40]+' / '+rb[:40]); rec['pass']=False
    else: rec['pass']=False
    res.append(rec); json.dump(res,open(X+'/fz_results_%d.json'%seed,'w'),indent=1,default=str)
    print(name,'PASS' if rec['pass'] else 'FAIL',cfg,rec.get('restart',''),'' if rec['pass'] else rec.get('msg',''),flush=True)
