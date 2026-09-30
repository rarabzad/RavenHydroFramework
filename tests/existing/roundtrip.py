"""Round trip: a Raven-generated Liard model (no rivers) re-imported as an existing model must reproduce the original run."""
import os,sys; sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
import os,shutil,subprocess,re,sys
audit=tc.audit; import pandas as pd
tc.need_mf6(); X=tc.work('extmodel'); os.makedirs(X,exist_ok=True); B=tc.LIARD; ENV=tc.ENV; RV=tc.RAVEN
def base(d):
    shutil.rmtree(d,ignore_errors=True); os.makedirs(d+'/out')
    for f in os.listdir(B):
        p=B+'/'+f
        if os.path.isfile(p) and f.split('.')[-1] in ('rvi','rvh','rvp','rvt','rvc','geojson','asc'): shutil.copy(p,d)
        elif f in ('data_obs','data_forcing'):
            try: os.symlink(p,d+'/'+f)
            except OSError: shutil.copytree(p,d+'/'+f)
    t=open(d+'/Liard.rvi').read(); t=re.sub(r':Duration\s+\S+(\s+\S+)?',':Duration 365',t); open(d+'/Liard.rvi','w').write(t)
# 1. generated run without rivers
g=X+'/rt_generated'; base(g)
rvg=open(B+'/Liard.rvg').read(); rvg='\n'.join(l for l in rvg.split('\n') if not l.startswith(':RiverGeometry')); open(g+'/Liard.rvg','w').write(rvg)
subprocess.run([RV,'Liard','-o','out/'],cwd=g,capture_output=True,env=ENV)
# 2. its MODFLOW model, without Raven's packages, becomes the existing model
m=X+'/rt_model'; shutil.rmtree(m,ignore_errors=True); shutil.copytree(g+'/out/mf6',m)
for f in ('gwf.rch','gwf.drn','linkage.cache','mfsim.lst','gwf.lst','gwf.hds','gwf.cbc'):
    if os.path.exists(m+'/'+f): os.remove(m+'/'+f)
nam=open(m+'/gwf.nam').read(); nam='\n'.join(l for l in nam.split('\n') if 'RCH_RAVEN' not in l and 'DRN_RAVEN' not in l); open(m+'/gwf.nam','w').write(nam)
e=X+'/rt_external'; base(e); shutil.copy(g+'/out/mf6/grid_cells.geojson',e+'/cells.geojson')
open(e+'/Liard.rvg','w').write(''':HRUGeometry       Liard_HRUs_standin.geojson  HRU_ID
:MinCellCoverage   0.25
:MF6Simulation     %s/mfsim.nam
:MF6GridFile       cells.geojson cell
:MF6StartDate      1985-10-01
:SeepageToRaven    ON
'''%m)
r=subprocess.run([RV,'Liard','-o','out/'],cwd=e,capture_output=True,text=True,env=ENV)
print('existing-model run:',(re.findall(r'Exiting Gracefully: (.*)',r.stdout + r.stderr) or ['?'])[-1])
a=pd.read_csv(g+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]']); b=pd.read_csv(e+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]'])
for c in ['Raven recharge [m3/d]','MF6 recharge [m3/d]','seepage to Raven [m3/d]','storage release [m3/d]']:
    print('  %-26s max relative difference %.1e'%(c,(a[c]-b[c]).abs().max()/a[c].abs().max()))
ha=pd.read_csv(g+'/out/Hydrographs.csv'); hb=pd.read_csv(e+'/out/Hydrographs.csv'); col='SUB_43 [m3/s]'
print('  outlet flow SUB_43              max relative difference %.1e'%((ha[col]-hb[col]).abs().max()/ha[col].abs().max()))
A=audit.audit(e+'/out'); print('  audit: recharge %.0e residual %.1e failed %s'%(abs(A['recharge: (Raven - MF6 - dStore - toCHD)/Raven']),A['MF6 cumulative residual / gross flow'],A['unconverged steps']))
