import os,sys; sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
import os,shutil,subprocess,re,sys,datetime,pandas as pd
audit=tc.audit
tc.need_mf6(); X=tc.ext_models(); ENV=tc.ENV; RV=tc.RAVEN
def prep(d,start,dur,rvc=None,crs=True):
    subprocess.run([sys.executable,tc.SETUP_EXT,d,X+'/md'],capture_output=True)
    t=open(d+'/Liard.rvi').read(); t=re.sub(r':StartDate\s+\S+\s+\S+',':StartDate %s 00:00:00'%start,t); t=re.sub(r':Duration\s+\S+',':Duration %d'%dur,t); open(d+'/Liard.rvi','w').write(t)
    if crs: g=open(d+'/Liard.rvg').read(); open(d+'/Liard.rvg','w').write(re.sub(r':MF6GridFile.*',':MF6CRS            EPSG:32610',g))
    if rvc: shutil.copy(rvc,d+'/Liard.rvc')
def run(d):
    r=subprocess.run([RV,'Liard','-o','out/'],cwd=d,capture_output=True,text=True,env=ENV); return (re.findall(r'Exiting Gracefully: (.*)',r.stdout + r.stderr) or ['?'])[-1]
c=X+'/rs_cont'; prep(c,'1985-10-01',200); print('continuous:',run(c))
C=pd.read_csv(c+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]'])
for split in (30,150):
    a=X+'/rs_a%d'%split; b=X+'/rs_b%d'%split
    prep(a,'1985-10-01',split); ra=run(a)
    start=(datetime.date(1985,10,1)+datetime.timedelta(days=split)).isoformat()
    prep(b,start,200-split,a+'/out/solution.rvc'); rb=run(b)
    if 'Successful' not in ra+rb: print('split',split,'FAILED',ra,rb); continue
    B=pd.read_csv(b+'/out/GWBudget.csv').dropna(subset=['Raven recharge [m3/d]']); m=C.merge(B,on='date',suffixes=('_c','_r'))
    dev=max(((m[x+'_c']-m[x+'_r']).abs()/(m[x+'_c'].abs().max()+1e3)).max() for x in ['seepage to Raven [m3/d]','storage release [m3/d]','aquifer discharge to rivers [m3/d]','MF6 recharge [m3/d]','other model packages [m3/d]'])
    w=m['other model packages [m3/d]_c']; jump=w.diff().abs().max()
    A=audit.audit(b+'/out')
    print('split at day %3d (%s): %d days compared, max deviation %.1e; wells/boundaries flow change seen in restarted run: %s; audit recharge %.0e residual %.1e failed %s'%(
        split,start,len(m),dev,'yes' if m['other model packages [m3/d]_r'].diff().abs().max()>0.5*jump and jump>1000 else 'no',abs(A['recharge: (Raven - MF6 - dStore - toCHD)/Raven']),A['MF6 cumulative residual / gross flow'],A['unconverged steps']))
    print('   renumbered first period block note:', [l.strip() for l in open(b+'/out/mf6/gwf_liard.wel').read().split('\n') if 'in effect' in l][:1])
