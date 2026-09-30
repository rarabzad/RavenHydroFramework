"""Automated audit of a Raven-MODFLOW 6 coupled run directory (<run>/out)."""
import pandas as pd, numpy as np, sys, os, re
def rd(f):
    d=pd.read_csv(f); d.columns=[c.strip() for c in d.columns]; return d
def audit(run,dt=1.0):
    o=(run+'/') if os.path.exists(run+'/GWBudget.csv') else run+'/out/'; R={}
    b=rd(o+'GWBudget.csv')
    tt=b['time [d]'].dropna().values
    if len(tt)>1: dt=float(sorted(abs(tt[1:]-tt[:-1]))[len(tt)//2])
    R['steps']=len(b); R['rows_complete']=int(b.notna().all(axis=1).sum())
    b=b.dropna()
    R['time step [d]']=dt; R['dates_contiguous']=bool((abs(b['time [d]'].diff().dropna()-dt)<1e-6).all()) and b['date'].is_unique
    rr=b['Raven recharge [m3/d]'].sum()*dt; mr=b['MF6 recharge [m3/d]'].sum()*dt
    st=b['recharge delay storage [m3]'].iloc[-1] if 'recharge delay storage [m3]' in b else 0.0
    chdv=b['recharge onto CHD cells [m3/d]'] if 'recharge onto CHD cells [m3/d]' in b else 0*b['MF6 recharge [m3/d]']
    s0=(st if False else 0.0)
    if 'recharge delay storage [m3]' in b:
        s0=b['recharge delay storage [m3]'].iloc[0]-(b['Raven recharge [m3/d]'].iloc[0]-b['MF6 recharge [m3/d]'].iloc[0]-chdv.iloc[0])*dt
    R['recharge: (Raven - MF6 - dStore - toCHD)/Raven']=(rr-mr-(st-s0)-chdv.sum()*dt)/max(rr,1)
    R['recharge delay store at end [m3]']=st
    flows=b['MF6 recharge [m3/d]']+b['storage release [m3/d]']-b['seepage to Raven [m3/d]']+b['river leakage to aquifer [m3/d]']\
          -b['aquifer discharge to rivers [m3/d]']-b['water-table ET [m3/d]']+b['wells actual [m3/d]']+b['GHB [m3/d]']+b['CHD [m3/d]']
    if 'reservoir seepage delivered [m3/d]' in b: flows=flows+b['reservoir seepage delivered [m3/d]']; gross_res=b['reservoir seepage delivered [m3/d]'].abs()
    else: gross_res=0*flows
    if 'surface-store leakage to aquifer [m3/d]' in b: flows=flows+b['surface-store leakage to aquifer [m3/d]']-b['aquifer discharge to surface stores [m3/d]']
    gross=(b['MF6 recharge [m3/d]'].abs()+b['seepage to Raven [m3/d]']+b['river leakage to aquifer [m3/d]']+b['aquifer discharge to rivers [m3/d]']
           +b['water-table ET [m3/d]']+b['wells actual [m3/d]'].abs()+b['GHB [m3/d]'].abs()+b['CHD [m3/d]'].abs()+b['storage release [m3/d]'].abs())
    gross=gross+gross_res
    if 'other model packages [m3/d]' in b: flows=flows+b['other model packages [m3/d]']; gross=gross+b['other model packages [m3/d]'].abs()  #existing-model mode: the model's own packages
    if 'surface-store leakage to aquifer [m3/d]' in b:
        gross=gross+b['surface-store leakage to aquifer [m3/d]']+b['aquifer discharge to surface stores [m3/d]']
        R['surface-store: applied - (discharge-leakage) (rel)']=(b['surface-store exchange applied to Raven [m3]']-(b['aquifer discharge to surface stores [m3/d]']-b['surface-store leakage to aquifer [m3/d]'])*dt).abs().sum()/max((b['surface-store leakage to aquifer [m3/d]']+b['aquifer discharge to surface stores [m3/d]']).sum()*dt,1)
    R['MF6 cumulative residual / gross flow']=flows.sum()/gross.sum()
    if 'reservoir seepage awaiting delivery [m3]' in b and b['reservoir seepage sent [m3/d]'].abs().sum()>0:
        #reservoir seepage reaches MODFLOW one step later: sent(t)*dt = awaiting(t-1); MODFLOW takes what is sent unless the cell is dry
        import numpy as _np
        tt=b['time [d]'].values; dts=_np.diff(_np.concatenate([[0.0],tt])); sv=b['reservoir seepage sent [m3/d]'].values*dts
        aw=b['reservoir seepage awaiting delivery [m3]'].values; dl=b['reservoir seepage delivered [m3/d]'].values*dts
        R['reservoir: sent(t) - awaiting(t-1) (rel)']=float(_np.max(_np.abs(sv[1:]-aw[:-1]))/max(_np.max(_np.abs(sv)),1e-12)) if len(sv)>1 else 0.0
        R['reservoir: delivered / sent']=float(dl.sum()/sv.sum()) if sv.sum()!=0 else 1.0
        if 'reservoir gain not supplied by MODFLOW [m3]' in b:
            shv=b['reservoir gain not supplied by MODFLOW [m3]'].fillna(0).values   #gains MODFLOW could not supply, taken back from the reach
            R['reservoir: (delivered + not supplied) / sent']=float((dl.sum()+shv.sum())/sv.sum()) if sv.sum()!=0 else 1.0
    R['MF6 worst step error [%]']=b['MF6 balance error [%]'].abs().max()
    R['unconverged steps']=int((b['converged']==0).sum()); R['max outer iterations']=int(b['outer iterations'].max())
    R['seepage returned vs MF6 seepage (rel)']=abs(b['returned to Raven [m3]'].sum()-b['seepage to Raven [m3/d]'].sum()*dt)/max(b['seepage to Raven [m3/d]'].sum()*dt,1)
    net=(b['aquifer discharge to rivers [m3/d]']-b['river leakage to aquifer [m3/d]']).sum()*dt
    cc='river loss carried to next step [m3]' if 'river loss carried to next step [m3]' in b else 'river loss deficit [m3]'
    if cc.startswith('river loss carried'):
        c=b[cc].fillna(0); cp=c.shift(1).fillna(0); n=(b['aquifer discharge to rivers [m3/d]']-b['river leakage to aquifer [m3/d]'])*dt; ok=b['converged'].astype(int)>=1
        sh=b['reservoir gain not supplied by MODFLOW [m3]'].fillna(0) if 'reservoir gain not supplied by MODFLOW [m3]' in b else 0.0   #taken from the reach in the same step
        r=(b['river exchange applied to Raven [m3]']-n-c+cp+sh)[ok]; R['river: applied - net MF6 - carry change (rel)']=r.sum()/max(n.abs().sum(),1)
        R['safeguard-limited steps']=int((~ok).sum())
        R['river loss carried at end [m3]']=float(c.iloc[-1]); R['river loss carried max [m3]']=float(c.max())
    else:
        R['river: applied - net MF6 - deficit (rel)']=(b['river exchange applied to Raven [m3]'].sum()-net-b[cc].sum())/max(abs(net),b['river leakage to aquifer [m3/d]'].sum()*dt,1)
        R['river loss deficit [m3]']=b[cc].sum()
    if b['wells specified [m3/d]'].abs().sum()>0:
        R['wells delivered / specified']=b['wells actual [m3/d]'].sum()/b['wells specified [m3/d]'].sum()
    if b['CHD [m3/d]'].abs().sum()>0: R['CHD mean flow [m3/d]']=b['CHD [m3/d]'].mean()
    w=rd(o+'WatershedStorage.csv'); R['Raven MB error max [mm]']=w['MB Error [mm]'].abs().max()
    h=rd(o+'GWHRUState.csv'); H=h[[c for c in h.columns if c.startswith('head_')]].values; D=h[[c for c in h.columns if c.startswith('wtdepth_')]].values
    R['HRU heads finite']=bool(np.isfinite(H).all()); R['max water table above land (HRU mean) [m]']=float(-D.min())
    if os.path.exists(o+'GWHeads.csv'):
        g=rd(o+'GWHeads.csv'); R['GWHeads rows (incl. t=0)']=len(g); R['GWHeads starts at t=0']=bool(g['time [d]'].iloc[0]==0)
    # Raven_errors.txt: in the output folder, or (ensemble members, runs audited from their case folder) next to it
    up=os.path.dirname(o.rstrip('/')); cand=[o+'Raven_errors.txt',up+'/Raven_errors.txt',up+'/out/Raven_errors.txt',up+'/output/Raven_errors.txt']
    ef=[f for f in cand if os.path.exists(f)]
    if ef: e=open(ef[0]).read(); R['unconverged-guard warnings']=e.count('unconverged step'); R['ERROR lines']=len(re.findall(r'^ERROR',e,re.M))
    return R
if __name__=='__main__':
    if len(sys.argv)<2: sys.exit('usage: python3 tools_audit.py OUTPUT_FOLDER [OUTPUT_FOLDER ...]')
    for r in sys.argv[1:]:
        if not (os.path.exists(r+'/GWBudget.csv') or os.path.exists(r+'/out/GWBudget.csv')):
            sys.exit('%s: no GWBudget.csv - the run did not finish (see Raven_errors.txt there)' % r)
    rows={r.rstrip('/').split('/')[-1]:audit(r) for r in sys.argv[1:]}
    df=pd.DataFrame(rows); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',60)
    print(df.to_string(float_format=lambda v:'%.3g'%v))
