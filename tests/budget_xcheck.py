# Copyright (c) 2026 Rezgar Arabzadeh -- Raven-MODFLOW 6 groundwater coupling
# SPDX-License-Identifier: Artistic-2.0
import pandas as pd, sys, re
# Compares the coupler's GWBudget.csv with MODFLOW's own list-file budget, package by package (by package name, as
# MODFLOW prints it), for the last step.
# Usage: python3 budget_xcheck.py OUTPUT_FOLDER [OUTPUT_FOLDER ...]   (Raven output folders of generated models)
if len(sys.argv)<2: sys.exit(__doc__ if __doc__ else 'usage: budget_xcheck.py OUTPUT_FOLDER [...]')
runs={d:d for d in sys.argv[1:]}
pairs=[('RCH_RAVEN',lambda b:b['MF6 recharge [m3/d]']),('DRN_RAVEN',lambda b:-b['seepage to Raven [m3/d]']),
       ('RIV_RAVEN',lambda b:b['river leakage to aquifer [m3/d]']-b['aquifer discharge to rivers [m3/d]']),('EVT_RAVEN',lambda b:-b['water-table ET [m3/d]']),
       ('WEL_RAVEN',lambda b:b['wells actual [m3/d]']),('GHB_RAVEN',lambda b:b['GHB [m3/d]']),('CHD_RAVEN',lambda b:b['CHD [m3/d]']),
       ('SWX_RAVEN',lambda b:b['surface-store leakage to aquifer [m3/d]']-b['aquifer discharge to surface stores [m3/d]']),
       ('RES_RAVEN',lambda b:b['reservoir seepage delivered [m3/d]'])]
def list_budget(lst):
    """rates of the last volume budget in a MODFLOW 6 list file, by package name: {name: (in, out)} [L3/T]"""
    txt=open(lst).read(); k=txt.rfind('VOLUME BUDGET FOR ENTIRE MODEL'); blk=txt[k:txt.find('TOTAL OUT',k)+200]
    rates={}; side=None
    for l in blk.split('\n'):
        s=l.strip()
        if s.startswith('IN:'): side=0
        elif s.startswith('OUT:'): side=1
        m=re.match(r'\s*\S+\s*=\s*[-\d.Ee+]+\s+\S+\s*=\s*([-\d.Ee+]+)\s+(\S+)\s*$',l)
        if m and side is not None and 'TOTAL' not in l:
            r=rates.setdefault(m.group(2).upper(),[0.0,0.0]); r[side]+=float(m.group(1))
    return rates
worst=0; ncmp=0
for name,d in runs.items():
    L=list_budget(d+'/mf6/gwf.lst')
    b=pd.read_csv(d+'/GWBudget.csv').rename(columns=lambda c:c.strip()).dropna(subset=['MF6 recharge [m3/d]']).iloc[-1]
    out=[]
    for pk,f in pairs:
        if pk not in L: continue
        mf=L[pk][0]-L[pk][1]; ours=f(b); rel=abs(mf-ours)/max(abs(mf),1.0); worst=max(worst,rel); ncmp+=1
        out.append('%s %.6g/%.6g (%.1e)'%(pk.replace('_RAVEN',''),ours,mf,rel))
    print('%-14s'%name,' | '.join(out))
if ncmp==0: sys.exit('no package of the list-file budget could be matched: nothing compared')
print('%d package budgets compared; WORST relative difference coupler vs MODFLOW list-file budget: %.2e'%(ncmp,worst))
