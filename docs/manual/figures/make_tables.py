import sys, os; sys.path.insert(0,'/home/claude/p4'); import audit, pandas as pd, numpy as np
runs=[('20-year Liard','/home/claude/liard_gw2/out_final'),('Pumping','/home/claude/p3/R2obs'),('Fixed head, delay','/home/claude/p3/T3'),
      ('Restart','/home/claude/p3/T3b'),('Ensemble','/home/claude/p3/T5/out/ens_2'),('Wetlands','/home/claude/p3/T7'),('Reservoir','/home/claude/p3/T8'),('Soil-less HRUs','/home/claude/p3/T9'),('Drying (extreme)','/home/claude/liard_dry/out_final')]
R={}
for n,p in runs:
    if os.path.exists(p+'/GWBudget.csv') or os.path.exists(p+'/out/GWBudget.csv'):
        try: R[n]=audit.audit(p)
        except Exception as e: print('skip',n,e)
def f(v,kind):
    if v is None or (isinstance(v,float) and np.isnan(v)): return '--'
    if kind=='e': return '0' if v==0 else '$%.1e$'%v
    if kind=='p': return '%.3f'%v
    if kind=='i': return '%d'%v
    if kind=='b': return 'yes' if v else 'no'
    return '%.3g'%v
rows=[('steps','steps','i'),('Recharge balance (rel.)','recharge: (Raven - MF6 - dStore - toCHD)/Raven','e'),('MF6 cumulative residual','MF6 cumulative residual / gross flow','e'),
      ('MF6 worst step [\\%]','MF6 worst step error [%]','p'),('Unconverged steps','unconverged steps','i'),('River bookkeeping (rel.)','river: applied - net MF6 - carry change (rel)','e'),
      ('River loss carried at end [m$^3$]','river loss carried at end [m3]','e'),('Seepage returned (rel.)','seepage returned vs MF6 seepage (rel)','e'),('Surface-store bookkeeping (rel.)','surface-store: applied - (discharge-leakage) (rel)','e'),('Raven balance error [mm]','Raven MB error max [mm]','p'),
      ('Dates contiguous','dates_contiguous','b')]
names=list(R)
import math
def sci(v):
    if v is None or (isinstance(v,float) and math.isnan(v)): return '--'
    if v==0: return '0'
    e=int(math.floor(math.log10(abs(v)))); m=v/10**e
    return '$%s%.1f\\times10^{%d}$'%('-' if m<0 else '',abs(m),e)
def pct(v):
    return '--' if v is None else ('$%s$'%sci(v).strip('$') if abs(v)<0.001 and v!=0 else '%.3g'%v)
cols=[('Steps','steps',lambda v:'{:,}'.format(int(v)).replace(',','\\,')),('Recharge','recharge: (Raven - MF6 - dStore - toCHD)/Raven',sci),
      ('\\shortstack[r]{MODFLOW\\\\residual}','MF6 cumulative residual / gross flow',sci),('\\shortstack[r]{Worst\\\\step [\\%]}','MF6 worst step error [%]',pct),
      ('\\shortstack[r]{Failed\\\\steps}','unconverged steps',lambda v:'%d'%v),('River','river: applied - net MF6 - carry change (rel)',sci),('\\shortstack[r]{Raven\\\\balance [mm]}','Raven MB error max [mm]',lambda v:'%.2f'%v)]
out=[r'\begin{table}[!htbp]\centering',
     r'\caption{Automated audit of the final code (\cmd{tools\_audit.py}). Relative errors are shares of the water moved. In every run the seepage returned to Raven equals MODFLOW\textquotesingle s seepage exactly, the surface-store bookkeeping is exact, no river loss is carried at the end and the output dates are contiguous. The one failed step is the rewetting day of the deliberately extreme drying test (Section~\ref{sec:behaviour}).}\label{tab:audit}',
     r'\footnotesize\setlength{\tabcolsep}{2.8pt}\begin{tabular}{l'+'r'*len(cols)+r'}\toprule\hdr',
     'Run & '+' & '.join(c[0] for c in cols)+r'\\\midrule']
for n in names:
    out.append(n+' & '+' & '.join(fn(R[n].get(key)) if R[n].get(key) is not None else '--' for _,key,fn in cols)+r'\\')
out.append(r'\bottomrule\end{tabular}\end{table}')
open('/home/claude/deliver/3_overleaf_manual/tables/audit_table.tex','w').write('\n'.join(out)+'\n'); print('table runs:',names)
