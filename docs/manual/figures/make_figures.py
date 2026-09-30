import json, re, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, flopy
from matplotlib.patches import Polygon as MPoly
from matplotlib.collections import PatchCollection
plt.rcParams.update({'font.size':9,'figure.dpi':150})
G='/home/claude/liard_gw2/'; RUN=G+'out_final/'; M=RUN+'mf6/'
dis=open(M+'gwf.dis').read(); lat0,lon0=map(float,re.search(r'lat=([-\d.]+) lon=([-\d.]+)',dis).groups())
x0=float(re.search(r'XORIGIN (\S+)',dis).group(1)); y0=float(re.search(r'YORIGIN (\S+)',dis).group(1))
sim=flopy.mf6.MFSimulation.load(sim_ws=M,verbosity_level=0,load_only=['dis']); d=sim.get_model().dis
top=d.top.array; idom=d.idomain.array; nr,nc=top.shape; cs=2000.; R=6371007.181; p0=np.radians(lat0)
def inv(x,y):
    rho=np.hypot(x,y); c=2*np.arcsin(rho/(2*R)); phi=np.arcsin(np.cos(c)*np.sin(p0)+y*np.sin(c)*np.cos(p0)/rho)
    lam=np.arctan2(x*np.sin(c),rho*np.cos(p0)*np.cos(c)-y*np.sin(p0)*np.sin(c)); return np.degrees(lam)+lon0,np.degrees(phi)
def cellpoly(r,c):
    xs=[x0+c*cs,x0+(c+1)*cs,x0+(c+1)*cs,x0+c*cs]; ys=[y0+(nr-r)*cs,y0+(nr-r)*cs,y0+(nr-r-1)*cs,y0+(nr-r-1)*cs]
    return np.array([inv(a,b) for a,b in zip(xs,ys)])
# ---- Figure: grid, HRUs, rivers
gj=json.load(open(G+'Liard_HRUs_standin.geojson')); col={'VALLEY_ALLUVIUM':'#9ecae1','UPLAND_TILL':'#fdae6b','[NONE]':'#e0e0e0'}
fig,ax=plt.subplots(figsize=(7.2,5.2))
for f in gj['features']:
    g=f['geometry']; polys=[g['coordinates']] if g['type']=='Polygon' else g['coordinates']
    for P in polys: ax.add_patch(MPoly(np.array(P[0]),closed=True,fc=col[f['properties']['AQUIFER_PROFILE']],ec='white',lw=0.4))
act=[MPoly(cellpoly(r,c),closed=True) for r in range(nr) for c in range(nc) if idom[0,r,c]>0]
ax.add_collection(PatchCollection(act,facecolor='none',edgecolor='k',lw=0.15,alpha=0.5))
riv=[tuple(map(int,l.split()[:3])) for l in open(M+'gwf.riv') if re.match(r'\s+\d+ \d+ \d+ ',l)]
ax.add_collection(PatchCollection([MPoly(cellpoly(r-1,c-1),closed=True) for _,r,c in riv],facecolor='#08519c',edgecolor='none',alpha=0.8))
for f in json.load(open(G+'Liard_rivers.geojson'))['features']:
    a=np.array(f['geometry']['coordinates']); ax.plot(a[:,0],a[:,1],color='#2171b5',lw=0.6)
allc=np.vstack([cellpoly(r,c) for r in (0,nr-1) for c in (0,nc-1)]); ax.set_xlim(allc[:,0].min(),allc[:,0].max()); ax.set_ylim(allc[:,1].min(),allc[:,1].max())
from matplotlib.patches import Patch
ax.legend(handles=[Patch(fc=col['VALLEY_ALLUVIUM'],label='VALLEY_ALLUVIUM HRUs'),Patch(fc=col['UPLAND_TILL'],label='UPLAND_TILL HRUs'),Patch(fc=col['[NONE]'],label='no groundwater (glacier)'),
                   Patch(fc='none',ec='k',label='active column (2 km)'),Patch(fc='#08519c',label='river cell')],loc='lower left',fontsize=7)
ax.set_xlabel('longitude'); ax.set_ylabel('latitude'); ax.set_aspect(1/np.cos(np.radians(lat0))); fig.tight_layout(); fig.savefig('fig_grid_liard.png'); plt.close()
# ---- Figure: water-table depth map (end of run)
hd=flopy.utils.HeadFile(M+'gwf.hds'); H=hd.get_data(); dep=np.where(idom[0]>0,top-H[0],np.nan)
fig,ax=plt.subplots(figsize=(6,4.6)); im=ax.imshow(np.clip(dep,-1,20),cmap='viridis_r',extent=[x0/1e3,(x0+nc*cs)/1e3,y0/1e3,(y0+nr*cs)/1e3])
plt.colorbar(im,ax=ax,label='water-table depth below land surface [m]'); ax.set_xlabel('x [km] (equal-area projection)'); ax.set_ylabel('y [km]')
ax.set_title('Liard example: water-table depth, last saved head (%s)'%str(hd.get_times()[-1])+' d',fontsize=8); fig.tight_layout(); fig.savefig('fig_wtdepth.png'); plt.close()
# ---- Figure: groundwater budget (monthly means)
b=pd.read_csv(RUN+'GWBudget.csv').rename(columns=lambda c:c.strip()).dropna(); b['date']=pd.to_datetime(b['date']); m=b.set_index('date').resample('MS').mean(numeric_only=True)/1e6
fig,ax=plt.subplots(2,1,figsize=(7.2,5),sharex=True)
ax[0].plot(m.index,m['MF6 recharge [m3/d]'],label='recharge'); ax[0].plot(m.index,m['seepage to Raven [m3/d]'],label='seepage to Raven')
ax[0].plot(m.index,m['aquifer discharge to rivers [m3/d]'],label='discharge to rivers'); ax[0].plot(m.index,m['river leakage to aquifer [m3/d]'],label='river leakage')
ax[0].plot(m.index,m['storage release [m3/d]'],label='storage release',lw=0.8); ax[0].set_ylabel('10$^6$ m$^3$/d (monthly mean)'); ax[0].legend(fontsize=7,ncol=3)
ax[1].semilogy(b['date'],b['MF6 balance error [%]'].abs().clip(lower=1e-9),lw=0.4); ax[1].set_ylabel('|balance error| [% of gross flow]'); ax[1].axhline(0.1,color='r',lw=0.8,ls='--')
fig.tight_layout(); fig.savefig('fig_budget.png'); plt.close()
# ---- Figure: hydrographs subbasin 43 (last 3 years)
q0=pd.read_csv('/home/claude/liard_test/out_base/Hydrographs.csv'); q1=pd.read_csv(RUN+'Hydrographs.csv')
c=[x for x in q1.columns if x.startswith('SUB_43') and 'observed' not in x][0]; o=[x for x in q1.columns if x.startswith('SUB_43') and 'observed' in x]
t=pd.to_datetime(q1['date']); s=t>='2002-10-01'
fig,ax=plt.subplots(figsize=(7.2,3))
if o: ax.plot(t[s],pd.to_numeric(q1[o[0]],errors='coerce')[s],'k.',ms=1.5,label='observed')
ax.plot(t[s],q0[c][s],lw=0.8,label='Raven (uncoupled)'); ax.plot(t[s],q1[c][s],lw=0.8,label='Raven + MODFLOW 6')
ax.set_ylabel('outflow SUB_43 [m$^3$/s]'); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig('fig_hydrograph.png'); plt.close()
# ---- Figure: pumping test (Phase 3 runs R1, R2)
P='/home/claude/p3/'; h1=pd.read_csv(P+'R1/out/GWHeads.csv').rename(columns=lambda c:c.strip()); h2=pd.read_csv(P+'R2/out/GWHeads.csv').rename(columns=lambda c:c.strip())
b1=pd.read_csv(P+'R1/out/GWBudget.csv').rename(columns=lambda c:c.strip()).dropna(); b2=pd.read_csv(P+'R2/out/GWBudget.csv').rename(columns=lambda c:c.strip()).dropna()
tt=pd.to_datetime(h1['date']); fig,ax=plt.subplots(1,2,figsize=(7.4,2.9))
for cc,lab in zip(h1.columns[2:],['OW_NEAR (pumped cell)','OW_FAR (12 km)']): ax[0].plot(tt,h1[cc]-h2[cc],label=lab)
ax[0].axvline(pd.Timestamp('1987-01-01'),color='k',ls=':',lw=0.8); ax[0].set_ylabel('drawdown [m]'); ax[0].legend(fontsize=7)
net=lambda b:b['seepage to Raven [m3/d]']+b['aquifer discharge to rivers [m3/d]']-b['river leakage to aquifer [m3/d]']
t2=pd.to_datetime(b2['date']); cap=(net(b1)-net(b2)).rolling(30).mean(); st=(b2['storage release [m3/d]']-b1['storage release [m3/d]']).rolling(30).mean()
ax[1].plot(t2,-b2['wells actual [m3/d]']/1e3,'k',lw=0.8,label='pumping (actual)'); ax[1].plot(t2,cap/1e3,label='captured discharge'); ax[1].plot(t2,st/1e3,label='storage')
ax[1].set_ylabel('10$^3$ m$^3$/d (30-d mean)'); ax[1].legend(fontsize=7); fig.autofmt_xdate(); fig.tight_layout(); fig.savefig('fig_pumping.png'); plt.close()
# ---- Figure: connectivity (drying test)
cn=pd.read_csv('/home/claude/liard_dry/out_final/GWConnectivity.csv').rename(columns=lambda c:c.strip()); tc=pd.to_datetime(cn['date'])
fig,ax=plt.subplots(figsize=(7.2,2.8)); ax.plot(tc,cn['connected wet clusters'],'k',lw=0.8); ax.set_ylabel('connected saturated clusters')
a2=ax.twinx(); a2.plot(tc,cn['wet fraction VALLEY_ALLUVIUM'],color='tab:blue',lw=0.8); a2.set_ylabel('wet fraction (valley)',color='tab:blue'); fig.tight_layout(); fig.savefig('fig_connectivity.png'); plt.close()
# ---- Figure: calibration convergence
L=open('/home/claude/calib/calib_log.txt').read().split('\n'); objs=[]; best=[]
for l in L:
    mm=re.match(r'eval (\d+).*obj ([\d.]+)',l)
    if mm: objs.append(float(mm.group(2)))
bb=np.minimum.accumulate(objs); fig,ax=plt.subplots(figsize=(6,2.8)); ax.semilogy(objs,'.',ms=3,color='gray',label='evaluation'); ax.semilogy(bb,'k',label='best so far')
ax.axhline(0.0802,color='tab:red',ls='--',lw=0.8,label='objective at true parameters'); ax.axvline(30,color='k',ls=':',lw=0.8)
ax.set_xlabel('model evaluation (DDS)'); ax.set_ylabel('sum of head RMSE [m]'); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig('fig_calibration.png'); plt.close()
print('figures done')
