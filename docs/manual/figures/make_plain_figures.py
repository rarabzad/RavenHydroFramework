"""fig_placing.pdf: the synthetic Liard model's cells in the model's own coordinates and placed on the map.
Reads the cells of a coupled run (grid_cells.geojson) and the Liard HRUs; needed only to regenerate the figure.
The other plain-language figures are drawn in TikZ (tikz_*.tex; tikz_linking.tex from make_tikz_linking.py)."""
import json, math, re, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt.rcParams.update({'font.family':'TeX Gyre Heros','font.size':9,'axes.titlesize':9.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8})
TEAL='#1D6F5A'; BAND='#8FB0D6'; FILL='#DCE7F4'

def coupled_hrus(rvh):
    t=open(rvh).read(); hdr=[l for l in t.split('\n') if ':Attributes' in l and 'AQUIFER_PROFILE' in l][0].replace(',',' ').split()[1:]
    ia=hdr.index('AQUIFER_PROFILE')+1; ids=set()
    for l in t.split(':HRUs')[1].split(':EndHRUs')[0].split('\n'):
        p=l.replace(',',' ').split()
        if p and p[0].isdigit() and len(p)>ia and p[ia]!='[NONE]': ids.add(int(p[0]))
    return ids

def fig_placing(cells='/home/claude/extmodel/case_crs/out/mf6/grid_cells.geojson',hrus='/home/claude/liard_gw2/Liard_HRUs_standin.geojson',
                rvh='/home/claude/liard_gw2/Liard.rvh',ncol=98,nrow=117,cs=1.5):
    g=json.load(open(cells))['features']; H=json.load(open(hrus))['features']; keep=coupled_hrus(rvh)
    loc=[];geo=[];fc=[]
    for f in g:
        if not f['properties']['active']: continue
        r,c=divmod(f['properties']['cell']-1,ncol)
        x0,y0=c*cs,(nrow-1-r)*cs; loc.append([[x0,y0],[x0+cs,y0],[x0+cs,y0+cs],[x0,y0+cs]])
        geo.append(np.array(f['geometry']['coordinates'][0])[:,:2]); fc.append(BAND if (r%10==0 or c%10==0) else FILL)
    fig,(a,b)=plt.subplots(1,2,figsize=(6.3,3.9),gridspec_kw=dict(wspace=0.38))
    a.add_collection(PolyCollection(loc,facecolors=fc,edgecolors='none'))
    L=np.vstack(loc); a.set_xlim(L[:,0].min()-4,L[:,0].max()+4); a.set_ylim(L[:,1].min()-4,L[:,1].max()+4); a.set_aspect('equal')
    a.set_xlabel('distance along the rows (km)'); a.set_ylabel('distance along the columns (km)')
    a.set_title('(a) As MODFLOW describes it',loc='left')
    b.add_collection(PolyCollection(geo,facecolors=fc,edgecolors='none'))
    for f in H:
        if int(f['properties']['HRU_ID']) not in keep: continue
        gg=f['geometry']; ps=[gg['coordinates']] if gg['type']=='Polygon' else gg['coordinates']
        for p in ps:
            r=np.array(p[0])[:,:2]; b.plot(r[:,0],r[:,1],color=TEAL,lw=0.55)
    G=np.vstack(geo); b.set_xlim(G[:,0].min()-0.08,G[:,0].max()+0.08); b.set_ylim(G[:,1].min()-0.05,G[:,1].max()+0.05)
    b.set_aspect(1/math.cos(math.radians(G[:,1].mean()))); b.set_xlabel('longitude (°)'); b.set_ylabel('latitude (°)')
    b.set_title('(b) Placed on the map by Raven',loc='left')
    fig.legend(handles=[Patch(fc=FILL,ec='none',label='active model cells'),Patch(fc=BAND,ec='none',label='every tenth row and column'),
                        Line2D([0],[0],color=TEAL,lw=1,label='coupled HRUs')],loc='lower center',ncol=3,frameon=False,fontsize=8.5,bbox_to_anchor=(0.5,0.0))
    fig.subplots_adjust(bottom=0.17,top=0.93,left=0.1,right=0.98)
    fig.savefig('fig_placing.pdf'); plt.close(fig)

if __name__=='__main__':
    fig_placing(); print('fig_placing.pdf written')
