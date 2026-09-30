"""Liard Raven model coupled to the synthetic existing MF6 model. Usage: setup_case.py DEST MODELDIR [extra rvg lines]"""
import os,sys; sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
import os,shutil,sys,re
B=tc.LIARD; dest=sys.argv[1]; md=sys.argv[2]; extra=sys.argv[3] if len(sys.argv)>3 else ''
shutil.rmtree(dest,ignore_errors=True); os.makedirs(dest)
for f in os.listdir(B):
    p=B+'/'+f
    if os.path.isfile(p) and f.split('.')[-1] in ('rvi','rvh','rvp','rvt','rvc','geojson','asc'): shutil.copy(p,dest)
    elif f in ('data_obs','data_forcing'):
        try: os.symlink(p,dest+'/'+f)
        except OSError: shutil.copytree(p,dest+'/'+f)   # (no symlinks, e.g. Windows without privileges)
shutil.copy(md+'/ext_cells.geojson',dest)
t=open(dest+'/Liard.rvi').read(); t=re.sub(r':Duration\s+\S+(\s+\S+)?',':Duration 731',t); open(dest+'/Liard.rvi','w').write(t)
open(dest+'/Liard.rvg','w').write(f''':HRUGeometry       Liard_HRUs_standin.geojson  HRU_ID
:MinCellCoverage   0.25
:MF6Simulation     {md}/mfsim.nam
:MF6GridFile       ext_cells.geojson CELL_ID
:MF6StartDate      1985-10-01
:RavenTakesOver    RCH_MODEL EVT_MODEL
:StreamPackage     RIV_STREAMS AUX SUBBASIN
:SeepageToRaven    ON
{extra}''')
os.makedirs(dest+'/out')
print('case ready:',dest)
