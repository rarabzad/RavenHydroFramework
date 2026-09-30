import os,shutil,subprocess,re,math,pandas as pd,sys
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(HERE)); import testconfig as tc; audit=tc.audit
T=500.0; Q=5000.0; a=Q/(2*math.pi*T)*math.log(2)
tc.need_mf6(); ENV=tc.ENV; RV=tc.RAVEN; OUT=os.path.dirname(tc.work('thiem','x'))
base=open(HERE+'/base.rvg').read()
V=[('regular 1000 m',':GridCellSize 1000',''),('regular 500 m',':GridCellSize 500',''),('regular 250 m',':GridCellSize 250',''),('regular 125 m',':GridCellSize 125',''),
   ('regular 1000 m, variable spacing to 125 m',':GridCellSize 1000',':GridRefinement WELLS 125\n'),
   ('regular 250 m, rotated 17 deg',':GridCellSize 250',':GridRotation 17\n'),
   ('quadtree 1000 m to 62.5 m',':GridCellSize 1000',':GridType QUADTREE\n:GridRefinement WELLS 62.5\n'),
   ('HRU mesh 1000 m, 125 m at wells',':GridCellSize 1000',':GridType HRU_MESH\n:GridRefinement WELLS 125\n'),
   ('quadtree 500 m to 62.5 m',':GridCellSize 500',':GridType QUADTREE\n:GridRefinement WELLS 62.5\n'),
   ('quadtree 250 m to 62.5 m',':GridCellSize 250',':GridType QUADTREE\n:GridRefinement WELLS 62.5\n'),
   ('HRU mesh 500 m',':GridCellSize 500',':GridType HRU_MESH\n:GridRefinement WELLS 125\n'),
   ('HRU mesh 250 m',':GridCellSize 250',':GridType HRU_MESH\n:GridRefinement WELLS 125\n'),
   ('HRU mesh 250 m without XT3D',':GridCellSize 250',':GridType HRU_MESH\n:GridRefinement WELLS 125\n:FlowCorrection NONE\n'),
   ('HRU mesh 250 m with XT3D',':GridCellSize 250',':GridType HRU_MESH\n:GridRefinement WELLS 125\n:FlowCorrection XT3D\n'),
   ('HRU mesh 250 m lattice only',':GridCellSize 250',':GridType HRU_MESH\n'),
   ('regular 250 m via DISU',':GridCellSize 250',':LayerConnection ELEVATION\n'),
   ('imported (quadtree cells)',':GridCellSize 1000',':GridType FILE\n:GridFile qt_cells.geojson cell\n')]
rows=[]
for name,csl,extra in V:
    d=OUT+'/v_'+re.sub(r'[^a-z0-9]+','_',name.lower()).strip('_')
    shutil.rmtree(d,ignore_errors=True); os.makedirs(d+'/out')
    for f in ('Thiem.rvi','Thiem.rvh','Thiem.rvp','Thiem.rvt','Thiem.rvc','hru.geojson','edge.geojson'): shutil.copy(HERE+'/'+f,d)
    if 'imported' in name: shutil.copy(OUT+'/qt_cells.geojson',d)
    open(d+'/Thiem.rvg','w').write(base.replace(':GridCellSize 1000',csl)+extra)
    r=subprocess.run([RV,'Thiem','-o','out/'],cwd=d,capture_output=True,text=True,env=ENV)
    ex=(re.findall(r'Exiting Gracefully: (.*)',r.stdout + r.stderr) or ['?'])[-1]
    if 'Successful' not in ex: print(name,'FAILED',ex); continue
    if 'quadtree' in name: shutil.copy(d+'/out/mf6/grid_cells.geojson',OUT+'/qt_cells.geojson')
    h=pd.read_csv(d+'/out/GWHeads.csv').iloc[-1]
    c=[x for x in h.index if x.startswith('OW')]; hv=dict((int(x.split('OW')[1].split()[0]),h[x]) for x in c)
    e1=(hv[1200]-hv[600])/a-1; e2=(hv[2400]-hv[1200])/a-1
    s=open(d+'/out/GWModelSummary.txt').read(); ncell=re.findall(r'(\d+) (?:cells|cols)',s)
    A=audit.audit(d+'/out')
    rows.append((name,e1,e2,A.get('MF6 cumulative residual / gross flow'),A.get('unconverged steps')))
    print('%-42s dh(600-1200) err %+7.2f%%   dh(1200-2400) err %+7.2f%%   residual %.1e  failed %s'%(name,100*e1,100*e2,A.get('MF6 cumulative residual / gross flow'),A.get('unconverged steps')))
print('analytic head difference per doubling of distance: %.4f m'%a)
