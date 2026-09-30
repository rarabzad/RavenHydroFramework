import json,numpy as np,math,flopy,sys
B=5000/(2*math.pi*500.0); lat0=55.0; lon0=-120.0; k=111194.9
for d in sys.argv[1:]:
    g=json.load(open(d+'/out/mf6/grid_cells.geojson')); C=[];S=[]
    for f in g['features']:
        c=np.array(f['geometry']['coordinates'][0][:-1]); x=(c[:,0]-lon0)*math.cos(math.radians(lat0))*k; y=(c[:,1]-lat0)*k
        C.append((x.mean(),y.mean())); S.append(math.sqrt(f['properties']['area_m2']))
    C=np.array(C); S=np.array(S); h=np.array(flopy.utils.HeadFile(d+'/out/mf6/gwf.hds').get_alldata()[-1]).ravel()[:len(C)]
    act=np.abs(h)<1e5; iw=np.argmin(np.where(act,h,1e30)); v=C-C[iw]; r=np.hypot(*v.T); ang=np.degrees(np.arctan2(v[:,1],v[:,0]))
    out=[]
    for a0 in (-45,45,135,225):
        a=(ang-a0)%360; m=act&(r>=4*S)&(r<=3000)&(a<90)
        X=np.vstack([np.ones(m.sum()),np.log(r[m])]).T; c,_,_,_=np.linalg.lstsq(X,h[m],rcond=None); out.append(100*(c[1]/B-1))
    print('%-34s slope error by direction  E %+5.2f%%  N %+5.2f%%  W %+5.2f%%  S %+5.2f%%'%(d[2:],*out))
