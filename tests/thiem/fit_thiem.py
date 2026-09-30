import json,glob,os,sys,math,numpy as np,flopy,re
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
B=5000/(2*math.pi*500.0); lat0=55.0; lon0=-120.0; k=111194.9
out=[]
for d in sorted(glob.glob(os.path.join(tc.WORK,'thiem','v_*'))):  # runs of run_thiem.py
    g=json.load(open(d+'/out/mf6/grid_cells.geojson'))
    cen=[]; size=[]
    for f in g['features']:
        c=np.array(f['geometry']['coordinates'][0][:-1]); x=(c[:,0]-lon0)*math.cos(math.radians(lat0))*k; y=(c[:,1]-lat0)*k
        A=0.5*np.sum(x*np.roll(y,-1)-np.roll(x,-1)*y); cx=np.sum((x+np.roll(x,-1))*(x*np.roll(y,-1)-np.roll(x,-1)*y))/(6*A); cy=np.sum((y+np.roll(y,-1))*(x*np.roll(y,-1)-np.roll(x,-1)*y))/(6*A)
        cen.append((cx,cy)); size.append(math.sqrt(abs(A)))
    cen=np.array(cen); size=np.array(size); n=len(cen)
    hf=flopy.utils.HeadFile(d+'/out/mf6/gwf.hds',text='head')
    h=np.array(hf.get_data(totim=hf.get_times()[-1])).ravel()[:n]
    act=(h>-1e20)&(h<1e20)&(np.abs(h)<1e5)
    iw=np.argmin(np.where(act,h,1e30)); x0,y0=cen[iw]
    r=np.hypot(cen[:,0]-x0,cen[:,1]-y0)
    sel=act&(r>=4*size)&(r<=3000)&(r>0)
    X=np.vstack([np.ones(sel.sum()),np.log(r[sel])]).T; coef,res,_,_=np.linalg.lstsq(X,h[sel],rcond=None)
    rms=math.sqrt(np.mean((X@coef-h[sel])**2))
    name=os.path.basename(d)[2:]
    out.append((name,coef[1]/B-1,rms,sel.sum()))
    print('%-40s slope error %+6.2f%%   fit RMS %.4f m   cells used %d'%(name,100*(coef[1]/B-1),rms,sel.sum()))
