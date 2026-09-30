"""Built-in projections (CGWCRS) against PROJ/pyproj: 2000 random points per coordinate system.
Build first: g++ -std=c++11 -O2 -I../../src crstest.cpp ../../src/GWGeometry.cpp -o crstest   (pyproj is needed only for this test)"""
import numpy as np, pyproj, subprocess
rng=np.random.default_rng(3)
cases=[('EPSG:32610',(-126,-120),(40,72)),('EPSG:32609',(-132,-126),(50,62)),('EPSG:32718',(-78,-72),(-50,-5)),('EPSG:26910',(-126,-120),(45,60)),
       ('EPSG:3005',(-139,-114),(48,60)),('EPSG:5070',(-124,-67),(25,49)),('EPSG:3978',(-141,-52),(42,83)),('EPSG:3347',(-141,-52),(42,83)),('ESRI:102001',(-141,-52),(42,83))]
worst=0
for code,lo,la in cases:
    lon=rng.uniform(*lo,2000); lat=rng.uniform(*la,2000)
    x,y=pyproj.Transformer.from_crs(4326,code,always_xy=True).transform(lon,lat)
    out=subprocess.run(['./crstest',code.replace('ESRI:','EPSG:')],input='\n'.join('%.6f %.6f'%p for p in zip(x,y)),capture_output=True,text=True).stdout.split()
    r=np.array(out,float).reshape(-1,2); e=np.max(np.hypot((r[:,0]-lon)*np.cos(np.radians(lat))*111320,(r[:,1]-lat)*110574)); worst=max(worst,e)
    print('%-12s largest position difference %.2e m'%(code,e))
print('worst %.1e m'%worst)
