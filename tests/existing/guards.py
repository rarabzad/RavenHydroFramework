import os,sys; sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
import os,shutil,subprocess,re,json
tc.need_mf6(); X=tc.ext_models(); ENV=tc.ENV; RV=tc.RAVEN
def case(name,rvg_edit=None,model_edit=None,rvi_edit=None,grid_edit=None,expect_ok=False):
    md=X+'/md'
    if model_edit:
        md=X+'/g_model_'+name; shutil.rmtree(md,ignore_errors=True); shutil.copytree(X+'/md',md); model_edit(md)
    d=X+'/g_'+name; subprocess.run([sys.executable,tc.SETUP_EXT,d,md],capture_output=True)
    t=open(d+'/Liard.rvi').read(); t=re.sub(r':Duration\s+\S+',':Duration 31',t)
    if rvi_edit: t=rvi_edit(t)
    open(d+'/Liard.rvi','w').write(t)
    if rvg_edit: g=open(d+'/Liard.rvg').read(); open(d+'/Liard.rvg','w').write(rvg_edit(g))
    if grid_edit: grid_edit(d+'/ext_cells.geojson')
    r=subprocess.run([RV,'Liard','-o','out/'],cwd=d,capture_output=True,text=True,env=ENV,timeout=280)
    msg=(re.findall(r'Exiting Gracefully: (.*)',r.stdout + r.stderr) or ['(no exit message; rc %d)'%r.returncode])[-1]
    ok=('Successful' in msg)==expect_ok
    print('%-4s %-34s %s'%('ok' if ok else 'FAIL',name,msg[:150]))
    return ok
def two_models(md):
    s=open(md+'/mfsim.nam').read(); s=re.sub(r'(?i)(\s+gwf6\s+\S+\s+gwf_liard)',r'\1\n  gwf6  gwf_liard.nam  gwf_copy',s,1); open(md+'/mfsim.nam','w').write(s)
def bad_period(md):
    s=open(md+'/ext.tdis').read(); s=s.replace('1.00000000       1       1.00000000','1.50000000       1       1.00000000',1)
    if '1.50000000' not in s: s=re.sub(r'(BEGIN perioddata\s*\n\s*)1\.0+\S*\s+1\s+1\.0+',r'\g<1>1.5 1 1.0',s,flags=re.I)
    open(md+'/ext.tdis','w').write(s)
def feet_grid(f):
    g=json.load(open(f))
    for ft in g['features']:
        c=ft['geometry']['coordinates'][0]; cx=sum(p[0] for p in c)/len(c); cy=sum(p[1] for p in c)/len(c)
        ft['geometry']['coordinates'][0]=[[cx+(p[0]-cx)*0.3048,cy+(p[1]-cy)*0.3048] for p in c]
    json.dump(g,open(f,'w'))
R=[]
R.append(case('valid_aux'                ,expect_ok=True))
R.append(case('valid_dominant_hru'       ,rvg_edit=lambda g:g.replace(':StreamPackage     RIV_STREAMS AUX SUBBASIN',':StreamPackage     RIV_STREAMS DOMINANT_HRU'),expect_ok=True))
R.append(case('undeclared_evt'           ,rvg_edit=lambda g:g.replace(':RavenTakesOver    RCH_MODEL EVT_MODEL',':RavenTakesOver    RCH_MODEL')))
R.append(case('takeover_of_wells'        ,rvg_edit=lambda g:g.replace('RCH_MODEL EVT_MODEL','RCH_MODEL EVT_MODEL WEL_PUMPS')))
R.append(case('stream_not_stream_type'   ,rvg_edit=lambda g:g.replace('RIV_STREAMS AUX SUBBASIN','WEL_PUMPS AUX SUBBASIN')))
R.append(case('unknown_package'          ,rvg_edit=lambda g:g.replace('RIV_STREAMS AUX SUBBASIN','RIV_NOPE AUX SUBBASIN')))
R.append(case('start_date_mismatch'      ,rvg_edit=lambda g:g.replace(':MF6StartDate      1985-10-01',':MF6StartDate      1985-10-02')))
R.append(case('longer_than_model'        ,rvi_edit=lambda t:re.sub(r':Duration\s+\S+',':Duration 800',t)))
R.append(case('period_not_whole_steps'   ,model_edit=bad_period))
R.append(case('two_models_unnamed'       ,model_edit=two_models))
R.append(case('grid_in_wrong_units'      ,grid_edit=feet_grid))
R.append(case('generated_option_left'    ,rvg_edit=lambda g:g+':GridType QUADTREE\n'))
R.append(case('both_link_methods'        ,rvg_edit=lambda g:g+':OverlapWeights\n  2736 1 0.5\n:EndOverlapWeights\n'))
R.append(case('no_link_method'           ,rvg_edit=lambda g:g.replace(':MF6GridFile       ext_cells.geojson CELL_ID\n','')))
print('%d of %d behave as required'%(sum(R),len(R)))
