import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, flopy, sys
d=sys.argv[1]
t=open(d+'gwf.disu').read()
def arr(a,b,skip): return t[t.index('  '+a+'\n'):t.index('  '+b+'\n')].split()[skip:]
iac=np.array(arr('IAC','JA',2),int); ja=np.array(arr('JA','IHC',2),int)-1; cl=np.array(arr('CL12','HWVA',2),float); hw=np.array(arr('HWVA','ANGLDEGX',2),float)
top=np.array(arr('TOP','BOT',4),float); bot=np.array(arr('BOT','AREA',4),float)
idm=np.array(t[t.index('  IDOMAIN\n'):t.index('END GRIDDATA')].split()[4:],int)
n=len(iac); K=10.0; start=np.concatenate([[0],np.cumsum(iac)])
look={}
for i in range(n):
    for k in range(start[i]+1,start[i+1]): look[(i,ja[k])]=k
rows=[];cols=[];vals=[]
for i in range(n):
    for k in range(start[i]+1,start[i+1]):
        j=ja[k]; C=K*(top[i]-bot[i])*hw[k]/(cl[k]+cl[look[(j,i)]])
        rows+=[i,i]; cols+=[i,j]; vals+=[-C,C]
A=sp.csr_matrix((vals,(rows,cols)),shape=(n,n)).tolil(); b=np.zeros(n)
chd=[int(l.split()[0])-1 for l in open(d+'gwf.chd').read().split('BEGIN PERIOD 1')[1].split('END PERIOD')[0].strip().split('\n')]
wel=[int(l.split()[0])-1 for l in open(d+'gwf.wel').read().split('BEGIN PERIOD 1')[1].split('END PERIOD')[0].strip().split('\n')]
act=idm>0
for c in chd: A[c,:]=0; A[c,c]=1; b[c]=80.0
b[wel[0]]+=5000.0
for c in np.where(~act)[0]: A[c,:]=0; A[c,c]=1; b[c]=0
h=spl.spsolve(A.tocsr(),b)
hm=np.array(flopy.utils.HeadFile(d+'gwf.hds').get_alldata()[-1]).ravel()
print('own solution vs MODFLOW: max |dh| = %.3e m over %d active cells; pumped cell own %.4f MODFLOW %.4f'%(np.max(np.abs(h[act]-hm[act])),act.sum(),h[wel[0]],hm[wel[0]]))
np.save(d+'own.npy',h)
