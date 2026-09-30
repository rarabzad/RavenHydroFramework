import os,sys; sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); import testconfig as tc
import ctypes,os,sys,time
d=sys.argv[1]; tc.need_mf6(); lib=ctypes.CDLL(tc.ENV['RAVEN_MF6_LIB']); os.chdir(d)
assert lib.initialize()==0,'initialize failed (see mfsim.lst)'
cur=ctypes.c_double(); end=ctypes.c_double(); lib.get_current_time(ctypes.byref(cur)); lib.get_end_time(ctypes.byref(end))
t0=time.time(); n=0
while cur.value<end.value:
    if lib.update()!=0: print('update failed at',cur.value); break
    lib.get_current_time(ctypes.byref(cur)); n+=1
lib.finalize(); print('steps',n,'end time',cur.value,'wall %.1f s'%(time.time()-t0))
