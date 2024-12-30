#!/usr/bin/env python
# 
# Written by Satoshi Ono
#
# make_hist
#
# Usage: make_hist
#
import numpy as np
import kkkit
import rangeinfo

width=100.0
skipsteps=10000
g = 1.0/(1.3806488e-23 * 6.02214129e23 * 1.0e-3)

######################################################################
# Read current V-McMD parameters from temp_s
ri = rangeinfo.RangeInfo('../temp_s')
ri.read_ttp_inp()
nvst = len(ri.vs_order)
vinterval = ri.vinterval

######################################################################
# v_distrib
# Pass 1: Find ene_min and ene_max
np.seterr(divide='ignore')
ene_min, ene_max=1.e99, -1.e99
ene=np.fromfile('md.vmcmd.ene',dtype='float64')
ene_min=min(min(ene[skipsteps:]),ene_min)
ene_max=max(max(ene[skipsteps:]),ene_max)

# Histogram properties
imin = int((ene_min-width)/width)
imax = int((ene_max+width)/width)
bin_index= np.arange(imin,imax+1)*width
hist=np.zeros((nvst,(imax-imin)))
ene_count=np.zeros(nvst,dtype=int)

# Pass 2: Make histograms
vstep=[]
vstate=[]
vstatenew=[]
tmp=np.loadtxt("md.vmcmd.log", dtype=int)
vstep.extend(tmp[:,0].tolist())
vstate.extend(tmp[:,1].tolist())
vlast=vstep[-1]
for i in range(0,vlast,vinterval):
    for j in range(0,vinterval):
        vstatenew.append(vstate[int(i/vinterval)])
vstatenew.append(vstate[int(vlast/vinterval)])
ene=ene[skipsteps:]
vstate=np.array(vstatenew[skipsteps:])

for i in range(nvst):
    ene3=ene[np.where(vstate==i)]
    y, ind = np.histogram(ene3+width/2,bins=bin_index)
    hist[i] += y
    ene_count[i] += ene3.size

ind2 =np.reshape(ind, (1, ind.shape[0])) # 2D array
# exclude last bin
np.savetxt("md.vmcmd.hst",np.concatenate([ind2[:,:-1],hist]).T,fmt='%10d')
