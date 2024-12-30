#!/usr/bin/env python
#
# Written by Satoshi Ono
#

import numpy as np
import mdtraj as md

np.random.seed(seed=123456)

sel=5000
newtraj=[]

t=md.load('300K_fit.xtc',top='md_Ref.gro')

if t.n_frames < sel:
    printf("Traj has < {} flames. exitting..".format(sel))
    exit(1)

idx=np.random.choice(t.n_frames,sel,replace=False)
idx=np.sort(idx)

for i in idx:
    newtraj.append(t[i])
newtraj=md.join(newtraj)
newtraj.save_xtc("300K_fit_5000.xtc")
