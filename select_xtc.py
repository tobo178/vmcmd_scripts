#!/usr/bin/env python
#
# Written by Satoshi Ono
#

import argparse
import numpy as np
import mdtraj as md

parser = argparse.ArgumentParser(description='select')
parser.add_argument('-f', '--traj', dest='traj', nargs='?',
                    required=True, help='test')
parser.add_argument('-s', '--top', dest='top', nargs='?',
                    required=True)
parser.add_argument('-n', '--num', dest='num', nargs='?',type=int,
                    required=True)
parser.add_argument('-o', '--out', dest='out', nargs='?',
                    required=True)
parser.add_argument('-q', '--quiet', action='store_true',default=False)
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

np.random.seed(seed=123456)

sel=args.num
newtraj=[]

t=md.load(args.traj,top=args.top)

if t.n_frames < sel:
    printf("Traj has < {} flames. exitting..".format(sel))
    exit(1)

idx=np.random.choice(t.n_frames,sel,replace=False)
idx=np.sort(idx)

for i in idx:
    newtraj.append(t[i])
newtraj=md.join(newtraj)
newtraj.save_xtc(args.out)
