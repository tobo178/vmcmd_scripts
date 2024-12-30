#!/usr/bin/env python
# 
# Written by Satoshi Ono
#
# gen_cano prototype
#
# Usage: gen_cano.py -c md9 -o gen_cano.md9 -T 300
#
#

import argparse
import os
from pathlib import Path
import numpy as np
from numpy.polynomial import polynomial as P
import kkkit
import rangeinfo
#import matplotlib.pyplot as plt

# 
g = 1.0/(1.3806488e-23 * 6.02214129e23 * 1.0e-3)
rgas = 8.314462145 / 1000.0

parser = argparse.ArgumentParser(description='Test')
parser.add_argument('-c', '--current', dest='curr_md', nargs='?',
                    required=True, type=str,
                    help='Directory for current MD')
parser.add_argument('-o', '--out', dest='gen_cano', nargs='?',
                    type=str, help='Directory for gen_cano')
parser.add_argument('-T', '--temp', dest='temp', nargs='?',
                    type=str, help='Temperature for gen_cano')
parser.add_argument('-q', '--quiet', action='store_true',default=False)
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

RT= rgas * float(args.temp)

# First, check if the output directory exist.
if Path(args.gen_cano).exists():
    print("ERROR: {} exist. Stopping...".format(args.gen_cano))
    exit(1)

######################################################################
# Read current V-McMD parameters from temp_s
ri = rangeinfo.RangeInfo('{}/temp_s'.format(args.curr_md))
ri.read_ttp_inp()
nvst = len(ri.vs_order)
vinterval = ri.vinterval

######################################################################
# fit_pmc
ene=[]
pdd=[]

for i in range(nvst):
    with open('{0}/e{1:d}.pdd'.format(args.curr_md, i+1),'r') as f:
        ene_min=(ri.vs_range[i])[0]
        ene_max=(ri.vs_range[i])[1]
        ene_in, pdd_in = np.loadtxt(f, unpack=True)
        ene_min_ind=np.where(ene_in==ene_min)[0][0]
        ene_max_ind=np.where(ene_in==ene_max)[0][0]+1
        ene.append(list(ene_in[ene_min_ind:ene_max_ind]))
        pdd.append(list(pdd_in[ene_min_ind:ene_max_ind]))

fl=[]
fdl=[]
for i in range(nvst):
    ndeg=ri.vs_order[i][0]
    c, stat = P.polyfit(ene[i],pdd[i],ndeg-1,full=True)
    fd = P.polyder(c)
    fdl.append(list(fd))

######################################################################
# fit_mix
buf = ""
buf += '{0:d}\n'.format(nvst)
buf += '{0:d}\n'.format(vinterval)

for i in range(nvst):
    buf += '{0:.1f} {1:.1f}\n'.format(ri.vs_range[i][0], ri.vs_range[i][1])
    buf += '{0:.2f} {1:.2f}\n'.format(ri.vs_range[i][2], ri.vs_range[i][3])

for j in range(nvst):
    ndeg=ri.vs_order[j][0]
    coeff1=np.array(ri.vs_params[j],dtype=float)
    coeff1=coeff1[:-2]
    coeff2=np.array(fdl[j],dtype=float)
    iord=max(ndeg,len(fdl[j])-1)
    buf += '{0:d}\n'.format(iord)
    coeff3=P.polyadd(coeff1,coeff2)
    for i in range(len(coeff3)):
        buf += '{0:22.15E}\n'.format(coeff3[i])
    buf += '{0:22.15E}\n'.format(P.polyval(ri.vs_range[j][0],fdl[j]))
    buf += '{0:22.15E}\n'.format(P.polyval(ri.vs_range[j][1],fdl[j]))
buf += '{0:.2f}\n'.format(ri.temp)

p = Path(args.gen_cano)
p.mkdir()
with open(p.joinpath('temp_s'), 'w') as f:
    f.write(buf)
if args.verbose:
    print(buf, end="")

######################################################################
# Read temp_s for gen_cano
ri2 = rangeinfo.RangeInfo('{}/temp_s'.format(args.gen_cano))
ri2.read_ttp_inp()
nvst2 = len(ri2.vs_order)

ene = np.array([])
pdd = np.array([])

offset=1000.0
for i in range(nvst2):
    ene_min = ri2.vs_range[i][0]
    ene_max = ri2.vs_range[i][1]
    # Expand the energy range
    if i == 0:
        ene_min -= offset
    elif i == nvst2-1:
        ene_max += offset
    x = np.arange(ene_min, ene_max+1,10.0)
    coeff=np.array(ri2.vs_params[i],dtype=float)
    coeff=coeff[:-2]
    y=P.polyval(x,coeff)
    ene=np.concatenate((ene,x))
    pdd=np.concatenate((pdd,y))
 
c, stat = P.polyfit(ene,pdd,ri2.vs_order[0][0],full=True)

ene = np.arange(ri2.vs_range[0][0]-offset, ri2.vs_range[nvst2-1][1]+offset, 20.0)
cano = -ene/RT + P.polyval(ene, P.polyint(c))
cano = cano - np.max(cano)
np.savetxt('{0}/P_E_T_{1}.dat'.format(args.gen_cano, args.temp),
           np.transpose(np.array([ene,cano])), fmt='%.10e')
