#!/usr/bin/env python
# 
# Written by Satoshi Ono
#
# for_next 
#
# Usage: for_next.py -c md1 -n md2
#
# use for_first.py to making the first md1 from ini*
#

import argparse
import os
from pathlib import Path
import numpy as np
import kkkit
import rangeinfo
from vmcmd_config import *

# jupyter or CLI
try:
    # noinspection PyUnresolvedReferences
    if get_ipython().__class__.__name__ == 'ZMQInteractiveShell':
        from tqdm import tnrange as trange
        from tqdm import tqdm_notebook as tqdm
    else:
        raise RuntimeError
except (NameError, RuntimeError):
    from tqdm import trange
    from tqdm import tqdm
#
conf=vmcmd_config()

md_num=conf.md_num
width=100.0
skipsteps=10001
# ndeg = 8
# nvst=8
#interval=20000
g = 1.0/(1.3806488e-23 * 6.02214129e23 * 1.0e-3)
wd=os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description='Test')
parser.add_argument('-c', '--current', dest='curr_md', nargs='?',
                    required=True, type=str,
                    help='Directory for current MD')
parser.add_argument('-n', '--next', dest='next_md', nargs='?',
                    type=str, help='Directory for Next MD')
parser.add_argument('-H', '--histogram', action='store_true', default=False,
                    help='Calculate Histogram(s)')
parser.add_argument('-s', '--submit', action='store_true', default=False,
                    help='Submit jobs')
parser.add_argument('-N', '--name', dest='job_name', type=str,
                    help='Job name')
parser.add_argument('-q', '--quiet', action='store_true',default=False)
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

# First, check if the output directory exist.
if Path(args.next_md).exists():
    print("ERROR: {} exist. Stopping...".format(args.next_md))
    exit(1)

######################################################################
# Read current V-McMD parameters from temp_s
ri = rangeinfo.RangeInfo('{}/temp_s'.format(args.curr_md))
ri.read_ttp_inp()
nvst = len(ri.vs_order)
vinterval = ri.vinterval

######################################################################
# v_distrib
# Pass 1: Find ene_min and ene_max
np.seterr(divide='ignore')
ene_min, ene_max=1.e99, -1.e99
p = Path(args.curr_md)
if list(p.glob('**/*.pdd'))  and not args.histogram:
    print("pdd file exist. Skipping")
else:
    for s in tqdm(sorted(p.glob("**/md.vmcmd.ene")),desc='Ene(Pass1)'):
        tmp=np.fromfile(str(s),dtype='float64')
        ene_min=min(min(tmp[skipsteps:]),ene_min)
        ene_max=max(max(tmp[skipsteps:]),ene_max)

    # Histogram properties
    imin = int((ene_min-width)/width)
    imax = int((ene_max+width)/width)
    bin_index= np.arange(imin,imax+1)*width
    hist=np.zeros((nvst,(imax-imin)))
    ene_count=np.zeros(nvst,dtype=int)

    # Pass 2: Make histograms
    for s in tqdm(sorted(p.glob("**/md.vmcmd.ene")),desc='Ene(Pass2)'):
        vstep=[]
        vstate=[]
        vstatenew=[]
        tmp=np.loadtxt("{}/md.vmcmd.log".format(s.parent), dtype=int)
        vstep.extend(tmp[:,0].tolist())
        vstate.extend(tmp[:,1].tolist())
        ene=np.fromfile(str(s))
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

    for i in range(nvst):
        y2=hist[i]/ene_count.sum()
        y2=np.log(y2)
        np.savetxt('{0}/e{1:d}.pdd'.format(args.curr_md,i+1), 
                   np.transpose(np.array([ind[np.where(y2!=-np.inf)],
                                          y2[np.where(y2!=-np.inf)]])))

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
    c=np.polyfit(ene[i],pdd[i],ndeg-1)
    f=np.poly1d(c)
    fd=np.polyder(f)
    fl.append(list(f))
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
    coeff1=np.zeros((ndeg+1))
    coeff2=np.zeros((ndeg+1))
    for i in range(ndeg+1):
        coeff1[i]=ri.vs_params[j][i]
    for i in range(len(fdl[j])):
        coeff2[i]=list(reversed(fdl[j]))[i]
    iord=max(ndeg,len(fdl[j])-1)
    buf += '{0:d}\n'.format(iord)
    for i in range(iord+1):
        buf += '{0:22.15E}\n'.format(coeff1[i]+coeff2[i])
    buf += '{0:22.15E}\n'.format(np.poly1d(fdl[j])(ri.vs_range[j][0]))
    buf += '{0:22.15E}\n'.format(np.poly1d(fdl[j])(ri.vs_range[j][1]))
buf += '{0:.2f}\n'.format(ri.temp)

p=Path(args.next_md)
p.mkdir()
with open(p.joinpath( 'temp_s'), 'w') as f:
    f.write(buf)
if args.verbose:
    print(buf,end="")

######################################################################
# files for next_md

# for plot_ene.gp
ll=[]
for i in range(len(ri.vs_range)):
    ll.append(ri.vs_range[i][0])
    ll.append(ri.vs_range[i][1])
ll_uniq=list(set(ll))
ll_uniq.sort()

# Slurm
with open(p.joinpath('multi.sh'), 'w') as f:
    f.write('''#!/bin/bash
#SBATCH --output=%a/md.stdout
#SBATCH --array=1-{md_num}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
source {gmx}
export BIN={wd}
export PATH=$HOME/miniforge3/bin:$PATH
export PYTHONPATH=$BIN:$PYTHONPATH
mkdir $SLURM_ARRAY_TASK_ID
cd $SLURM_ARRAY_TASK_ID
iprev=`tail -1 ../../{PREV}/$SLURM_ARRAY_TASK_ID/md.vmcmd.log |awk '{{print $2}}'`
echo $iprev > start.vert
echo $(((RANDOM<<15)|RANDOM)) >> start.vert
export OMP_NUM_THREADS=1
gmx grompp -c ../../{PREV}/$SLURM_ARRAY_TASK_ID/md.tpr \
  -t ../../{PREV}/$SLURM_ARRAY_TASK_ID/md.cpt \
  -p ../../topol.top \
  -f ../md.mdp -o md.tpr -po md.out.mdp 
gmx mdrun -nov -deffnm md -ntomp 1
sleep 5
if [ -f md.vmcmd.ene ] ; then
  python $BIN/make_hist.py
#  rm md.vmcmd.ene
fi
'''.format(md_num=md_num, PREV=args.curr_md, wd=wd, gmx=conf.gmx))

# # job submitter
# jname = p.cwd().parts[-1]+'_'+args.next_md
# with open(p.joinpath('01_mkdir_and_submit.sh'), 'w') as f:
#     f.write('''#!/bin/bash
# # mkdir out
# # for ii in `seq 1 {0}`
# # do
# #   mkdir $ii
# # done
# sbatch -J "{1}" ./multi.sh
# exit'''.format(md_num, jname))
# os.chmod(p.joinpath('01_mkdir_and_submit.sh'), 0o755)

with open(p.joinpath('md.mdp'), 'w') as f:
    f.write('''; mdp file for TTP-V-McMD
gen_vel             = no
continuation        = yes
do-vmcmd            = yes
vmcmd-param-file    = ../temp_s
vmcmd-start-file    = start.vert
integrator          = md
dt                  = 0.002
nsteps              = {NSTEPS}
nstxout             = 0
nstvout             = 0
nstfout             = 0
nstlog              = 5000
nstenergy           = {NSTENERGY}
nstcalcenergy       = 1
nstxout-compressed  = {NSTXOUT}
compressed-x-grps   = Protein
; energygrps          = Protein non-Protein
nstlist             = 10
ns_type             = grid
coulombtype         = PME 
rvdw                = 1.0 
rlist               = 1.0
rcoulomb            = 1.0 
cutoff-scheme       = Verlet
pbc                 = xyz
tcoupl              = v-rescale
tc-grps             = Protein non-Protein
tau_t               = 0.1 0.1
ref_t               = 300 300
constraints         = h-bonds
;
'''.format(NSTEPS=1000000, NSTENERGY=10000, NSTXOUT=10000))

# plot_log.gp
with open(p.joinpath('plot_log.gp'), 'w') as f:
    f.write('''#!/usr/bin/gnuplot
reset
start  = 1
end    = {md_num}
offset = 12
set yrange[-1:8]
do for [i = start : end : offset] {{
   plot for [j=0:offset-1] sprintf("%d/md.vmcmd.log",i+j) ev 1 u 2 t sprintf("%dth",i+j) w l
 pause -1
}}
'''.format(md_num=md_num))
os.chmod(p.joinpath('plot_log.gp'), 0o755)
