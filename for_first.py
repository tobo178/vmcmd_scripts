#!/usr/bin/env python
#
# Written by Satoshi Ono
#
# for_first
# Calculate Probability Density Distribution (pdd) from
# series of temparatures.
# Usage: for_first.py -i ini* -n md1
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
skipsteps=10000
ndeg = 8
#nvst=8
vinterval=20000
g = 1.0/(1.3806488e-23 * 6.02214129e23 * 1.0e-3)
wd=os.path.dirname(os.path.abspath(__file__))

def calc_range(emin,emax,bin_width,r):
    """Calculate potential energy range from these dictionaries."""
    r7= {0: (0.00, 0.10, 1.0, 1.0),
         1: (0.05, 0.20, 1.0, 1.0),
         2: (0.10, 0.30, 1.0, 1.0),
         3: (0.20, 0.45, 1.0, 1.0),
         4: (0.30, 0.60, 1.0, 1.0),
         5: (0.45, 0.80, 1.0, 1.0),
         6: (0.60, 1.00, 1.0, 1.0)}
    r8= {0: (0.00, 0.20, 1.0, 1.0),
         1: (0.10, 0.30, 1.0, 1.0),
         2: (0.20, 0.40, 1.0, 1.0),
         3: (0.30, 0.50, 1.0, 1.0),
         4: (0.40, 0.60, 1.0, 1.0),
         5: (0.50, 0.70, 1.0, 1.0),
         6: (0.60, 0.80, 1.0, 1.0),
         7: (0.70, 1.00, 1.0, 1.0)}
    if r == 7:
        rr = r7
    else:
        rr = r8

    e_width=emax-emin
    vs_range={}
    for i in range(len(rr)):
        ratio_min = float(rr[i][0])
        ratio_max = float(rr[i][1])
        ene_min = int((emin + e_width * ratio_min)/bin_width)*bin_width
        ene_max = int((emin + e_width * ratio_max)/bin_width)*bin_width
        vs_range[i]=(ene_min, ene_max, rr[i][2], rr[i][3])
    return vs_range

parser = argparse.ArgumentParser(description='Test')
parser.add_argument('-i', '--initial', dest='init_md', nargs='+',
                    required=True, type=str,
                    help='Directories for initial MDs')
parser.add_argument('-n', '--next', dest='next_md', nargs='?',
                    type=str, help='Directory for Next MD')
parser.add_argument('-H', '--histogram', action='store_true', default=False,
                    help='Calculate Histogram(s)')
parser.add_argument('-s', '--submit', action='store_true', default=False,
                    help='Submit jobs')
parser.add_argument('-N', '--name', dest='job_name', type=str,
                    help='Job name', default='SOL_mdn')
parser.add_argument('-r', '--range', dest='range', nargs='?', type=int,
                    help='Range')
parser.add_argument('-q', '--quiet', action='store_true',default=False)
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

# First, check if the output directory exist.
if Path(args.next_md).exists():
    print("ERROR: {} exist. Stopping...".format(args.next_md))
    exit(1)

######################################################################
# distrib
np.seterr(divide='ignore')
for dir in tqdm(args.init_md, disable=args.quiet):
    p = Path(dir)
    if ( (p / 'e1.pdd').is_file() and not args.histogram):
        if not args.quiet:
            print('Warning: {}/e1.pdd exist. Skipping...'.format(str(p)))
        continue
    ene = []
    for s in sorted(p.glob("**/md.vmcmd.ene")):
        tmp=np.fromfile(str(s),dtype='float64')
        tmp_list=tmp[skipsteps:].tolist()
        ene.extend(tmp_list)
    pot=np.array(ene)
    ene_min, ene_max= np.amin(pot[:]), np.amax(pot[:])
    
    imin = int((ene_min-width)/width)
    imax = int((ene_max+width)/width)
    bin_index= np.arange(imin,imax+1)*width
    y,ind = np.histogram(pot+width/2,bins=bin_index)

    y2=y/pot.size
    y2=np.log(y2)
    np.savetxt('{}/e1.pdd'.format(dir),
               np.transpose(np.array([ind[np.where(y2!=-np.inf)],
                                      y2[np.where(y2!=-np.inf)]])))

######################################################################
# derv_den_Pc
dir_num=0
for dir in args.init_md:
    dir_num += 1
    ri=rangeinfo.RangeInfo('{}/temp_s1'.format(dir))
    ri.read_ttp_inp()
    temp=g/float(ri.vs_params[0][0])
    rt=temp*8.314462145/1000.0
    ene, pdd = np.loadtxt('{}/e1.pdd'.format(dir), unpack=True)
    pddmax=np.amax(pdd)
    pdd2=pdd[np.where(pdd>=pddmax-1)]
    ene2=ene[np.where(pdd>=pddmax-1)]
    den=ene2/rt + pdd2
    # Assume ene2.diff's are same, equals width
    grad=np.gradient(den,ene2[1]-ene2[0],edge_order=2)
    np.savetxt('{}/dden.dat'.format(dir),
               np.transpose(np.array([ene2,grad])),fmt='%.10e')

######################################################################
# fit_dden
ene=np.array([])
dden=np.array([])
for dir in args.init_md:
    inp = np.loadtxt('{}/dden.dat'.format(dir))
    ene=np.append(ene,inp[:,0])
    dden=np.append(dden,inp[:,1])

c=np.polyfit(ene,dden,ndeg)
f=np.poly1d(c)
fd=np.polyder(f)

ene_min = np.round(np.min(ene)/width)*width
ene_max = np.round(np.max(ene)/width)*width

rr = calc_range(ene_min, ene_max, width, args.range)

buf = ""
buf += '{0:d}\n'.format(len(rr))
buf += '{0:d}\n'.format(vinterval)

for i in range(len(rr)):
    buf += '{0:.1f} {1:.1f}\n'.format(rr[i][0], rr[i][1])
    buf += '{0:.2f} {1:.2f}\n'.format(rr[i][2], rr[i][3])
for i in range(len(rr)):
    buf += '{0:d}\n'.format(ndeg)
    for i in range(np.size(c)-1,-1,-1):
        buf +='{0:22.15E}\n'.format(c[i])
    buf += '{0:22.15E}\n'.format(fd(ene_min))
    buf += '{0:22.15E}\n'.format(fd(ene_max))
buf += '300.0\n'

p=Path(args.next_md)
p.mkdir()
with open(p.joinpath( 'temp_s'), 'w') as f:
    f.write(buf)
if args.verbose:
    print(buf,end="")

# for plot_ene.gp
ll=[]
for i in range(len(rr)):
    ll.append(rr[i][0])
    ll.append(rr[i][1])
ll_uniq=list(set(ll))
ll_uniq.sort()

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
prev=$(( (RANDOM % {dir_num}) + 1 ))
iprev=`python $BIN/get_last_state.py -f ../../md0/$prev/$SLURM_ARRAY_TASK_ID/md.edr -t ../temp_s`
echo $iprev > start.vert
echo $(((RANDOM<<15)|RANDOM)) >> start.vert
export OMP_NUM_THREADS=1
gmx grompp -c ../../md0/$prev/$SLURM_ARRAY_TASK_ID/md.tpr \
  -t ../../md0/$prev/$SLURM_ARRAY_TASK_ID/md.cpt \
  -p ../../topol.top \
  -f ../md.mdp -o md.tpr -po md.out.mdp 
gmx mdrun -nov -deffnm md -ntomp 1
sleep 5
if [ -f md.vmcmd.ene ] ; then
  python $BIN/make_hist.py
#  rm md.vmcmd.ene
fi
'''.format(md_num=md_num, wd=wd, dir_num=dir_num, gmx=conf.gmx))

# jname = p.cwd().parent.parts[-1]+'_'+args.next_md
# with open(p.joinpath('01_mkdir_and_submit.sh'), 'w') as f:
#     f.write('''#!/bin/bash
# export BIN={0}
# #mkdir out
# #for ii in `seq 1 {1}`
# #do
# #  mkdir $ii
# #done
# #ln -s $BIN/get_last_state.py .
# sbatch -J "{2}" ./multi.sh
# exit'''.format(wd, md_num, jname))
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
