#!/usr/bin/env python
#
# Written by Satoshi Ono
#
# make_md0.py -n md0_No -N name -d job_id
#

import argparse
import os
from pathlib import Path
import numpy as np
from vmcmd_config import *

conf=vmcmd_config()
md_num=conf.md_num
temp_s = [0.08,0.1,0.12,0.15,0.17,0.19,0.21,0.24,0.27,0.3,0.34,0.37,0.4,0.431]
#np.random.seed(seed=123456)
idx=np.random.choice(1048565, md_num,replace=False)
idx=np.sort(idx)
idx2=np.random.choice(1048565, md_num,replace=False)

parser = argparse.ArgumentParser(description='run_md0')
parser.add_argument('-N', '--name', dest='name', nargs='?',
                    required=True, type=str, help='Name for queue')
parser.add_argument('-n', '--md0_No', dest='no', nargs='?',
                    required=True, type=int, help='md0 number')
parser.add_argument('-d', '--dependency', dest='job_id', nargs='?',
                    required=True, type=int, help='Job ID')
args = parser.parse_args()

p = Path(".")
# job submitter
#jname = p.cwd().parts[-1]+'_'+args.next_md
#jname = p.cwd().parts[-2]
jname=args.name

with open (p.joinpath('temp_s1'), 'w') as f:
    f.write('''1
10000000
-2000000.0  2000000.0
0.0  0.0
1
{TEMPS}
0.0
0.0
0.0
300.0
'''.format(TEMPS=temp_s[args.no-1]))

if args.no == 1:
    for i in range(md_num):
        os.mkdir(p.joinpath('{N}'.format(N=i+1)))
        with open (p.joinpath('{N}/md.mdp'.format(N=i+1)), 'w') as f:
            f.write(''';
gen_seed            = {RANDOM}
gen_vel             = yes
gen_temp            = 300
annealing           = single single
annealing-npoints   = 2 2
annealing-temp      = 0 300 0 300
annealing-time      = 0 1 0 1
integrator          = md
dt                  = 0.002
nsteps              = 500000
nstxout             = 0
nstvout             = 0
nstfout             = 0
nstlog              = 10000
nstenergy           = 10000
nstcalcenergy       = 1
nstxout-compressed  = 10000
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
do-vmcmd            = yes
vmcmd-param-file    = ../temp_s1
vmcmd-start-file    = start.vert
'''.format(RANDOM=idx[i]))
        with open (p.joinpath('{N}/start.vert'.format(N=i+1)), 'w') as f:
            f.write('''0
{IDX2}
'''.format(IDX2=idx2[i]))

        with open (p.joinpath('cano.sh'), 'w') as f:
            f.write('''#!/bin/bash -l
#SBATCH --output=%a/md.stdout
#SBATCH --array=1-{MD_NUM}
#SBATCH --job-name={JNAME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
source $HOME/gmx2024.4_vmcmd_cpu/bin/GMXRC.bash
export OMP_NUM_THREADS=1
#export GMX_DISABLE_GPU_DETECTION=1
cd $SLURM_ARRAY_TASK_ID
gmx grompp -c ../../../em1/npt.gro -p ../../../topol.top -f md.mdp -o md.tpr -po md.out.mdp
gmx mdrun -nov -deffnm md -ntomp 1
'''.format(MD_NUM=md_num,JNAME=jname+"_"+str(args.no)))

#
else:
    for i in range(md_num):
        os.mkdir(p.joinpath('{N}'.format(N=i+1)))
        with open (p.joinpath('{N}/start.vert'.format(N=i+1)), 'w') as f:
            f.write('''0
{IDX2}
'''.format(IDX2=idx2[i]))
    with open (p.joinpath('md.mdp'), 'w') as f:
            f.write(''';
gen_vel             = no
continuation        = yes
integrator          = md
dt                  = 0.002
nsteps              = 100000
nstxout             = 0
nstvout             = 0
nstfout             = 0
nstlog              = 10000
nstenergy           = 10000
nstcalcenergy       = 1
nstxout-compressed  = 10000
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
do-vmcmd            = yes
vmcmd-param-file    = ../temp_s1
vmcmd-start-file    = start.vert
''')
    with open (p.joinpath('cano.sh'), 'w') as f:
        f.write('''#!/bin/bash -l
#SBATCH --output=%a/md.stdout
#SBATCH --array=1-{MD_NUM}
#SBATCH --job-name={JNAME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --dependency=aftercorr:{JNAME_PRE}
source $HOME/gmx2024.4_vmcmd_cpu/bin/GMXRC.bash
export OMP_NUM_THREADS=1
#export GMX_DISABLE_GPU_DETECTION=1
cd $SLURM_ARRAY_TASK_ID
gmx grompp -c ../../{PRE}/$SLURM_ARRAY_TASK_ID/md.tpr -t ../../{PRE}/$SLURM_ARRAY_TASK_ID/md.cpt -p ../../../topol.top -f ../md.mdp -o md.tpr -po md.out.mdp
gmx mdrun -nov -deffnm md -ntomp 1
'''.format(MD_NUM=md_num,PRE=args.no-1,JNAME=jname+"_"+str(args.no),JNAME_PRE=args.job_id))
#'''.format(MD_NUM=md_num,PRE=args.no-1,JNAME=jname+"_"+str(args.no),JNAME_PRE=jname+"_"+str(args.no-1)))
